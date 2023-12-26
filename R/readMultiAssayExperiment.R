#' Read a MultiAssayExperiment from disk
#'
#' Read a \linkS4class{MultiAssayExperiment} from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{stageObject}} method for \linkS4class{MultiAssayExperiment} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to internal \code{\link{altReadObject}} calls.
#' 
#' @return A \linkS4class{MultiAssayExperiment} object.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' # Mocking up an MAE
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#' se <- SummarizedExperiment(list(counts=mat))
#' 
#' library(MultiAssayExperiment)
#' mae <- MultiAssayExperiment(list(gene=se))
#'
#' # Staging it:
#' tmp <- tempfile()
#' dir.create(tmp)
#' info <- stageObject(mae, tmp, "dataset")
#'
#' # Loading it back in:
#' loadMultiAssayExperiment(info, tmp)
#' 
#' @export
#' @aliases loadMultiAssayExperiment
#' @import alabaster.base alabaster.se
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom jsonlite fromJSON
#' @importFrom rhdf5 H5Fopen H5Fclose H5Gopen H5Gclose
readMultiAssayExperiment <- function(path, metadata, ...) {
    # Constructing the experiments.
    edir <- file.path(path, "experiments")
    if (file.exists(edir)) {
        exp.names <- fromJSON(file.path(edir, "names.json"))
        experiments <- list()

        for (e in seq_along(exp.names)) {
            ex <- readObject(file.path(edir, e-1L), ...)
            if (is.null(colnames(ex))) {
                warning("generating non-NULL column names for experiment ", e)
                colnames(ex) <- seq_len(ncol(ex))
            } else if (anyDuplicated(colnames(ex))) {
                warning("generating unique column names for experiment ", e)
                colnames(ex) <- make.unique(colnames(ex))
            }
            experiments[[e]] <- ex
        }

        names(experiments) <- exp.names
    }

    # Obtaining the sample data.
    sd <- readObject(file.path(path, "sample_data"), ...)
    if (is.null(rownames(sd))) {
        warning("generating non-NULL sample names for the sample data")
        rownames(sd) <- seq_len(nrow(sd))
    } else if (anyDuplicated(rownames(sd))) {
        warning("generating unique sample names for the sample data")
        rownames(sd) <- make.unique(rownames(sd))
    }

    # Constructing the sample mapping.
    fhandle <- H5Fopen(file.path(path, "sample_map.h5"), "H5F_ACC_RDONLY")
    on.exit(H5Fclose(fhandle))
    ghandle <- H5Gopen(fhandle, "multi_sample_dataset")
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    all.primary <- list()
    all.assay <- list()
    all.colname <- list()

    for (e in seq_along(experiments)) {
        curex <- experiments[[e]]
        ix <- h5_read_vector(ghandle, as.character(e-1L))
        all.assay[[e]] <- rep(names(experiments)[e], ncol(curex))
        all.colname[[e]] <- colnames(curex)
        all.primary[[e]] <- rownames(sd)[ix + 1L]
    }

    sm <- DataFrame(
        primary=unlist(all.primary), 
        assay=factor(unlist(all.assay)), 
        colname=unlist(all.colname)
    )

    mae <- MultiAssayExperiment(experiments, colData=sd, sampleMap=sm) 
    readMetadata(mae, metadata.path=file.path(path, "other_data"), mcols.path=NULL, ...)
}

##################################
######### OLD STUFF HERE #########
##################################

#' @export
loadMultiAssayExperiment <- function(ds.info, project, experiments=NULL, BPPARAM=NULL, include.nested=TRUE) {
    # Choosing the experiments to load.
    all.experiments <- ds.info$dataset$experiments
    keep <- .choose_experiments(experiments, all.experiments)
    all.experiments <- all.experiments[keep]

    if (is.null(BPPARAM)) {
        all.exps <- lapply(all.experiments, .load_experiment, project=project)
    } else {
        all.exps <- BiocParallel::bplapply(all.experiments, .load_experiment, project=project, BPPARAM=BPPARAM)

        # If the experiment loading caused new packages to be loaded in the
        # workers, we want to ensure those class-defining packages are
        # available; otherwise the MAE constructor will complain. It seems
        # sufficient to run some S4 methods that trigger package loading.
        lapply(all.exps, function(x) colnames(x[[1]])) 
    }
    all.exps <- do.call(c, all.exps)

    # Getting the sample mapping.
    map.info <- acquireMetadata(project, ds.info$dataset$sample_mapping$resource$path)
    mapping.raw <- .loadObject(map.info, project)
    mapping <- DataFrame(
        assay = factor(mapping.raw$experiment, names(all.exps)), # https://github.com/waldronlab/MultiAssayExperiment/issues/290#issuecomment-879206815
        primary = mapping.raw$sample,
        colname = mapping.raw$column
    )

    # Getting the subject data; this had better be a DataFrame.
    subject.info <- acquireMetadata(project, ds.info$dataset$sample_data$resource$path)
    coldata <- .loadObject(subject.info, project, include.nested=include.nested) 

    mae <- MultiAssayExperiment(all.exps, sampleMap=mapping, colData=coldata)
    .restoreMetadata(mae, mcol.data=NULL, meta.data=ds.info$dataset$other_data, project)
}

#' @import alabaster.base
.load_experiment <- function(exp.details, project) {
    meta <- acquireMetadata(project, exp.details$resource$path)
    output <- list(.loadObject(meta, project=project))
    names(output) <- exp.details$name
    output
}

.choose_experiments <- function(experiments, all.experiments) {
    if (is.null(experiments)) {
        experiments <- seq_along(all.experiments)
    } else if (is.character(experiments)) {
        all.names <- vapply(all.experiments, function(x) x$name, "")
        m <- match(experiments, all.names)
        if (any(lost <- is.na(m))) {
            stop("cannot find '", experiments[lost][1], "' in the available experiments")
        }
        experiments <- m
    } else {
        if (any(experiments <= 0 | experiments > length(all.experiments))) {
            stop("'experiments' must be positive and no greater than the total number of experiments (", length(all.experiments), ")")
        }
    }
    experiments
}
