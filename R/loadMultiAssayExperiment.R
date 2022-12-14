#' Load a MultiAssayExperiment
#'
#' Load a dataset as a \linkS4class{MultiAssayExperiment}, based on the metadata generated by the corresponding \code{\link{stageObject}} method.
#'
#' @param ds.info Named list containing the metadata for this object.
#' @param project Any argument accepted by the acquisition functions, see \code{?\link{acquireFile}}. 
#' By default, this should be a string containing the path to a staging directory.
#' @param experiments Character or integer vector specifying the subset of experiments to load.
#' If \code{NULL}, all experiments are loaded.
#' @param BPPARAM A BiocParallelParam object indicating how loading should be parallelized across multiple experiments.
#' If \code{NULL}, loading is done serially.
#' @param include.nested Logical scalar indicating whether to include nested DataFrames in the \code{colData} of the output.
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
#' @import alabaster.base alabaster.se
#' @importFrom MultiAssayExperiment MultiAssayExperiment
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
