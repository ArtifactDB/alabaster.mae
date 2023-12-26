#' Save a MultiAssayExperiment
#'
#' Save a \linkS4class{MultiAssayExperiment} to its on-disk representation.
#' 
#' @param x A \linkS4class{MultiAssayExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::saveObject
#'
#' @return The contents of \code{x} are saved into a \code{path}, and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Mocking up an MAE
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#' 
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(list(counts=mat))
#' 
#' library(MultiAssayExperiment)
#' mae <- MultiAssayExperiment(list(gene=se))
#'
#' # Saving it:
#' tmp <- tempfile()
#' saveObject(mae, tmp)
#' 
#' @export
#' @rdname saveMultiAssayExperiment
#' @aliases stageObject,MultiAssayExperiment-method
#' @import methods alabaster.base alabaster.se
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment colData sampleMap
#' @importFrom S4Vectors DataFrame
#' @importFrom jsonlite toJSON
#' @importFrom rhdf5 H5Fcreate H5Gcreate H5Fclose H5Gclose
setMethod("saveObject", "MultiAssayExperiment", function(x, path, ...) {
    dir.create(path, showWarnings=FALSE)
    saveObjectFile(
        path, 
        "multi_sample_dataset", 
        list(multi_sample_dataset=list(version="1.0"))
    )

    cur.exps <- experiments(x)
    if (length(cur.exps)) {
        if (anyDuplicated(names(cur.exps))) {
            stop("experiment names should be unique")
        }

        edir <- file.path(path, "experiments")
        dir.create(edir, showWarnings=FALSE)
        write(toJSON(names(cur.exps)), file=file.path(edir, "names.json"))

        for (e in seq_along(cur.exps)) {
            cur.exp <- cur.exps[[e]]
            if (anyDuplicated(colnames(cur.exp))) {
                stop("column names of 'experiments(<", class(x)[1], ">)[[", e, "]]' should be unique")
            }

            ename <- as.character(e-1L)
            tryCatch({
                saveObject(cur.exp, file.path(edir, ename), ...)
            }, error=function(e) {
                stop("failed to stage 'experiments(<", class(x)[1], ">)[[", e, "]]'\n  - ", e$message)
            })
        }
    }

    sdata <- colData(x)
    if (anyDuplicated(rownames(sdata))) {
        stop("rownames of 'colData(<", class(x)[1], ">)' should be unique")
    }
    tryCatch({
        saveObject(sdata, file.path(path, "sample_data"), ...)
    }, error=function(e) {
        stop("failed to stage 'colData(<", class(x)[1], ">)'\n  - ", e$message)
    })

    fhandle <- H5Fcreate(file.path(path, "sample_map.h5"), "H5F_ACC_TRUNC");
    on.exit(H5Fclose(fhandle))
    ghandle <- H5Gcreate(fhandle, "multi_sample_dataset")
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    map <- sampleMap(x)
    for (e in seq_along(cur.exps)) {
        keep <- map$assay == names(cur.exps)[e]
        colnames <- map$colname[keep]
        samples <- map$primary[keep]

        if (anyDuplicated(samples)) {
            stop("duplicated samples detected in 'sampleMap(<", class(x)[1], ">)'")
        }
        i <- match(samples, rownames(sdata))
        if (anyNA(i)) {
            stop("samples in 'sampleMap(<", class(x)[1], ">)' are not present in 'colData(<", class(x)[1], ">)'")
        }

        if (anyDuplicated(colnames)) {
            stop("duplicated column names detected in 'sampleMap(<", class(x)[1], ">)'")
        }
        exp.colnames <- colnames(cur.exps[[e]])
        j <- match(exp.colnames, colnames)
        if (anyNA(j)) {
            stop("column names in 'experiments(<", class(x)[1], ">)[[", e, "]]' are not present in 'sampleMap(<", class(x)[1], ">)'")
        }

        h5_write_vector(ghandle, as.character(e - 1L), (i - 1L)[j], type="H5T_NATIVE_UINT32")
    }

    saveMetadata(x, metadata.path=file.path(path, "other_data"), mcols.path=NULL, ...)
})

##################################
######### OLD STUFF HERE #########
##################################

#' @export
setMethod("stageObject", "MultiAssayExperiment", function(x, dir, path, child=FALSE, sm.name="sample_mapping", sd.name="sample_data", meta.name="dataset.json") {
    dir.create(file.path(dir, path), showWarnings=FALSE)

    # Saving the sample mappings.
    sm <- sampleMap(x)
    df <- DataFrame(
        sample = as.character(sm$primary),
        experiment = as.character(sm$assay),
        column = as.character(sm$colname)
    )

    map.info <- tryCatch({
        meta <- .stageObject(df, dir, file.path(path, sm.name), child=TRUE)
        .writeMetadata(dir=dir, meta)
    }, error=function(e) {
        stop("failed to stage 'sampleMap(<", class(x)[1], ">)'\n  - ", e$message)
    })

    # Saving the sample data.
    cd <- colData(x)
    if (is.null(rownames(cd)) || anyDuplicated(rownames(cd))) {
        stop("row names of 'colData(<", class(x)[1], ">)' should be non-NULL and unique")
    }

    sdpath <- file.path(path, sd.name)
    sd.info <- tryCatch({
        meta <- .stageObject(cd, dir, sdpath, child=TRUE)
        .writeMetadata(dir=dir, meta)
    }, error=function(e) {
        stop("failed to stage 'colData(<", class(x)[1], ">)'\n  - ", e$message)
    })

    # Saving the experiments.
    exp.info <- .stage_experiments(x, dir, path)

    # Saving other metadata.
    meta.info <- .processMetadata(x, dir, path, "metadata")

    list(
        `$schema`="dataset/v1.json",
        path=file.path(path, meta.name),
        dataset=list(
            experiments=exp.info,
            sample_mapping=list(resource=map.info),
            sample_data=list(resource=sd.info),
            other_data=meta.info
        ),
        is_child=child
    )
})

#' @import alabaster.base 
#' @importFrom MultiAssayExperiment experiments 
#' @importMethodsFrom alabaster.se stageObject
.stage_experiments <- function(x, dir, path) {
    exp.names <- names(experiments(x))
    if (anyDuplicated(exp.names)) {
        stop("detected duplicated experiment names in a ", class(x)[1], " object")
    }
    if (any(exp.names == "")) {
        stop("detected empty experiment name in a ", class(x)[1], " object")
    }

    collected <- list()
    for (i in seq_along(exp.names)) {
        newname <- file.path(path, paste0("experiment-", i))
        exp.name <- exp.names[i]
        se <- x[[exp.name]]
        if (is.null(colnames(se)) || anyDuplicated(colnames(se))) {
            stop("column names of '<", class(x)[1], ">[[", i, "]]' should be non-NULL and unique")
        }

        processed <- tryCatch({
            meta <- .stageObject(se, dir=dir, path=newname, child=TRUE)
            meta$path <- file.path(newname, "experiment.json")
            .writeMetadata(meta=meta, dir=dir)
        }, error=function(e) {
            stop("failed to stage experiment '", exp.name, "'\n  - ", e$message)
        })

        collected[[i]] <- list(name=exp.name, resource=processed)
    }

    collected
}
