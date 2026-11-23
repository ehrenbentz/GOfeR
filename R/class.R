#' Results of a GOFER run
#'
#' Holds every result table of a run, a summary, the genes of significant
#' categories, the settings and run information. Accessors return plain data
#' frames, which can be viewed with \code{View()} in RStudio.
#'
#' @slot results Named list of result tables. Score runs: one table per
#'   score table and method (\code{"<table>|<method>"}). Module runs: one
#'   table per module (\code{"<network>|<module>"}).
#' @slot summary Number of tests and significant terms per table, method and
#'   direction (or per module and method).
#' @slot category_genes For every term significant at \code{report_fdr} in
#'   a score table or network: all universe genes annotated to the term
#'   (semicolon-separated).
#' @slot settings All option values of the run.
#' @slot info GO release, annotation statistics, output directory and, for
#'   module runs, the recurrence tables.
#' @exportClass GOFERResult
methods::setClass("GOFERResult",
    slots = c(results = "list", summary = "data.frame", category_genes = "data.frame",
              settings = "list", info = "list"))

methods::setMethod("show", "GOFERResult", function(object) {
    cat("GOFERResult (", object@info$mode, ")\n", sep = "")
    cat("  GO release:", object@info$go_version, " namespace:", object@settings$namespace, "\n")
    cat("  Result tables:", length(object@results), "\n")
    cat("  Output directory:", object@info$output_dir, "\n")
    cat("  Significant at FDR <=", object@settings$report_fdr, "(first rows of summary(x)):\n")
    print(utils::head(object@summary, 10), row.names = FALSE)
    invisible(object)
})

#' Accessors for GOFERResult objects
#'
#' Every result table is one multiple-testing family. Score runs: one table
#' per score table, method and direction. Module runs: one table per
#' network, module and method.
#'
#' @param x,object A \code{GOFERResult}.
#' @param ... For \code{results()}: the parts of a table name in order
#'   (score runs: table, method, direction; module runs: network, module,
#'   method). Fewer parts return all matching tables combined, with a
#'   \code{family} column. No parts: the names of all tables.
#' @param table For \code{recurrence()}: the network; NULL lists networks.
#' @return A data frame, or the names of the available tables.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      cpus = 1, use_cache = FALSE, output_dir = file.path(tempdir(), "example_accessors"))
#' results(res)
#' head(results(res, "example_deseq2", "mwu", "down"))
#' both <- results(res, "example_deseq2", "mwu")
#' table(both$family)
#' @name GOFERResult-accessors
NULL

#' @rdname GOFERResult-accessors
#' @export
results <- function(x, ...) {
    parts <- c(...)
    if (length(parts) == 0) return(names(x@results))
    key <- paste(parts, collapse = "|")
    if (key %in% names(x@results)) return(x@results[[key]])
    keys <- names(x@results)[startsWith(names(x@results), paste0(key, "|"))]
    if (!length(keys))
        stop("No result table '", key, "'. Available: ", paste(names(x@results), collapse = ", "), call. = FALSE)
    out <- lapply(keys, function(k) cbind(family = k, x@results[[k]]))
    cols <- unique(unlist(lapply(out, names)))
    out <- lapply(out, function(df) { df[setdiff(cols, names(df))] <- NA; df[cols] })
    do.call(rbind, out)
}

#' @rdname GOFERResult-accessors
#' @export
methods::setMethod("summary", "GOFERResult", function(object, ...) object@summary)

#' @rdname GOFERResult-accessors
#' @export
category_genes <- function(x) x@category_genes

#' @rdname GOFERResult-accessors
#' @export
settings <- function(x) x@settings

#' @rdname GOFERResult-accessors
#' @export
recurrence <- function(x, table = NULL) {
    rec <- x@info$recurrence
    if (is.null(rec)) stop("Recurrence tables exist only for module runs", call. = FALSE)
    if (is.null(table)) return(names(rec))
    rec[[table]]
}
