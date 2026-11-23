## Logging, parallel backend, multiple testing, table writing.

.log_env <- new.env(parent = emptyenv())

.log <- function(...) {
    line <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
    message(line)
    if (!is.null(.log_env$file)) write(line, file = .log_env$file, append = TRUE)
    invisible(line)
}

## cpus NULL: all physical cores. backend auto: forked workers on Linux and
## macOS, separate R sessions on Windows. Separate-session workers are kept
## for the rest of the R session once started, so later runs with the same
## number of workers do not start them again (see stop_workers()).
.worker_env <- new.env(parent = emptyenv())

.make_bpparam <- function(cpus, backend) {
    n <- cpus
    if (is.null(n)) {
        n <- parallel::detectCores(logical = FALSE)
        if (is.na(n)) n <- parallel::detectCores()
        if (is.na(n)) n <- 1L
    }
    if (n <= 1 || backend == "serial") return(BiocParallel::SerialParam())
    if (backend == "snow" || (backend == "auto" && .Platform$OS.type == "windows")) {
        bp <- .worker_env$snow
        if (!is.null(bp) && BiocParallel::bpisup(bp) && BiocParallel::bpnworkers(bp) == n) return(bp)
        stop_workers()
        bp <- BiocParallel::SnowParam(workers = n, type = "SOCK")
        .worker_env$snow <- bp
        return(bp)
    }
    BiocParallel::MulticoreParam(workers = n)
}

#' Stop the worker sessions kept by GOFER
#'
#' On Windows (or with \code{backend = "snow"}) GOFER starts separate R
#' sessions as workers the first time a run needs them and keeps them for
#' later runs in the same R session. They end with the R session; this
#' function ends them earlier.
#'
#' @return \code{NULL}, invisibly.
#' @examples
#' stop_workers()
#' @export
stop_workers <- function() {
    bp <- .worker_env$snow
    if (!is.null(bp) && BiocParallel::bpisup(bp)) BiocParallel::bpstop(bp)
    .worker_env$snow <- NULL
    invisible(NULL)
}

.bh <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- !is.na(p)
    out[ok] <- stats::p.adjust(p[ok], method = "BH")
    out
}

.first <- function(x) if (length(x)) x[1] else NA

## Tab-separated table; numbers that are not whole numbers are written with
## four significant digits.
.write_tsv <- function(df, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    cols <- lapply(df, function(x) {
        if (is.double(x) && any(is.finite(x) & x != round(x))) x <- signif(x, 4)
        y <- as.character(x)
        y[is.na(x)] <- "NA"
        if (is.character(x)) y <- gsub("\t", " ", y, fixed = TRUE)
        y
    })
    lines <- if (nrow(df)) do.call(paste, c(unname(cols), sep = "\t")) else character(0)
    writeLines(c(paste(names(df), collapse = "\t"), lines), path)
}
