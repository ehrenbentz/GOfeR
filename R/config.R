## Configuration files: one option per line, "key = value", '#' starts a
## comment line. Values: TRUE / FALSE, numbers, text, comma-separated lists,
## NULL for not set.

.read_config <- function(path) {
    if (!file.exists(path)) stop("Config file '", path, "' does not exist", call. = FALSE)
    lines <- readLines(path, warn = FALSE)
    keys <- .option_keys()
    values <- list()
    for (n in seq_along(lines)) {
        line <- trimws(lines[n])
        if (!nzchar(line) || startsWith(line, "#")) next
        if (!grepl("=", line, fixed = TRUE))
            stop("Config file '", path, "', line ", n, ": expected 'key = value'", call. = FALSE)
        key  <- trimws(sub("=.*$", "", line))
        text <- trimws(sub("^[^=]*=", "", line))
        key  <- gsub("-", "_", key, fixed = TRUE)
        if (!key %in% keys)
            stop("Config file '", path, "', line ", n, ": unknown option '", key, "'", call. = FALSE)
        if (key %in% names(values))
            stop("Config file '", path, "', line ", n, ": option '", key, "' is given twice", call. = FALSE)
        values[key] <- list(.coerce_option(key, text, paste0("config file line ", n)))
    }
    values
}

## Defaults < config file < arguments of the call. Returns the settings and
## where each value came from.
.resolve_settings <- function(call_args, config) {
    settings <- .default_settings()
    origin <- stats::setNames(rep("default", length(settings)), names(settings))
    if (!is.null(config)) {
        cfg <- .read_config(config)
        for (key in names(cfg)) {
            settings[key] <- list(cfg[[key]])
            origin[key] <- "config file"
        }
    }
    unknown <- setdiff(names(call_args), names(settings))
    if (length(unknown) > 0)
        stop("Unknown option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    for (key in names(call_args)) {
        settings[key] <- list(.coerce_option(key, call_args[[key]], "function argument"))
        origin[key] <- "function argument"
    }
    list(settings = settings, origin = origin)
}

.wrap_comment <- function(text) {
    paste0("# ", strwrap(text, width = 76))
}

.render_settings <- function(settings, origin = NULL, applies = c("scores", "modules"), header = character(0)) {
    opts <- .options()
    width <- max(nchar(.option_keys())) + 1
    out <- if (length(header)) paste0("# ", header) else character(0)
    for (section in c("basic", "advanced")) {
        out <- c(out, "", strrep("#", 70),
                 paste0("# ", if (section == "basic") "Basic options" else "Advanced options"),
                 strrep("#", 70))
        for (opt in opts) {
            if (opt$section != section) next
            if (!(opt$applies == "both" || opt$applies %in% applies)) next
            out <- c(out, "", .wrap_comment(opt$help))
            if (!is.null(origin))
                out <- c(out, paste0("# [set by: ", origin[[opt$key]], "; default: ",
                                     .format_value(opt$default), "]"))
            out <- c(out, paste0(formatC(opt$key, width = -width), "= ",
                                 .format_value(settings[[opt$key]])))
        }
    }
    c(out, "")
}

#' Write the default configuration file
#'
#' Writes a configuration file with every option at its default value and a
#' description of each option. Edit the values and pass the file to
#' \code{enrich_scores()} or \code{enrich_modules()} with \code{config =}.
#' Precedence: defaults < config file < function arguments.
#'
#' @param path File to write.
#' @return \code{path}, invisibly.
#' @examples
#' f <- file.path(tempdir(), "gofer_config.txt")
#' write_config(f)
#' head(readLines(f), 20)
#' @export
write_config <- function(path = "gofer_config.txt") {
    header <- c("GOFER configuration file", "",
                "One option per line: key = value. Lines starting with '#' are comments.",
                "TRUE / FALSE for switches; comma-separated values for lists; NULL for not set.",
                "Precedence: built-in defaults < this file < arguments of the function call.")
    writeLines(.render_settings(.default_settings(), header = header), path)
    invisible(path)
}
