## Parsed ontologies and propagated annotations are kept on disk in the
## user cache directory, keyed by the MD5 sum of the input file and the
## settings that change the result, so later runs skip parsing.

.cache_dir <- function() {
    d <- file.path(tools::R_user_dir("GOFER", which = "cache"), "parsed")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
}

.cache_read <- function(key) {
    f <- file.path(.cache_dir(), paste0(key, ".rds"))
    if (!file.exists(f)) return(NULL)
    tryCatch(readRDS(f), error = function(e) NULL)
}

.cache_write <- function(key, value) {
    f <- file.path(.cache_dir(), paste0(key, ".rds"))
    tmp <- paste0(f, ".tmp", Sys.getpid())
    ok <- tryCatch({ saveRDS(value, tmp, compress = FALSE); TRUE }, error = function(e) FALSE)
    if (ok) file.rename(tmp, f) else unlink(tmp)
    invisible(ok)
}

## The ontology used last is kept in memory for the rest of the R session.
.session_ontology <- new.env(parent = emptyenv())

.ontology_key <- function(go_file, namespace, propagate_part_of) {
    paste("ontology", unname(tools::md5sum(go_file)), namespace,
          if (propagate_part_of) "partof" else "isa",
          as.character(utils::packageVersion("GOFER")), "f2", sep = "_")
}

.cached_ontology <- function(go_file, namespace, propagate_part_of, use_cache) {
    key <- .ontology_key(go_file, namespace, propagate_part_of)
    if (identical(.session_ontology$key, key)) return(.session_ontology$ont)
    on.exit({
        .session_ontology$key <- key
        .session_ontology$ont <- ont
    })
    ont <- if (use_cache) .cache_read(key) else NULL
    if (is.null(ont)) {
        obo  <- .read_obo(go_file, namespace)
        hier <- .build_hierarchy(obo, propagate_part_of)
        ont  <- list(obo = obo, hier = hier)
        if (use_cache) .cache_write(key, ont)
        .log("GO release ", obo$data_version, ": ", length(obo$current), " current ", namespace,
             " terms, ", hier$n_edges, " edges (parsed and cached)")
    } else {
        .log("GO release ", ont$obo$data_version, ": ", length(ont$obo$current), " current ", namespace,
             " terms, ", ont$hier$n_edges, " edges (from cache)")
    }
    ont
}

.cached_annotations <- function(annotations, run) {
    if (!run$s$use_cache || !(is.character(annotations) && length(annotations) == 1 && file.exists(annotations))) {
        return(.load_annotations(annotations, run$obo, run$hier))
    }
    key <- paste("annotations", unname(tools::md5sum(annotations)),
                 .ontology_key(run$go_file, run$s$namespace, run$s$propagate_part_of), sep = "_")
    ann <- .cache_read(key)
    if (is.null(ann)) {
        ann <- .load_annotations(annotations, run$obo, run$hier)
        .cache_write(key, ann)
    } else {
        .log("Annotations: ", length(ann$genes), " genes with terms in ", run$s$namespace, " (from cache)")
    }
    ann
}
