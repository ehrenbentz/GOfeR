## GO ontology: download, OBO parsing, term resolution, ancestors and levels.

.go_url <- "https://current.geneontology.org/ontology/go.obo"
.go_max_age_days <- 30

## Path of the OBO file to use. A missing file is downloaded; a file older
## than 30 days is replaced by the current release unless retain is TRUE.
## The download goes to a temporary file first, so a failed download keeps
## an existing file.
.ensure_go_obo <- function(go_file, retain) {
    if (is.null(go_file)) {
        dir <- tools::R_user_dir("GOFER", which = "cache")
        go_file <- file.path(dir, "go.obo")
    }
    exists <- file.exists(go_file)
    if (exists) {
        age <- as.numeric(difftime(Sys.time(), file.mtime(go_file), units = "days"))
        if (retain || age < .go_max_age_days) {
            .log("Using GO file '", go_file, "' (", round(age), " days old)")
            return(go_file)
        }
        .log("GO file '", go_file, "' is ", round(age), " days old; downloading the current release")
    } else {
        .log("GO file '", go_file, "' not found; downloading ", .go_url)
    }
    dir.create(dirname(go_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = ".go_obo_", tmpdir = dirname(go_file))
    ok <- tryCatch({
        utils::download.file(.go_url, tmp, mode = "wb", quiet = TRUE)
        TRUE
    }, error = function(e) {
        .log("Download failed: ", conditionMessage(e))
        FALSE
    })
    if (ok && file.exists(tmp) && file.size(tmp) > 0) {
        file.rename(tmp, go_file)
        .log("Downloaded the GO ontology to '", go_file, "'")
        return(go_file)
    }
    if (file.exists(tmp)) unlink(tmp)
    if (exists) {
        .log("Using the existing GO file '", go_file, "'")
        return(go_file)
    }
    stop("Could not download the GO ontology. Download go.obo from ", .go_url,
         " and pass it with go_file.", call. = FALSE)
}

## All [Term] stanzas of an OBO file; namespace_code selects the current
## (non-obsolete) terms of one namespace. Other stanza types are skipped.
.read_obo <- function(path, namespace_code) {
    ns_name <- c(BP = "biological_process", MF = "molecular_function",
                 CC = "cellular_component")[[namespace_code]]

    lines     <- trimws(readLines(path, warn = FALSE))
    is_header <- startsWith(lines, "[") & endsWith(lines, "]")
    stanza_no <- cumsum(is_header)
    first     <- which(is_header)[1]
    header    <- lines[seq_len(first - 1)]
    dv        <- header[startsWith(header, "data-version:")]
    data_version <- if (length(dv)) sub("^data-version:\\s*", "", dv[1]) else NA_character_

    term_stanzas <- which(lines[is_header] == "[Term]")
    keep <- !is_header & (stanza_no %in% term_stanzas) & nzchar(lines) & !startsWith(lines, "!")
    st   <- stanza_no[keep]
    ln   <- lines[keep]
    tag  <- sub(":.*$", "", ln)
    val  <- sub("^[^:]*:\\s*", "", ln)
    val  <- sub(" ! .*$", "", val)

    id_lines <- tag == "id"
    id_of    <- stats::setNames(val[id_lines], st[id_lines])
    term_of  <- unname(id_of[as.character(st)])
    all_ids  <- unique(val[id_lines])

    namespace <- stats::setNames(val[tag == "namespace"], term_of[tag == "namespace"])
    term_name <- stats::setNames(val[tag == "name"], term_of[tag == "name"])
    obsolete  <- unique(term_of[tag == "is_obsolete" & val == "true"])

    alt_map <- stats::setNames(term_of[tag == "alt_id"], val[tag == "alt_id"])
    alt_map <- alt_map[!duplicated(names(alt_map))]

    repl_from <- term_of[tag == "replaced_by"]
    repl_to   <- val[tag == "replaced_by"]
    n_repl    <- table(repl_from)
    single    <- names(n_repl)[n_repl == 1]
    repl_map  <- stats::setNames(repl_to[repl_from %in% single], repl_from[repl_from %in% single])

    current <- sort(setdiff(names(namespace)[namespace == ns_name], obsolete))

    rel_val   <- val[tag == "relationship"]
    rel_child <- term_of[tag == "relationship"]
    is_partof <- startsWith(rel_val, "part_of ")

    list(data_version = data_version,
         namespace    = namespace_code,
         current      = current,
         all_ids      = all_ids,
         term_name    = term_name,
         obsolete     = obsolete,
         alt_map      = alt_map,
         repl_map     = repl_map,
         is_a         = data.frame(child = term_of[tag == "is_a"],
                                   parent = sub("\\s.*$", "", val[tag == "is_a"])),
         part_of      = data.frame(child = rel_child[is_partof],
                                   parent = sub("^part_of\\s+(\\S+).*$", "\\1", rel_val[is_partof])),
         subsets      = data.frame(term = term_of[tag == "subset"], subset = val[tag == "subset"]))
}

## Ancestor matrix A[term, ancestor] (1 = ancestor or the term itself),
## level of every term (longest path from a root; root = 0) and the edges as
## term indices with their semantic weight (is_a 0.8, part_of 0.6; is_a when
## both relations link the same pair).
.build_hierarchy <- function(obo, propagate_part_of) {
    edges <- cbind(obo$is_a, weight = rep(0.8, nrow(obo$is_a)))
    if (propagate_part_of) edges <- rbind(edges, cbind(obo$part_of, weight = rep(0.6, nrow(obo$part_of))))
    edges <- edges[edges$child %in% obo$current & edges$parent %in% obo$current, ]
    edges <- edges[!duplicated(edges[, c("child", "parent")]), ]

    n  <- length(obo$current)
    ci <- match(edges$child, obo$current)
    pi <- match(edges$parent, obo$current)

    P <- Matrix::sparseMatrix(i = ci, j = pi, x = 1, dims = c(n, n))
    A <- methods::as(Matrix::Diagonal(n, x = 1), "CsparseMatrix")
    repeat {
        A_new <- A + A %*% P
        A_new@x[] <- 1
        A_new <- Matrix::drop0(A_new)
        if (Matrix::nnzero(A_new) == Matrix::nnzero(A)) break
        A <- A_new
    }
    dimnames(A) <- list(obo$current, obo$current)

    level <- rep(0, n)
    repeat {
        cand <- tapply(level[pi] + 1, ci, max)
        idx  <- as.integer(names(cand))
        new_level <- level
        new_level[idx] <- pmax(level[idx], cand)
        if (all(new_level == level)) break
        level <- new_level
    }
    names(level) <- obo$current
    list(ancestors = A, level = level, n_edges = nrow(edges),
         edges = data.frame(child = ci, parent = pi, weight = edges$weight))
}

## Map GO IDs to current terms of the namespace: alt_id to the primary ID,
## obsolete terms with a single replaced_by to the replacement; anything
## else that is not a current term of the namespace becomes NA.
.resolve_terms <- function(tokens, obo) {
    primary <- tokens
    alt_hit <- !(primary %in% obo$all_ids) & primary %in% names(obo$alt_map)
    primary[alt_hit] <- obo$alt_map[primary[alt_hit]]

    obs    <- primary %in% obo$obsolete
    rep_to <- obo$repl_map[primary[obs]]
    rep_alt <- !is.na(rep_to) & !(rep_to %in% obo$all_ids) & rep_to %in% names(obo$alt_map)
    rep_to[rep_alt] <- obo$alt_map[rep_to[rep_alt]]
    rep_to[rep_to %in% obo$obsolete] <- NA
    primary[obs] <- rep_to

    primary[!(primary %in% obo$current)] <- NA
    unname(primary)
}
