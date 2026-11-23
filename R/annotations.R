## Gene-to-GO annotations: a file with one gene per line (gene <TAB> GO IDs
## separated by ';', ',', '|' or spaces), or a data frame with a gene column
## and a GO column of the same form. Annotations are propagated to all
## ancestors of the namespace.

.annotation_pairs <- function(annotations) {
    if (is.character(annotations) && length(annotations) == 1) {
        if (!file.exists(annotations))
            stop("Annotation file '", annotations, "' does not exist", call. = FALSE)
        lines <- readLines(annotations, warn = FALSE)
        lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
        has_tab <- grepl("\t", lines, fixed = TRUE)
        gene  <- ifelse(has_tab, sub("\t.*$", "", lines), sub("\\s.*$", "", lines))
        field <- ifelse(has_tab, sub("^[^\t]*\t([^\t]*).*$", "\\1", lines),
                        sub("^\\S+\\s*", "", lines))
    } else if (is.data.frame(annotations)) {
        if (ncol(annotations) < 2)
            stop("Annotation data frame needs a gene column and a GO column", call. = FALSE)
        gene  <- as.character(annotations[[1]])
        field <- as.character(annotations[[2]])
    } else {
        stop("annotations must be a file path or a data frame", call. = FALSE)
    }
    gene <- trimws(gene)
    tok  <- strsplit(trimws(field), "[;,|[:space:]]+")
    data.frame(gene = rep(gene, lengths(tok)), id = unlist(tok, use.names = FALSE))
}

.load_annotations <- function(annotations, obo, hier) {
    pairs <- .annotation_pairs(annotations)
    pairs <- pairs[startsWith(pairs$id, "GO:"), ]
    uniq  <- unique(pairs$id)
    resolved <- stats::setNames(.resolve_terms(uniq, obo), uniq)

    in_ns_or_other <- uniq %in% obo$all_ids | uniq %in% names(obo$alt_map)
    stats <- list(
        distinct_ids  = length(uniq),
        unknown_ids   = sum(!in_ns_or_other),
        replaced_ids  = sum(uniq %in% obo$obsolete & !is.na(resolved)),
        dropped_obsolete_ids = sum(uniq %in% obo$obsolete & is.na(resolved)),
        ids_in_namespace = sum(!is.na(resolved)))

    term <- unname(resolved[pairs$id])
    ok   <- !is.na(term)
    pairs <- unique(data.frame(gene = pairs$gene[ok], term = term[ok]))
    genes <- unique(pairs$gene)
    D <- Matrix::sparseMatrix(i = match(pairs$gene, genes),
                              j = match(pairs$term, obo$current),
                              x = 1, dims = c(length(genes), length(obo$current)))
    G <- D %*% hier$ancestors
    G@x[] <- 1
    dimnames(G) <- list(genes, obo$current)
    .log("Annotations: ", length(genes), " genes with terms in ", obo$namespace,
         "; GO IDs: ", stats$distinct_ids, " distinct, ", stats$replaced_ids,
         " obsolete replaced, ", stats$dropped_obsolete_ids, " obsolete dropped, ",
         stats$unknown_ids, " not in the ontology")
    list(genes = genes, G = G, stats = stats)
}
