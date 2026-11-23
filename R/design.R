## Design of one test: categories x universe genes.
## Order: match scored genes to annotations (every score row is one gene),
## level filter, size filter, universe, optional merging of redundant
## categories.

.build_design <- function(gene_ids, values, ann, hier, s, label) {
    n_scored  <- length(gene_ids)
    gi        <- match(gene_ids, ann$genes)
    annotated <- which(!is.na(gi))
    if (length(annotated) == 0)
        stop("No gene IDs of '", label, "' are in the annotations (example: '",
             gene_ids[1], "'; annotation example: '", ann$genes[1], "')", call. = FALSE)
    Gs    <- ann$G[gi[annotated], , drop = FALSE]
    sizes <- Matrix::colSums(Gs)
    n_annot <- length(annotated)

    limits <- c(if (!is.null(s$max_genes)) s$max_genes,
                if (!is.null(s$max_genes_fraction)) floor(n_annot * s$max_genes_fraction))
    max_size <- if (length(limits)) min(limits) else n_annot

    level    <- hier$level[colnames(Gs)]
    level_ok <- rep(TRUE, length(level))
    if (!is.null(s$min_level)) level_ok <- level_ok & level >= s$min_level
    if (!is.null(s$max_level)) level_ok <- level_ok & level <= s$max_level
    keep <- sizes > 0 & level_ok & sizes >= s$min_genes & sizes <= max_size

    memb <- Matrix::t(Gs[, keep, drop = FALSE])
    memb <- methods::as(memb, "CsparseMatrix")
    universe_rows <- annotated
    if (identical(s$universe, "retained")) {
        in_cat <- Matrix::colSums(memb) > 0
        memb <- memb[, in_cat, drop = FALSE]
        universe_rows <- annotated[in_cat]
    }
    term_ids <- colnames(Gs)[keep]
    design <- list(membership = memb,
                   term_ids   = term_ids,
                   merged     = rep("", length(term_ids)),
                   level      = unname(level[keep]),
                   values     = values[universe_rows],
                   genes      = gene_ids[universe_rows],
                   n_scored   = n_scored,
                   n_annotated = n_annot,
                   n_universe = length(universe_rows))
    if (isTRUE(s$reduce_redundancy) && length(term_ids) > 1) {
        design <- .merge_redundant(design, s$cluster_cut_height)
    }
    design$sizes <- as.integer(Matrix::rowSums(design$membership))
    if (length(design$term_ids) == 0)
        .log(label, ": no GO category passes the filters (category sizes ", s$min_genes, " to ",
             max_size, " genes in a universe of ", n_annot, "); no tests")
    .log(label, ": universe ", design$n_universe, " genes (", n_scored, " scored, ",
         n_annot, " annotated); ", length(design$term_ids), " categories")
    design
}

## Categories that share most of their genes are merged by complete linkage
## on 1 - shared / larger, cut at cut_height. Clusters can only form inside
## groups of categories connected by pairs within the cut height, so overlaps
## are computed in chunks, only close pairs are kept, and each connected
## group is clustered on its own. A merged category holds the union of its
## members' genes and takes the ID of its largest member (ties: higher
## level, then ID).
.merge_redundant <- function(design, cut_height, chunk_size = 500) {
    memb   <- design$membership
    n      <- nrow(memb)
    sizes  <- Matrix::rowSums(memb)
    memb_t <- Matrix::t(memb)
    i <- integer(0); j <- integer(0); d <- numeric(0)
    for (a in seq(1, n, by = chunk_size)) {
        b  <- min(n, a + chunk_size - 1)
        ov <- methods::as(memb[a:b, , drop = FALSE] %*% memb_t, "TsparseMatrix")
        ii <- ov@i + a
        jj <- ov@j + 1L
        up <- ii < jj
        ii <- ii[up]; jj <- jj[up]
        dist  <- 1 - ov@x[up] / pmax(sizes[ii], sizes[jj])
        close <- dist <= cut_height
        i <- c(i, ii[close]); j <- c(j, jj[close]); d <- c(d, dist[close])
    }
    g <- igraph::make_empty_graph(n = n, directed = FALSE)
    if (length(i)) g <- igraph::add_edges(g, as.vector(rbind(i, j)))
    comp <- igraph::components(g)$membership

    cluster <- integer(n)
    next_id <- 0L
    for (grp in unique(comp)) {
        members <- which(comp == grp)
        if (length(members) == 1) {
            next_id <- next_id + 1L
            cluster[members] <- next_id
            next
        }
        m <- length(members)
        dense <- matrix(1, m, m)
        in_grp <- i %in% members
        li <- match(i[in_grp], members)
        lj <- match(j[in_grp], members)
        dense[cbind(li, lj)] <- d[in_grp]
        dense[cbind(lj, li)] <- d[in_grp]
        diag(dense) <- 0
        labels <- stats::cutree(stats::hclust(stats::as.dist(dense), method = "complete"), h = cut_height)
        cluster[members] <- next_id + labels
        next_id <- next_id + max(labels)
    }

    ord <- function(x) x[order(-sizes[x], -design$level[x], design$term_ids[x])]
    groups  <- lapply(split(seq_len(n), cluster), ord)
    rep_idx <- vapply(groups, `[`, integer(1), 1)
    o <- order(design$term_ids[rep_idx])
    groups  <- groups[o]
    rep_idx <- rep_idx[o]

    assign <- Matrix::sparseMatrix(i = rep(seq_along(groups), lengths(groups)),
                                   j = unlist(groups), x = 1, dims = c(length(groups), n))
    new_memb <- assign %*% memb
    new_memb@x[] <- 1
    .log("Redundancy merging at cut height ", cut_height, ": ", n, " -> ",
         length(groups), " categories")
    design$membership <- methods::as(new_memb, "CsparseMatrix")
    design$merged     <- vapply(groups, function(x) paste(design$term_ids[x[-1]], collapse = ";"), character(1))
    design$term_ids   <- design$term_ids[rep_idx]
    design$level      <- design$level[rep_idx]
    design
}
