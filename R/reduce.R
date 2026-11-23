## Term reduction: the significant terms of one result table (one
## multiple-testing family) are grouped into themes. No p-value is changed
## and no theme is named after its most significant term.

## Genes annotated to each term in the whole annotation table (after
## propagation), for the information content of terms.
.term_gene_counts <- function(ann) {
    n <- Matrix::colSums(ann$G)
    n[n > 0]
}

## Wang et al. (2007) similarity. The S-value of an ancestor t of term A is
## the largest product of edge weights along a path from A to t
## (S_A(A) = 1); sim(A, B) = sum over common ancestors of S_A + S_B, divided
## by the sum of all S-values of A and of B. Ancestors are visited from the
## deepest level upward, so every term is final before its parents are
## updated.
.wang_similarity <- function(idx, hier) {
    n_all <- length(hier$level)
    par_rows <- vector("list", n_all)
    sp <- split(seq_len(nrow(hier$edges)), hier$edges$child)
    par_rows[as.integer(names(sp))] <- sp
    anc_t <- methods::as(Matrix::t(hier$ancestors[idx, , drop = FALSE]), "CsparseMatrix")
    sv <- numeric(n_all)
    ii <- vector("list", length(idx)); jj <- ii; xx <- ii
    for (k in seq_along(idx)) {
        anc <- anc_t@i[seq.int(anc_t@p[k] + 1, anc_t@p[k + 1])] + 1L
        anc <- anc[order(-hier$level[anc])]
        sv[idx[k]] <- 1
        for (t in anc) {
            rows <- par_rows[[t]]
            if (is.null(rows)) next
            p <- hier$edges$parent[rows]
            sv[p] <- pmax(sv[p], hier$edges$weight[rows] * sv[t])
        }
        ii[[k]] <- rep(k, length(anc)); jj[[k]] <- anc; xx[[k]] <- sv[anc]
        sv[anc] <- 0
    }
    S <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx),
                              dims = c(length(idx), n_all))
    ind <- S
    ind@x[] <- 1
    num <- as.matrix(S %*% Matrix::t(ind))
    num <- num + t(num)
    total <- Matrix::rowSums(S)
    num / outer(total, total, "+")
}

## Similarity from information content, IC(t) = -log(genes annotated to t /
## annotated genes): Lin = 2 IC(MICA) / (IC(a) + IC(b)), where MICA is the
## most informative common ancestor; Schlicker relevance = Lin x
## (1 - p(MICA)).
.ic_similarity <- function(idx, hier, ic, measure) {
    anc_t <- methods::as(Matrix::t(hier$ancestors[idx, , drop = FALSE]), "CsparseMatrix")
    mica <- .mica_ic(anc_t@p, anc_t@i, ic)
    denom <- outer(ic[idx], ic[idx], "+")
    sim <- ifelse(denom > 0, 2 * mica / denom, 0)
    if (measure == "rel") sim <- sim * (1 - exp(-mica))
    diag(sim) <- 1
    sim
}

.union_genes <- function(term_ids, gene_sets) {
    length(unique(unlist(gene_sets[term_ids], use.names = FALSE)))
}

## Themes of one family by semantic clustering (wang, ic): complete linkage
## on 1 - similarity, cut so that every pair in a theme reaches the
## similarity threshold. A theme is described by the most specific common
## ancestor of its members (fewest descendants; ties by GO ID) and by its
## medoid (member with the highest mean similarity to the others).
.cluster_themes <- function(sig, method, ont, n_desc, gene_sets, ic, s) {
    obo <- ont$obo; hier <- ont$hier
    idx <- match(sig$term_id, obo$current)
    sim <- if (method == "wang") .wang_similarity(idx, hier) else .ic_similarity(idx, hier, ic, s$ic_measure)
    n <- length(idx)
    cl <- if (n == 1) 1L else stats::cutree(stats::hclust(stats::as.dist(1 - sim), method = "complete"),
                                             h = 1 - s$similarity_threshold)
    anc <- hier$ancestors[idx, , drop = FALSE]
    rows <- list(); mem <- list()
    for (g in sort(unique(cl))) {
        m <- which(cl == g)
        common <- which(Matrix::colSums(anc[m, , drop = FALSE]) == length(m))
        common <- common[order(n_desc[common], obo$current[common])][1]
        mean_sim <- if (length(m) > 1) (rowSums(sim[m, m, drop = FALSE]) - 1) / (length(m) - 1) else 1
        med <- m[order(-mean_sim, sig$term_id[m])][1]
        rows[[g]] <- data.frame(
            theme_name = unname(obo$term_name[obo$current[common]]),
            n_terms = length(m),
            n_genes = .union_genes(sig$term_id[m], gene_sets),
            min_fdr = min(sig$fdr[m]),
            common_ancestor_id = obo$current[common],
            common_ancestor_level = unname(hier$level[common]),
            broad = unname(hier$level[common]) <= 1,
            medoid_id = sig$term_id[med],
            medoid_name = sig$term_name[med])
        mem[[g]] <- data.frame(cluster = g, term_id = sig$term_id[m], term_name = sig$term_name[m],
                               n_genes = sig$n_genes[m], fdr = sig$fdr[m],
                               similarity_to_medoid = sim[m, med])
    }
    themes <- do.call(rbind, rows)
    themes$cluster <- sort(unique(cl))
    themes <- themes[order(-themes$n_terms, -themes$n_genes, themes$common_ancestor_id), , drop = FALSE]
    themes$theme_id <- sprintf("T%03d", seq_len(nrow(themes)))
    members <- do.call(rbind, mem)
    members$theme_id <- themes$theme_id[match(members$cluster, themes$cluster)]
    members$theme_name <- themes$theme_name[match(members$cluster, themes$cluster)]
    members <- members[order(members$theme_id, members$fdr, members$term_id),
                       c("theme_id", "theme_name", "term_id", "term_name", "n_genes", "fdr", "similarity_to_medoid")]
    themes <- themes[, c("theme_id", "theme_name", "n_terms", "n_genes", "min_fdr", "common_ancestor_id",
                         "common_ancestor_level", "broad", "medoid_id", "medoid_name")]
    rownames(themes) <- NULL; rownames(members) <- NULL
    list(themes = themes, members = members)
}

## Themes of one family by GO slim mapping: every significant term is
## assigned to each term of the subset among its ancestors (or itself); the
## root is not used; terms without such an ancestor form the theme "other".
.slim_themes <- function(sig, ont, gene_sets, s) {
    obo <- ont$obo; hier <- ont$hier
    available <- sort(unique(obo$subsets$subset))
    slim <- intersect(obo$subsets$term[obo$subsets$subset == s$slim_subset], obo$current)
    slim <- slim[hier$level[slim] > 0]
    if (!length(slim))
        stop("GO subset '", s$slim_subset, "' has no terms in namespace ", obo$namespace,
             " of the GO file. Subsets in the file: ", paste(available, collapse = ", "), call. = FALSE)
    idx <- match(sig$term_id, obo$current)
    hit <- methods::as(hier$ancestors[idx, match(slim, obo$current), drop = FALSE], "TsparseMatrix")
    members <- data.frame(theme_id = slim[hit@j + 1], term = hit@i + 1)
    none <- setdiff(seq_along(idx), members$term)
    if (length(none)) members <- rbind(members, data.frame(theme_id = "other", term = none))
    members$theme_name <- ifelse(members$theme_id == "other", "other", unname(obo$term_name[members$theme_id]))
    members$term_id <- sig$term_id[members$term]
    members$term_name <- sig$term_name[members$term]
    members$n_genes <- sig$n_genes[members$term]
    members$fdr <- sig$fdr[members$term]
    by_theme <- split(members, members$theme_id)
    themes <- do.call(rbind, lapply(by_theme, function(m) data.frame(
        theme_id = m$theme_id[1], theme_name = m$theme_name[1], n_terms = nrow(m),
        n_genes = .union_genes(m$term_id, gene_sets), min_fdr = min(m$fdr))))
    themes <- themes[order(themes$theme_id == "other", -themes$n_terms, -themes$n_genes, themes$theme_id), , drop = FALSE]
    members <- members[order(match(members$theme_id, themes$theme_id), members$fdr, members$term_id),
                       c("theme_id", "theme_name", "term_id", "term_name", "n_genes", "fdr")]
    rownames(themes) <- NULL; rownames(members) <- NULL
    list(themes = themes, members = members)
}

.result_ontology <- function(x) {
    s <- x@settings
    if (!file.exists(x@info$go_file))
        stop("GO file of this result not found: '", x@info$go_file, "'", call. = FALSE)
    ont <- .cached_ontology(x@info$go_file, s$namespace, s$propagate_part_of, s$use_cache)
    if (!identical(ont$obo$data_version, x@info$go_version))
        warning("The GO file has changed since this result was computed (", x@info$go_version,
                " then, ", ont$obo$data_version, " now)", call. = FALSE)
    ont
}

.gene_sets_by_table <- function(x) {
    cg <- x@category_genes
    if (!nrow(cg)) return(list())
    lapply(split(cg, cg$table), function(d) stats::setNames(strsplit(d$genes, ";", fixed = TRUE), d$term_id))
}

.family_stem <- function(key, s) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    if (!is.null(s$output_prefix)) parts <- c(s$output_prefix, parts)
    paste(c(parts, s$namespace), collapse = "_")
}

.reduce_object <- function(x) {
    s <- x@settings
    ont <- .result_ontology(x)
    n_desc <- Matrix::colSums(ont$hier$ancestors)
    ic <- rep(0, length(ont$obo$current))
    names(ic) <- ont$obo$current
    known <- intersect(names(x@info$term_gene_counts), ont$obo$current)
    ic[known] <- -log(x@info$term_gene_counts[known] / x@info$n_annotated_genes)
    gene_sets <- .gene_sets_by_table(x)
    themes <- x@info$themes
    if (is.null(themes)) themes <- list()
    n_fam <- 0
    for (key in names(x@results)) {
        df <- x@results[[key]]
        sig <- df[!is.na(df$fdr) & df$fdr <= s$report_fdr, , drop = FALSE]
        if (!nrow(sig)) next
        if (!all(sig$term_id %in% ont$obo$current))
            stop("Result terms are missing from the GO file '", x@info$go_file, "'", call. = FALSE)
        n_fam <- n_fam + 1
        gs <- gene_sets[[strsplit(key, "|", fixed = TRUE)[[1]][1]]]
        for (m in s$reduce_methods) {
            out <- if (m == "slim") .slim_themes(sig, ont, gs, s)
                   else .cluster_themes(sig, m, ont, n_desc, gs, ic, s)
            themes[[paste(key, m, sep = "|")]] <- out
            stem <- file.path(s$output_dir, "THEMES", toupper(m), .family_stem(key, s))
            .write_tsv(out$themes, paste0(stem, "_", m, "_themes.tsv"))
            .write_tsv(out$members, paste0(stem, "_", m, "_theme_members.tsv"))
        }
    }
    .log("Term reduction (", paste(s$reduce_methods, collapse = ", "), "): ", n_fam,
         " result tables with significant terms")
    x@info$themes <- themes
    x
}

## Options given to a function that works on an existing result: checked
## like run options and merged into the settings of the result.
.update_settings <- function(x, args) {
    s <- x@settings
    unknown <- setdiff(names(args), .option_keys())
    if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    for (key in names(args)) s[key] <- list(.coerce_option(key, args[[key]], "function argument"))
    x@settings <- s
    x
}

.post_process <- function(x) {
    s <- x@settings
    if (length(s$reduce_methods)) x <- .reduce_object(x)
    if (length(s$plots)) .plot_object(x)
    x
}

#' Group significant terms into themes
#'
#' Reduces the significant terms (FDR <= \code{report_fdr}) of every result
#' table to themes. No p-value is changed, and themes are not named after
#' their most significant term. The same step runs inside
#' \code{\link{enrich_scores}} and \code{\link{enrich_modules}} when
#' \code{reduce_methods} is set.
#'
#' \describe{
#'   \item{wang}{Semantic similarity from the GO graph (Wang et al. 2007;
#'     is_a 0.8, part_of 0.6). Independent of the annotations.}
#'   \item{ic}{Semantic similarity from the information content of terms in
#'     the annotation table of the run (\code{ic_measure}: Schlicker
#'     relevance or Lin).}
#'   \item{slim}{Mapping to the terms of a GO subset (\code{slim_subset});
#'     a term can belong to several themes.}
#' }
#' For wang and ic, terms are clustered by complete linkage, so every pair
#' of terms in a theme has a similarity of at least
#' \code{similarity_threshold}. A theme is named by the most specific
#' common ancestor of its members; \code{broad} marks themes whose common
#' ancestor is the root or one of its children. The medoid is the member
#' most similar to the other members.
#'
#' @param x A \code{GOFERResult}.
#' @param reduce_methods One or more of \code{"wang"}, \code{"ic"},
#'   \code{"slim"}.
#' @param ... Further options (\code{similarity_threshold},
#'   \code{ic_measure}, \code{slim_subset}, \code{report_fdr}).
#' @return The \code{GOFERResult} with themes added; see
#'   \code{\link{themes}}. Theme and member tables are written to the
#'   \code{THEMES} folder of the output directory.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      cpus = 1, use_cache = FALSE, report_fdr = 0.1,
#'                      output_dir = file.path(tempdir(), "example_themes"))
#' res <- reduce_terms(res, "wang")
#' themes(res)
#' head(themes(res, "example_deseq2", "mwu", "down", "wang"))
#' @export
reduce_terms <- function(x, reduce_methods = c("wang", "ic", "slim"), ...) {
    x <- .update_settings(x, c(list(reduce_methods = reduce_methods), list(...)))
    x@settings$output_dir <- x@info$output_dir
    .reduce_object(x)
}

.theme_table <- function(x, parts, what) {
    th <- x@info$themes
    if (is.null(th)) stop("No themes in this result; run reduce_terms() first", call. = FALSE)
    if (!length(parts)) return(names(th))
    key <- paste(parts, collapse = "|")
    if (!key %in% names(th))
        stop("No themes for '", key, "'. Available: ", paste(names(th), collapse = ", "), call. = FALSE)
    th[[key]][[what]]
}

#' Themes of a GOFERResult
#'
#' @param x A \code{GOFERResult} with themes (see \code{\link{reduce_terms}}).
#' @param ... The name of a result table (as in \code{\link{results}})
#'   followed by the reduction method, for example
#'   \code{"F_VNO", "mwu", "up", "wang"}. No arguments: the names of all
#'   theme tables.
#' @return \code{themes()}: one row per theme. \code{theme_members()}: one
#'   row per significant term with its theme.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      cpus = 1, use_cache = FALSE, report_fdr = 0.1, reduce_methods = "slim",
#'                      output_dir = file.path(tempdir(), "example_themes2"))
#' head(theme_members(res, "example_deseq2", "mwu", "down", "slim"))
#' @export
themes <- function(x, ...) .theme_table(x, c(...), "themes")

#' @rdname themes
#' @export
theme_members <- function(x, ...) .theme_table(x, c(...), "members")
