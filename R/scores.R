## Score tables: DESeq2 results (DESeqResults objects or tables written with
## write.table), or any table with gene IDs in the first column or in the row
## names.

.signed_score_types <- c("wald", "log2foldchange", "signed_pvalue", "pval_signedlog2fc", "raw")

## Separator: the first of tab, comma, space found in the header line.
## Tables written by write.table() with row names have one header field
## fewer than data fields.
.read_score_file <- function(path) {
    first2 <- readLines(path, n = 2, warn = FALSE)
    if (length(first2) == 0 || !nzchar(trimws(first2[1])))
        stop("Score file '", path, "' is empty", call. = FALSE)
    sep <- "\t"
    for (cand in c("\t", ",", " ")) if (grepl(cand, first2[1], fixed = TRUE)) { sep <- cand; break }
    n_head <- length(strsplit(first2[1], sep, fixed = TRUE)[[1]])
    n_data <- if (length(first2) > 1) length(strsplit(first2[2], sep, fixed = TRUE)[[1]]) else 0
    tab <- utils::read.table(path, sep = sep, header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE, quote = "\"", comment.char = "",
                             row.names = if (n_data == n_head + 1) 1 else NULL)
    if (n_data == n_head + 1) {
        tab <- data.frame(gene = rownames(tab), tab, check.names = FALSE, row.names = NULL)
    } else {
        names(tab)[1] <- "gene"
        tab$gene <- as.character(tab$gene)
    }
    tab
}

.as_score_table <- function(x, label) {
    if (inherits(x, "DataFrame") || inherits(x, "DESeqResults")) {
        x <- as.data.frame(x)
    }
    if (is.data.frame(x)) {
        if (!is.null(rownames(x)) && !identical(rownames(x), as.character(seq_len(nrow(x))))) {
            return(data.frame(gene = rownames(x), x, check.names = FALSE, row.names = NULL))
        }
        names(x)[1] <- "gene"
        x$gene <- as.character(x$gene)
        return(x)
    }
    if (is.character(x) && length(x) == 1 && file.exists(x)) return(.read_score_file(x))
    stop("Score input '", label, "' is not a file, a data frame or a DESeqResults object", call. = FALSE)
}

## Named list of score tables from a file, a directory, a data frame, a
## DESeqResults object, or a named list of these.
.collect_score_tables <- function(scores) {
    if (is.character(scores) && length(scores) == 1 && dir.exists(scores)) {
        files <- sort(list.files(scores, full.names = TRUE))
        files <- files[!startsWith(basename(files), ".") & !dir.exists(files)]
        if (length(files) == 0) stop("Directory '", scores, "' contains no files", call. = FALSE)
        scores <- stats::setNames(as.list(files), tools::file_path_sans_ext(basename(files)))
    } else if (is.character(scores)) {
        missing <- scores[!file.exists(scores)]
        if (length(missing)) stop("Score file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
        scores <- stats::setNames(as.list(scores), tools::file_path_sans_ext(basename(scores)))
    } else if (!is.list(scores) || is.data.frame(scores)) {
        scores <- list(scores = scores)
    }
    if (is.null(names(scores))) names(scores) <- rep("", length(scores))
    unnamed <- !nzchar(names(scores))
    names(scores)[unnamed] <- paste0("table_", which(unnamed))
    if (anyDuplicated(names(scores)))
        stop("Score table names must be unique: ", paste(unique(names(scores)[duplicated(names(scores))]), collapse = ", "), call. = FALSE)
    lapply(stats::setNames(names(scores), names(scores)),
           function(nm) .as_score_table(scores[[nm]], nm))
}

.numeric_column <- function(tab, column, role, label) {
    if (is.null(column)) stop(role, ": no column name given", call. = FALSE)
    if (!column %in% names(tab))
        stop(role, " '", column, "' not found in '", label, "'. Columns: ",
             paste(names(tab)[-1], collapse = ", "), call. = FALSE)
    v <- suppressWarnings(as.numeric(tab[[column]]))
    if (all(is.na(v))) stop(role, " '", column, "' in '", label, "' is not numeric", call. = FALSE)
    v
}

.pvalue_column <- function(tab, column, label) {
    p <- .numeric_column(tab, column, "P-value column", label)
    if (any(p < 0 | p > 1, na.rm = TRUE))
        stop("P-value column '", column, "' in '", label, "' has values outside [0, 1]", call. = FALSE)
    p
}

## Gene scores for mwu and pr (value), and p-values for fisher (pvalue).
.gene_scores <- function(tab, s, label) {
    if (s$respect_independent_filtering) {
        cols <- names(tab)[-1][grepl("pvalue|padj", tolower(names(tab)[-1]))]
        if (length(cols)) {
            keep <- stats::complete.cases(tab[, cols, drop = FALSE])
            .log(label, ": independent filtering removed ", sum(!keep), " of ", nrow(tab), " genes")
            tab <- tab[keep, , drop = FALSE]
        }
    }
    neglog <- function(p) -log10(pmax(p, 1e-300))
    value <- NULL
    if (any(s$methods %in% c("mwu", "pr"))) {
        value <- switch(s$score_type,
            pvalue  = neglog(.pvalue_column(tab, s$pvalue_column, label)),
            wald    = .numeric_column(tab, s$wald_column, "Wald statistic column", label),
            log2foldchange = .numeric_column(tab, s$fold_change_column, "Fold-change column", label),
            signed_pvalue = neglog(.pvalue_column(tab, s$pvalue_column, label)) *
                sign(.numeric_column(tab, s$fold_change_column, "Fold-change column", label)),
            pval_signedlog2fc = neglog(.pvalue_column(tab, s$pvalue_column, label)) *
                .numeric_column(tab, s$fold_change_column, "Fold-change column", label),
            raw = .numeric_column(tab, s$score_column, "Score column", label))
    }
    pvalue <- if ("fisher" %in% s$methods) .pvalue_column(tab, s$pvalue_column, label) else NULL
    dup <- table(tab$gene)
    dup <- dup[dup > 1]
    list(gene = tab$gene, value = value, pvalue = pvalue,
         duplicated = data.frame(gene = names(dup), n_rows = as.integer(dup)))
}
