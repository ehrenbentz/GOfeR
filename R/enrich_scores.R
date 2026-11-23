## Directions tested for mwu and pr, and the score type: when score_type is
## NULL, wald for directional tests and pvalue for enriched. The ranking by
## -log10(p) of a Wald test is the ranking by |Wald statistic|, and the
## ranking by signed_pvalue is the ranking by the Wald statistic.
.score_plan <- function(s) {
    d <- s$test_direction
    if (is.null(s$score_type)) s$score_type <- if (d == "enriched") "pvalue" else "wald"
    signed <- s$score_type %in% .signed_score_types
    if (signed && d == "enriched")
        stop("score_type ", s$score_type, " is signed; test_direction enriched needs an unsigned score (pvalue)", call. = FALSE)
    if (!signed && d != "enriched")
        stop("score_type pvalue is unsigned; test_direction ", d, " needs a signed score (wald or signed_pvalue)", call. = FALSE)
    s$directions <- if (d == "both") c("up", "down") else d
    s
}

.check_common <- function(s) {
    if (s$min_genes < 1) stop("min_genes must be >= 1", call. = FALSE)
    if (!is.null(s$max_genes) && s$max_genes < s$min_genes)
        stop("max_genes cannot be smaller than min_genes", call. = FALSE)
    if (!is.null(s$max_genes_fraction) && !(s$max_genes_fraction > 0 && s$max_genes_fraction <= 1))
        stop("max_genes_fraction must be in (0, 1], or NULL", call. = FALSE)
    if (!is.null(s$min_level) && !is.null(s$max_level) && s$min_level > s$max_level)
        stop("min_level cannot be greater than max_level", call. = FALSE)
    if (!(s$cluster_cut_height > 0 && s$cluster_cut_height < 1))
        stop("cluster_cut_height must be between 0 and 1", call. = FALSE)
    if (!is.null(s$cpus) && s$cpus < 1) stop("cpus must be >= 1, or NULL", call. = FALSE)
    .plot_colours(s)
    if (s$pr_resamples < 1000) stop("pr_resamples must be >= 1000", call. = FALSE)
    if (s$pr_max_resamples < s$pr_resamples)
        stop("pr_max_resamples must be >= pr_resamples", call. = FALSE)
    if (s$pr_tail_size < 10 || s$pr_tail_size >= s$pr_resamples)
        stop("pr_tail_size must be >= 10 and smaller than pr_resamples", call. = FALSE)
    for (key in c("namespace", "universe", "output_dir", "backend", "propagate_part_of",
                  "reduce_redundancy", "retain_go_obo", "write_gene_lists", "use_cache")) {
        if (is.null(s[[key]])) stop("Option '", key, "' cannot be NULL", call. = FALSE)
    }
}

## Shared start of a run: settings, output directory, log, workers,
## ontology and annotations.
.start_run <- function(call_args, config, applies) {
    resolved <- .resolve_settings(call_args, config)
    s <- resolved$settings
    if (is.null(s$pr_seed)) {
        s$pr_seed <- sample.int(1e6, 1)
        resolved$origin["pr_seed"] <- "drawn from R's random number generator"
    }
    .check_common(s)
    out_dir <- s$output_dir
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    .log_env$file <- file.path(out_dir, "run_log.txt")
    file.create(.log_env$file)
    .log("GOFER ", as.character(utils::packageVersion("GOFER")), "; output directory '",
         normalizePath(out_dir), "'")
    bp <- .make_bpparam(s$cpus, s$backend)
    go_file <- .ensure_go_obo(s$go_file, s$retain_go_obo)
    ont  <- .cached_ontology(go_file, s$namespace, s$propagate_part_of, s$use_cache)
    obo  <- ont$obo
    hier <- ont$hier
    list(s = s, origin = resolved$origin, bp = bp, obo = obo, hier = hier, go_file = go_file,
         applies = applies)
}

.finish_run <- function(run) {
    if (methods::is(run$bp, "MulticoreParam") && BiocParallel::bpisup(run$bp)) BiocParallel::bpstop(run$bp)
    header <- c("GOFER settings of this run", "",
                paste0("GOFER version: ", as.character(utils::packageVersion("GOFER"))),
                paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                paste0("GO file: ", normalizePath(run$go_file), " (data-version ", run$obo$data_version, ")"),
                "This file can be passed back as config = to repeat the run.")
    writeLines(.render_settings(run$s, run$origin, applies = run$applies, header = header),
               file.path(run$s$output_dir, "settings.txt"))
    .log("Settings written to '", file.path(run$s$output_dir, "settings.txt"), "'")
    .log_env$file <- NULL
}

## Columns shared by every result table.
.base_columns <- function(design, obo) {
    data.frame(term_id = design$term_ids,
               term_name = unname(obo$term_name[design$term_ids]),
               level = design$level,
               n_genes = design$sizes)
}

## Completes one family (one multiple-testing correction): FDR, universe
## size, merged terms; rows without a p-value are dropped. Column order:
## term columns, statistics, p-value, FDR, universe counts. Sorted by FDR,
## p-value and term.
.finish_family <- function(df, design, s) {
    df <- df[!is.na(df$p_value), , drop = FALSE]
    df$fdr <- .bh(df$p_value)
    df$n_universe <- rep(design$n_universe, nrow(df))
    if (isTRUE(s$reduce_redundancy)) df$merged_terms <- design$merged[match(df$term_id, design$term_ids)]
    last <- intersect(c("p_value", "p_value_method", "fdr", "n_universe", "n_sig_universe", "merged_terms"), names(df))
    df <- df[, c(setdiff(names(df), last), last), drop = FALSE]
    df <- df[order(df$fdr, df$p_value, df$term_id), , drop = FALSE]
    rownames(df) <- NULL
    df
}

## For each listed term, all universe genes annotated to it.
.term_genes <- function(design, terms) {
    rows <- which(design$term_ids %in% terms)
    if (length(rows) == 0) return(data.frame(term = character(0), genes = character(0)))
    sub <- methods::as(design$membership[rows, , drop = FALSE], "TsparseMatrix")
    genes <- tapply(design$genes[sub@j + 1], sub@i + 1, function(g) paste(sort(unique(g)), collapse = ";"))
    data.frame(term_id = design$term_ids[rows][as.integer(names(genes))], genes = unname(as.character(genes)))
}

#' GO enrichment of gene score tables
#'
#' Tests GO categories for enrichment in one or more tables of gene scores,
#' for example DESeq2 results. Every option has a default; a config file
#' (see \code{\link{write_config}}) overrides the defaults, and options given
#' as arguments override the config file. All result tables, the settings of
#' the run and a log are written to \code{output_dir}.
#'
#' @param scores A file, a directory of files, a data frame, a DESeqResults
#'   object, or a named list of these. Gene IDs are taken from the row names
#'   (DESeqResults, or tables written with \code{write.table}) or from the
#'   first column.
#' @param annotations Gene-to-GO annotations: a file with lines
#'   \code{gene <TAB> GO:...;GO:...}, or a data frame with a gene column and
#'   a GO column of the same form.
#' @param ... Options as \code{key = value}; see \code{\link{write_config}}
#'   for the full list.
#' @param config Optional configuration file.
#' @return A \code{GOFERResult} object (see \code{\link{GOFERResult-class}}).
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      methods = "mwu,fisher", cpus = 1, use_cache = FALSE,
#'                      output_dir = file.path(tempdir(), "example_scores"))
#' res
#' head(results(res, "example_deseq2", "mwu", "up"))
#' @export
enrich_scores <- function(scores, annotations, ..., config = NULL) {
    run <- .start_run(list(...), config, "scores")
    s <- run$s
    on.exit(.finish_run(run), add = TRUE)
    s <- .score_plan(s)
    run$s$score_type <- s$score_type
    directions <- s$directions
    if ("fisher" %in% s$methods && !(s$threshold > 0 && s$threshold <= 1))
        stop("threshold must be in (0, 1]", call. = FALSE)
    if (s$score_type == "raw" && is.null(s$score_column))
        stop("score_type raw requires score_column", call. = FALSE)

    ann <- .cached_annotations(annotations, run)
    tables <- .collect_score_tables(scores)
    results <- list(); pr_families <- list(); designs <- list()

    for (tab_name in names(tables)) {
        prefix <- if (is.null(s$output_prefix)) tab_name else paste0(s$output_prefix, "_", tab_name)
        gs <- .gene_scores(tables[[tab_name]], s, tab_name)
        if (nrow(gs$duplicated)) {
            .write_tsv(gs$duplicated, file.path(s$output_dir, "DUPLICATED_GENE_IDS",
                                                paste0(prefix, "_duplicated_gene_ids.tsv")))
            .log(tab_name, ": ", nrow(gs$duplicated), " gene IDs occur on more than one row; ",
                 "each row is analyzed as a separate gene")
        }

        if (any(s$methods %in% c("mwu", "pr"))) {
            ok <- !is.na(gs$value)
            if (s$score_type %in% c("wald", "raw") && all(gs$value[ok] >= 0))
                stop("All ", s$score_type, " scores in '", tab_name, "' are >= 0, so they carry no direction ",
                     "(DESeq2 LRT results have an unsigned stat column). Use score_type signed_pvalue, ",
                     "or test_direction enriched with score_type pvalue.", call. = FALSE)
            d <- .build_design(gs$gene[ok], gs$value[ok], ann, run$hier, s, tab_name)
            designs[[paste0(tab_name, "|score")]] <- d
            base <- .base_columns(d, run$obo)
            if ("mwu" %in% s$methods) {
                if (isTRUE(s$two_sided)) {
                    res <- .mwu_tests(d$membership, d$values, "two.sided")
                    df <- base
                    df$direction <- ifelse(res$delta_rank > 0, "up", ifelse(res$delta_rank < 0, "down", "none"))
                    df$delta_rank <- res$delta_rank
                    df$p_value <- res$two.sided
                    results[[paste0(tab_name, "|mwu|two_sided")]] <- .finish_family(df, d, s)
                } else {
                    alts <- c(up = "greater", down = "less", enriched = "greater")[directions]
                    res <- .mwu_tests(d$membership, d$values, unique(unname(alts)))
                    for (dir in directions) {
                        df <- base
                        df$delta_rank <- res$delta_rank
                        df$p_value <- res[[alts[[dir]]]]
                        results[[paste0(tab_name, "|mwu|", dir)]] <- .finish_family(df, d, s)
                    }
                }
            }
            if ("pr" %in% s$methods) {
                for (dir in directions) {
                    v <- if (dir == "down") -d$values else d$values
                    ap <- .average_precision(d$membership, v)
                    pr_families[[paste0(tab_name, "|pr|", dir)]] <- list(
                        ap = ap, sizes = d$sizes, n_universe = d$n_universe, base = base,
                        design_key = paste0(tab_name, "|score"))
                }
            }
        }
        if ("fisher" %in% s$methods) {
            ok <- !is.na(gs$pvalue)
            d <- .build_design(gs$gene[ok], gs$pvalue[ok], ann, run$hier, s, paste0(tab_name, " (fisher)"))
            designs[[paste0(tab_name, "|fisher")]] <- d
            significant <- d$values <= s$threshold
            df <- .base_columns(d, run$obo)
            if (sum(significant) == 0 || sum(significant) == length(significant)) {
                if (length(d$term_ids)) .log(tab_name, ": fisher not possible (all genes in one class)")
                df$p_value <- rep(NA_real_, nrow(df))
            } else {
                fr <- .fisher_tests(d$membership, significant)
                df$n_sig_genes <- fr$sig_genes
                df$enrichment_ratio <- fr$enrichment_ratio
                df$odds_ratio <- fr$odds_ratio
                df$p_value <- fr$pval
                df$n_sig_universe <- rep(sum(significant), nrow(df))
            }
            results[[paste0(tab_name, "|fisher|enriched")]] <- .finish_family(df, d, s)
        }
    }

    if (length(pr_families)) {
        pv <- .pr_families(pr_families, s, run$bp)
        for (k in names(pr_families)) {
            f <- pr_families[[k]]
            df <- f$base
            df$average_precision <- f$ap
            df$null_exceedances <- pv[[k]]$exceedances
            df$null_resamples <- pv[[k]]$resamples
            df$p_value <- pv[[k]]$p
            df$p_value_method <- pv[[k]]$method
            results[[k]] <- .finish_family(df, designs[[f$design_key]], s)
        }
    }

    summary_rows <- list(); sig_terms <- list()
    for (key in names(results)) {
        parts <- strsplit(key, "|", fixed = TRUE)[[1]]
        tab_name <- parts[1]; method <- parts[2]; dir <- parts[3]
        prefix <- if (is.null(s$output_prefix)) tab_name else paste0(s$output_prefix, "_", tab_name)
        df <- results[[key]]
        .write_tsv(df, file.path(s$output_dir, toupper(method), toupper(dir),
                                 paste0(prefix, "_", method, "_", dir, "_", s$namespace, ".tsv")))
        summary_rows[[key]] <- data.frame(table = tab_name, method = method, direction = dir,
                                          n_tests = nrow(df), n_universe = .first(df$n_universe),
                                          n_significant = sum(df$fdr <= s$report_fdr, na.rm = TRUE),
                                          fdr_cutoff = s$report_fdr)
        sig_terms[[tab_name]] <- union(sig_terms[[tab_name]], df$term_id[!is.na(df$fdr) & df$fdr <= s$report_fdr])
    }
    term_genes <- list()
    for (tab_name in names(sig_terms)) {
        d <- designs[[paste0(tab_name, "|score")]]
        if (is.null(d)) d <- designs[[paste0(tab_name, "|fisher")]]
        tg <- .term_genes(d, sig_terms[[tab_name]])
        prefix <- if (is.null(s$output_prefix)) tab_name else paste0(s$output_prefix, "_", tab_name)
        if (isTRUE(s$write_gene_lists) && nrow(tg))
            .write_tsv(tg, file.path(s$output_dir, "GENE_LISTS",
                                     paste0(prefix, "_significant_term_genes_", s$namespace, ".tsv")))
        if (nrow(tg)) term_genes[[tab_name]] <- data.frame(table = tab_name, tg)
    }
    gene_lists <- do.call(rbind, unname(term_genes))
    summary <- do.call(rbind, unname(summary_rows))
    .write_tsv(summary, file.path(s$output_dir, paste0("summary_", s$namespace, ".tsv")))
    .log("Done: ", length(results), " result tables")

    res <- methods::new("GOFERResult", results = results, summary = summary,
                 category_genes = if (is.null(gene_lists)) data.frame() else gene_lists,
                 settings = s,
                 info = list(go_file = run$go_file, go_version = run$obo$data_version,
                             term_gene_counts = .term_gene_counts(ann),
                             n_annotated_genes = length(ann$genes),
                             annotation_stats = ann$stats, mode = "scores",
                             output_dir = normalizePath(s$output_dir)))
    .post_process(res)
}
