## kME matrices: a file, a directory of files, a data frame, or a named list
## of these. Each holds gene IDs in the first column, the module assignment
## in module_column and one kME column per module (kme_prefix + module).
.collect_kme_tables <- function(kme) {
    if (is.character(kme) && length(kme) == 1 && dir.exists(kme)) {
        files <- sort(list.files(kme, full.names = TRUE))
        files <- files[!startsWith(basename(files), ".") & !dir.exists(files)]
        kme <- stats::setNames(as.list(files), tools::file_path_sans_ext(basename(files)))
    } else if (is.character(kme)) {
        missing <- kme[!file.exists(kme)]
        if (length(missing)) stop("kME file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
        kme <- stats::setNames(as.list(kme), tools::file_path_sans_ext(basename(kme)))
    } else if (is.data.frame(kme)) {
        kme <- list(modules = kme)
    }
    if (is.null(names(kme))) names(kme) <- rep("", length(kme))
    unnamed <- !nzchar(names(kme))
    names(kme)[unnamed] <- paste0("network_", which(unnamed))
    lapply(kme, function(x) {
        if (is.character(x)) x <- utils::read.delim(x, check.names = FALSE, stringsAsFactors = FALSE)
        x <- as.data.frame(x, check.names = FALSE)
        names(x)[1] <- "gene"
        x$gene <- as.character(x$gene)
        x
    })
}

#' GO enrichment of WGCNA co-expression modules
#'
#' Tests GO categories for every module of one or more WGCNA networks. Each
#' module is its own multiple-testing family. Methods (option
#' \code{module_methods}): \code{fisher} (module membership),
#' \code{signed_kme_mwu} and \code{unsigned_kme_mwu} (Mann-Whitney U on the
#' signed or absolute kME of all genes, upward), \code{signed_kme_pr} and
#' \code{unsigned_kme_pr} (average precision on the same rankings). A
#' recurrence table per network counts the modules in which each term is
#' significant.
#'
#' @param kme A kME matrix file (gene IDs in the first column, a module
#'   column, one \code{kME<module>} column per module), a directory of such
#'   files (one network each), a data frame, or a named list of these.
#' @param annotations Gene-to-GO annotations (see \code{\link{enrich_scores}}).
#' @param ... Options as \code{key = value}; see \code{\link{write_config}}.
#' @param config Optional configuration file.
#' @return A \code{GOFERResult} object.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_modules(file.path(ext, "example_kme.tsv"),
#'                       file.path(ext, "example_annotations.tsv"),
#'                       go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                       cpus = 1, use_cache = FALSE, output_dir = file.path(tempdir(), "example_modules"))
#' res
#' head(results(res, "example_kme", "blue", "signed_kme_mwu"))
#' @export
enrich_modules <- function(kme, annotations, ..., config = NULL) {
    run <- .start_run(list(...), config, "modules")
    s <- run$s
    on.exit(.finish_run(run), add = TRUE)
    ann <- .cached_annotations(annotations, run)
    tables <- .collect_kme_tables(kme)
    meths <- s$module_methods
    results <- list(); pr_families <- list(); designs <- list(); module_info <- list()

    for (net in names(tables)) {
        tab <- tables[[net]]
        if (!s$module_column %in% names(tab))
            stop("Module column '", s$module_column, "' not found in '", net, "'", call. = FALSE)
        kme_cols <- names(tab)[startsWith(names(tab), s$kme_prefix)]
        modules  <- setdiff(substring(kme_cols, nchar(s$kme_prefix) + 1), s$exclude_modules)
        if (!length(modules)) stop("No kME columns with prefix '", s$kme_prefix, "' in '", net, "'", call. = FALSE)
        d <- .build_design(tab$gene, seq_len(nrow(tab)), ann, run$hier, s, net)
        rows <- d$values
        designs[[net]] <- d
        base <- .base_columns(d, run$obo)
        for (m in modules) {
            k <- tab[[paste0(s$kme_prefix, m)]][rows]
            if (anyNA(k)) stop("Missing kME values for module '", m, "' in '", net, "'", call. = FALSE)
            member <- tab[[s$module_column]][rows] == m
            key0 <- paste0(net, "|", m)
            module_info[[key0]] <- list(network = net, module = m, members = sum(member))
            mbase <- base
            mbase$module <- rep(m, nrow(base))
            mbase$module_members <- rep(sum(member), nrow(base))
            if ("fisher" %in% meths) {
                fr <- .fisher_tests(d$membership, member)
                df <- mbase
                df$n_module_genes <- fr$sig_genes
                df$enrichment_ratio <- fr$enrichment_ratio
                df$odds_ratio <- fr$odds_ratio
                df$p_value <- fr$pval
                results[[paste0(key0, "|fisher")]] <- .finish_family(df, d, s)
            }
            for (mm in intersect(c("signed_kme_mwu", "unsigned_kme_mwu"), meths)) {
                v <- if (mm == "signed_kme_mwu") k else abs(k)
                res <- .mwu_tests(d$membership, v, "greater")
                df <- mbase
                df$delta_rank <- res$delta_rank
                df$p_value <- res$greater
                results[[paste0(key0, "|", mm)]] <- .finish_family(df, d, s)
            }
            for (pm in intersect(c("signed_kme_pr", "unsigned_kme_pr"), meths)) {
                v <- if (pm == "signed_kme_pr") k else abs(k)
                pr_families[[paste0(key0, "|", pm)]] <- list(
                    ap = .average_precision(d$membership, v), sizes = d$sizes,
                    n_universe = d$n_universe, base = mbase, network = net)
            }
        }
    }

    if (length(pr_families)) {
        pv <- .pr_families(pr_families, s, run$bp)
        for (fk in names(pr_families)) {
            f <- pr_families[[fk]]
            df <- f$base
            df$average_precision <- f$ap
            df$null_exceedances <- pv[[fk]]$exceedances
            df$null_resamples <- pv[[fk]]$resamples
            df$p_value <- pv[[fk]]$p
            df$p_value_method <- pv[[fk]]$method
            results[[fk]] <- .finish_family(df, designs[[f$network]], s)
        }
    }

    summary_rows <- list(); sig_terms <- list()
    for (key in names(results)) {
        parts <- strsplit(key, "|", fixed = TRUE)[[1]]
        net <- parts[1]; m <- parts[2]; meth <- parts[3]
        df <- results[[key]]
        prefix <- if (is.null(s$output_prefix)) net else paste0(s$output_prefix, "_", net)
        .write_tsv(df, file.path(s$output_dir, "MODULES", net, toupper(meth),
                                 paste0(prefix, "_", m, "_", meth, "_", s$namespace, ".tsv")))
        sig <- !is.na(df$fdr) & df$fdr <= s$report_fdr
        summary_rows[[key]] <- data.frame(network = net, module = m,
                                          module_members = module_info[[paste0(net, "|", m)]]$members,
                                          method = meth, n_tests = nrow(df), n_universe = .first(df$n_universe),
                                          n_significant = sum(sig), fdr_cutoff = s$report_fdr)
        sig_terms[[net]] <- union(sig_terms[[net]], df$term_id[sig])
    }
    summary <- do.call(rbind, unname(summary_rows))
    .write_tsv(summary, file.path(s$output_dir, paste0("summary_modules_", s$namespace, ".tsv")))

    recurrence <- list()
    for (net in names(tables)) {
        rec <- .base_columns(designs[[net]], run$obo)
        mods <- unique(vapply(names(results)[startsWith(names(results), paste0(net, "|"))],
                              function(k) strsplit(k, "|", fixed = TRUE)[[1]][2], character(1)))
        for (meth in meths) {
            keys <- paste0(net, "|", mods, "|", meth)
            sig <- vapply(keys, function(k) {
                df <- results[[k]]
                q <- df$fdr[match(rec$term_id, df$term_id)]
                !is.na(q) & q <= s$report_fdr
            }, logical(nrow(rec)))
            sig <- matrix(sig, nrow = nrow(rec))
            rec[[paste0(meth, "_n_modules")]] <- rowSums(sig)
            hit <- which(sig, arr.ind = TRUE)
            hit_mods <- character(nrow(rec))
            if (nrow(hit)) {
                joined <- tapply(mods[hit[, 2]], hit[, 1], paste, collapse = ";")
                hit_mods[as.integer(names(joined))] <- joined
            }
            rec[[paste0(meth, "_modules")]] <- hit_mods
        }
        n_cols <- paste0(meths, "_n_modules")
        rec <- rec[rowSums(rec[, n_cols, drop = FALSE]) > 0, , drop = FALSE]
        rec <- rec[order(-rec[[n_cols[1]]], rec$term_id), , drop = FALSE]
        recurrence[[net]] <- rec
        prefix <- if (is.null(s$output_prefix)) net else paste0(s$output_prefix, "_", net)
        .write_tsv(rec, file.path(s$output_dir, "MODULES", net, paste0(prefix, "_recurrence_", s$namespace, ".tsv")))
    }

    term_genes <- list()
    for (net in names(sig_terms)) {
        tg <- .term_genes(designs[[net]], sig_terms[[net]])
        prefix <- if (is.null(s$output_prefix)) net else paste0(s$output_prefix, "_", net)
        if (isTRUE(s$write_gene_lists) && nrow(tg))
            .write_tsv(tg, file.path(s$output_dir, "GENE_LISTS",
                                     paste0(prefix, "_significant_term_genes_", s$namespace, ".tsv")))
        if (nrow(tg)) term_genes[[net]] <- data.frame(table = net, tg)
    }
    gene_lists <- do.call(rbind, unname(term_genes))
    .log("Done: ", length(module_info), " modules in ", length(tables), " network(s)")

    res <- methods::new("GOFERResult", results = results, summary = summary,
                 category_genes = if (is.null(gene_lists)) data.frame() else gene_lists,
                 settings = s,
                 info = list(go_file = run$go_file, go_version = run$obo$data_version,
                             term_gene_counts = .term_gene_counts(ann),
                             n_annotated_genes = length(ann$genes),
                             annotation_stats = ann$stats, mode = "modules",
                             recurrence = recurrence, output_dir = normalizePath(s$output_dir)))
    .post_process(res)
}
