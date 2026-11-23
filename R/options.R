## Registry of all options. The configuration file template, the settings
## file written with every run and argument checking are generated from this
## table.

.opt <- function(key, type, default, section, applies, help, choices = NULL) {
    list(key = key, type = type, default = default, section = section,
         applies = applies, help = help, choices = choices)
}

.options <- function() {
    list(
        ## ---------------- Basic ----------------
        .opt("output_dir", "character", "GOFER_results", "basic", "both",
             "Output directory, created inside the working directory unless an absolute path is given"),
        .opt("output_prefix", "character", NULL, "basic", "both",
             "Prefix of result file names. NULL: the name of each score table or kME matrix"),
        .opt("namespace", "character", "BP", "basic", "both",
             "GO namespace: BP, MF or CC", c("BP", "MF", "CC")),
        .opt("methods", "character_vector", "mwu", "basic", "scores",
             "Tests for score tables, comma-separated: mwu (Mann-Whitney U on ranks), fisher (Fisher's exact test on genes with p <= threshold), pr (average precision with a resampled null)",
             c("mwu", "fisher", "pr")),
        .opt("module_methods", "character_vector", "signed_kme_mwu", "basic", "modules",
             "Tests for modules, comma-separated: fisher (module membership), signed_kme_mwu, signed_kme_pr, unsigned_kme_mwu, unsigned_kme_pr. Signed kME tests are for signed and signed hybrid networks; unsigned kME (|kME|) tests are valid only for unsigned networks, because they rank genes that follow the inverted eigengene together with genes that follow the eigengene",
             c("fisher", "signed_kme_mwu", "unsigned_kme_mwu", "signed_kme_pr", "unsigned_kme_pr")),
        .opt("test_direction", "character", "both", "basic", "scores",
             "mwu and pr direction. both: UP and DOWN, each corrected separately. up / down: one of the two. enriched: one test toward high unsigned scores, regardless of direction",
             c("both", "up", "down", "enriched")),
        .opt("score_type", "character", NULL, "basic", "scores",
             "Gene score for mwu and pr. NULL: wald for both, up and down; pvalue for enriched. wald: Wald statistic (signed). signed_pvalue: -log10(p) * sign(log2FC). log2foldchange. pval_signedlog2fc: -log10(p) * log2FC. pvalue: -log10(p), unsigned. raw: score_column as given (signed)",
             c("wald", "signed_pvalue", "log2foldchange", "pval_signedlog2fc", "pvalue", "raw")),
        .opt("threshold", "numeric", 0.05, "basic", "scores",
             "fisher: genes with pvalue_column <= threshold form the significant set"),
        .opt("min_genes", "integer", 5L, "basic", "both",
             "Smallest GO category tested (number of genes in the universe)"),
        .opt("max_genes_fraction", "numeric", 0.1, "basic", "both",
             "Largest GO category tested, as a fraction of the annotated scored genes. NULL: no limit"),
        .opt("cpus", "integer", NULL, "basic", "both",
             "Worker processes for parallel steps. NULL: all physical cores"),
        .opt("report_fdr", "numeric", 0.05, "basic", "both",
             "FDR used for the summary table and the module recurrence table (results are not filtered)"),
        ## ---------------- Advanced ----------------
        .opt("pvalue_column", "character", "pvalue", "advanced", "scores",
             "P-value column. padj ranks or thresholds on adjusted p-values (genes with NA padj leave the universe)"),
        .opt("fold_change_column", "character", "log2FoldChange", "advanced", "scores",
             "Fold-change column"),
        .opt("wald_column", "character", "stat", "advanced", "scores",
             "Wald statistic column (score_type wald)"),
        .opt("score_column", "character", NULL, "advanced", "scores",
             "Score column for score_type raw"),
        .opt("respect_independent_filtering", "logical", FALSE, "advanced", "scores",
             "Drop genes with NA in any pvalue or padj column (DESeq2 independent filtering)"),
        .opt("two_sided", "logical", FALSE, "advanced", "scores",
             "mwu: one two-sided test per category instead of one-sided tests"),
        .opt("universe", "character", "annotated", "advanced", "both",
             "Background genes. annotated: all scored genes with a GO annotation in the namespace. retained: scored genes in at least one tested category",
             c("annotated", "retained")),
        .opt("propagate_part_of", "logical", TRUE, "advanced", "both",
             "Propagate annotations through part_of relationships in addition to is_a"),
        .opt("go_file", "character", NULL, "advanced", "both",
             "GO ontology in OBO format. NULL: go.obo in the user cache directory of the package. Downloaded if missing and replaced by the current release when older than 30 days (see retain_go_obo)"),
        .opt("retain_go_obo", "logical", FALSE, "advanced", "both",
             "TRUE: use an existing go_file whatever its age"),
        .opt("min_level", "integer", NULL, "advanced", "both",
             "Lowest GO level tested (longest path from the root; root = 0)"),
        .opt("max_level", "integer", NULL, "advanced", "both",
             "Highest GO level tested"),
        .opt("max_genes", "integer", NULL, "advanced", "both",
             "Largest GO category tested, as a number of genes. With max_genes_fraction the smaller limit applies"),
        .opt("reduce_redundancy", "logical", FALSE, "advanced", "both",
             "Merge categories that share most of their genes (complete linkage on 1 - shared / larger)"),
        .opt("cluster_cut_height", "numeric", 0.25, "advanced", "both",
             "Dissimilarity at or below which categories are merged (reduce_redundancy TRUE)"),
        .opt("pr_resamples", "integer", 100000L, "advanced", "both",
             "Random gene sets per category size for the average-precision null (first stage)"),
        .opt("pr_max_resamples", "integer", 1000000L, "advanced", "both",
             "Adaptive resampling: category sizes with an observed value below pr_min_exceedances get ten times as many random sets, repeated up to this number"),
        .opt("pr_seed", "integer", 1L, "advanced", "both",
             "Seed of the random gene sets; the same seed gives the same p-values on every system and for any cpus or backend. NULL: a seed drawn from R's random number generator, so set.seed() before the run fixes it; the seed used is written to settings.txt"),
        .opt("pr_min_exceedances", "integer", 10L, "advanced", "both",
             "Observed values with fewer random sets at or above them are resampled further (pr_max_resamples); beyond that the p-value comes from the generalized Pareto tail"),
        .opt("pr_tail_size", "integer", 250L, "advanced", "both",
             "Number of largest null values per category size used to fit the generalized Pareto tail"),
        .opt("reduce_methods", "character_vector", "ic", "basic", "both",
             "Term reduction: group the significant terms of every result table into themes, comma-separated methods. wang: semantic similarity from the GO graph. ic: semantic similarity from the information content of the annotations. slim: GO slim mapping. NULL: no reduction",
             c("wang", "ic", "slim")),
        .opt("plots", "character_vector", NULL, "basic", "both",
             "Figures written to the output directory (PNG), comma-separated: dots (dot plot of the top terms per test), dendrogram (significant terms clustered by shared genes), graph (significant terms in the GO hierarchy; needs the Rgraphviz package), theme_dots (dot plot of the themes from reduce_methods), theme_graph (themes in the GO hierarchy; needs Rgraphviz). NULL: none",
             c("dots", "dendrogram", "graph", "theme_dots", "theme_graph")),
        .opt("similarity_threshold", "numeric", 0.4, "advanced", "both",
             "reduce_methods wang and ic: terms are grouped by complete linkage so that every pair in a theme has at least this similarity (0 to 1)"),
        .opt("ic_measure", "character", "rel", "advanced", "both",
             "reduce_methods ic: rel (Schlicker relevance) or lin (Lin)", c("rel", "lin")),
        .opt("slim_subset", "character", "goslim_generic", "advanced", "both",
             "reduce_methods slim: GO subset defined in the GO file (subsetdef lines), for example goslim_generic, goslim_pir, goslim_drosophila"),
        .opt("dots_top_terms", "integer", 10L, "advanced", "both",
             "plots dots and theme_dots: most significant terms (themes) shown per test (reduced automatically when more than dots_max_terms rows result)"),
        .opt("dots_max_terms", "integer", 60L, "advanced", "both",
             "plots dots and theme_dots: largest number of rows in one dot plot"),
        .opt("dendrogram_top_terms", "integer", 40L, "advanced", "both",
             "plots dendrogram: largest number of significant terms shown"),
        .opt("graph_top_terms", "integer", 10L, "advanced", "both",
             "plots graph and theme_graph: most significant terms (themes) shown with all their ancestors"),
        .opt("graph_edge_labels", "logical", TRUE, "advanced", "both",
             "plots graph: print the relation (is_a, part_of) on every edge"),
        .opt("graph_term_ids", "logical", TRUE, "advanced", "both",
             "plots graph: print the GO ID in every box"),
        .opt("color_up", "character", "#B2182B", "advanced", "both",
             "Figure colour of UP results (graph fill ramps from white to this colour with significance; dendrogram labels)"),
        .opt("color_down", "character", "#2166AC", "advanced", "both",
             "Figure colour of DOWN results"),
        .opt("color_enriched", "character", "#D7301F", "advanced", "both",
             "Figure colour of results without direction (enriched, two-sided, Fisher, modules)"),
        .opt("color_is_a", "character", "#555555", "advanced", "both",
             "plots graph: colour of is_a edges (solid)"),
        .opt("color_part_of", "character", "#2166AC", "advanced", "both",
             "plots graph: colour of part_of edges (dashed)"),
        .opt("dots_colors", "character_vector", "viridis", "advanced", "both",
             "plots dots: colour scale from low to high significance; the name of an R colour palette (see grDevices::hcl.pals(), for example viridis, inferno, Reds) or two or more comma-separated colours (names or #RRGGBB)"),
        .opt("plot_dpi", "integer", 600L, "advanced", "both",
             "Resolution of the PNG figures"),
        .opt("module_column", "character", "module_color", "advanced", "modules",
             "kME matrix column holding the module assignment of each gene"),
        .opt("kme_prefix", "character", "kME", "advanced", "modules",
             "Prefix of the kME columns; the module name follows the prefix"),
        .opt("exclude_modules", "character_vector", "grey", "advanced", "modules",
             "Modules not tested, comma-separated (exact names). NULL: test all"),
        .opt("use_cache", "logical", TRUE, "advanced", "both",
             "Keep the parsed ontology and propagated annotations in the user cache directory of the package, so later runs with the same files skip parsing"),
        .opt("write_gene_lists", "logical", TRUE, "advanced", "both",
             "Write, for every term significant in a result table, all universe genes annotated to it (term, semicolon-separated genes)"),
        .opt("backend", "character", "auto", "advanced", "both",
             "Parallel backend. auto: forked workers on Linux and macOS, separate R sessions on Windows. multicore, snow, serial",
             c("auto", "multicore", "snow", "serial"))
    )
}

.option_keys <- function() vapply(.options(), `[[`, character(1), "key")

.option <- function(key) {
    opts <- .options()
    opts[[match(key, .option_keys())]]
}

## Text form of a value, as written in config and settings files.
.format_value <- function(value) {
    if (is.null(value)) return("NULL")
    if (is.logical(value)) return(paste(ifelse(value, "TRUE", "FALSE"), collapse = ", "))
    paste(as.character(value), collapse = ", ")
}

## Convert and check one option value, given as text (config file) or as an
## R value (function argument).
.coerce_option <- function(key, value, source) {
    opt <- .option(key)
    where <- paste0("Option '", key, "' (", source, ")")
    if (is.null(value)) return(NULL)
    if (is.character(value) && length(value) == 1 && toupper(trimws(value)) == "NULL") return(NULL)
    out <- switch(opt$type,
        logical = {
            v <- if (is.character(value)) toupper(trimws(value)) else value
            if (identical(v, "TRUE") || isTRUE(v)) TRUE
            else if (identical(v, "FALSE") || identical(v, FALSE)) FALSE
            else stop(where, ": expected TRUE or FALSE, got '", value, "'", call. = FALSE)
        },
        integer = {
            v <- suppressWarnings(as.numeric(value))
            if (length(v) != 1 || is.na(v) || v != round(v))
                stop(where, ": expected an integer, got '", value, "'", call. = FALSE)
            as.integer(v)
        },
        numeric = {
            v <- suppressWarnings(as.numeric(value))
            if (length(v) != 1 || is.na(v))
                stop(where, ": expected a number, got '", value, "'", call. = FALSE)
            v
        },
        character = {
            if (length(value) != 1) stop(where, ": expected one value", call. = FALSE)
            trimws(as.character(value))
        },
        character_vector = {
            v <- trimws(unlist(strsplit(as.character(value), ",")))
            v[nzchar(v)]
        })
    if (!is.null(opt$choices)) {
        bad <- setdiff(out, opt$choices)
        if (length(bad) > 0)
            stop(where, ": '", paste(bad, collapse = ", "), "' is not one of ",
                 paste(opt$choices, collapse = ", "), call. = FALSE)
    }
    out
}

.default_settings <- function() {
    stats::setNames(lapply(.options(), `[[`, "default"), .option_keys())
}
