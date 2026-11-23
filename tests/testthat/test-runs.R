ext <- system.file("extdata", package = "GOFER")
go  <- file.path(ext, "example_go.obo")
ann <- file.path(ext, "example_annotations.tsv")

test_that("config file values are read and arguments take precedence", {
    f <- tempfile(fileext = ".txt")
    write_config(f)
    expect_equal(GOFER:::.read_config(f), GOFER:::.default_settings()[names(GOFER:::.read_config(f))])
    writeLines(c("min_genes = 7", "methods = mwu, fisher", "retain_go_obo = TRUE", "max_level = NULL"), f)
    r <- GOFER:::.resolve_settings(list(min_genes = 9), f)
    expect_equal(r$settings$min_genes, 9L)
    expect_equal(r$settings$methods, c("mwu", "fisher"))
    expect_true(r$settings$retain_go_obo)
    expect_equal(r$origin[["methods"]], "config file")
    expect_equal(r$origin[["min_genes"]], "function argument")
    expect_error(GOFER:::.resolve_settings(list(min_gene = 3), NULL), "Unknown option")
    expect_error(GOFER:::.resolve_settings(list(namespace = "XX"), NULL), "not one of")
})

test_that("score run writes results and settings, and settings reproduce the run", {
    out <- file.path(tempdir(), "test_scores")
    res <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                         methods = "mwu,fisher,pr", pr_resamples = 2000, pr_max_resamples = 20000,
                         pr_tail_size = 50, cpus = 1, output_dir = out)
    expect_s4_class(res, "GOFERResult")
    expect_true(file.exists(file.path(out, "settings.txt")))
    expect_true(file.exists(file.path(out, "MWU", "UP", "example_deseq2_mwu_up_BP.tsv")))
    mw <- results(res, "example_deseq2", "mwu")
    expect_setequal(unique(mw$family), c("example_deseq2|mwu|up", "example_deseq2|mwu|down"))
    expect_true(all(mw$fdr >= mw$p_value))
    expect_equal(settings(res)$score_type, "wald")
    tg <- category_genes(res)
    expect_equal(names(tg), c("table", "term_id", "genes"))
    expect_true(file.exists(file.path(out, "GENE_LISTS", "example_deseq2_significant_term_genes_BP.tsv")))
    en <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                        test_direction = "enriched", cpus = 1, output_dir = file.path(tempdir(), "test_enriched"))
    expect_equal(settings(en)$score_type, "pvalue")
    expect_error(enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                               score_type = "pvalue", cpus = 1, output_dir = file.path(tempdir(), "test_err")),
                 "unsigned")
    res2 <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann,
                          config = file.path(out, "settings.txt"),
                          output_dir = file.path(tempdir(), "test_scores_repeat"))
    expect_equal(results(res2, "example_deseq2", "pr"), results(res, "example_deseq2", "pr"))
})

test_that("precision-recall results do not depend on the number of workers", {
    skip_on_os("windows")
    a <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                       score_type = "wald", methods = "pr", pr_resamples = 3000, pr_max_resamples = 30000, pr_tail_size = 50,
                       cpus = 1, output_dir = file.path(tempdir(), "test_pr1"))
    b <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                       score_type = "wald", methods = "pr", pr_resamples = 3000, pr_max_resamples = 30000, pr_tail_size = 50,
                       cpus = 2, output_dir = file.path(tempdir(), "test_pr2"))
    expect_equal(results(a, "example_deseq2", "pr"), results(b, "example_deseq2", "pr"))
})

test_that("module default method is signed kME MWU", {
    res <- enrich_modules(file.path(ext, "example_kme.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                          cpus = 1, output_dir = file.path(tempdir(), "test_modules_default"))
    expect_equal(settings(res)$module_methods, "signed_kme_mwu")
    expect_equal(results(res), paste0("example_kme|", c("darkolivegreen", "purple", "blue"), "|signed_kme_mwu"))
})

test_that("module run tests every module except grey and writes recurrence tables", {
    out <- file.path(tempdir(), "test_modules")
    res <- enrich_modules(file.path(ext, "example_kme.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                          module_methods = "fisher,signed_kme_mwu,unsigned_kme_mwu",
                          cpus = 1, output_dir = out)
    expect_length(results(res), 9)
    expect_true(file.exists(file.path(out, "MODULES", "example_kme", "example_kme_recurrence_BP.tsv")))
    expect_true(file.exists(file.path(out, "MODULES", "example_kme", "FISHER", "example_kme_blue_fisher_BP.tsv")))
    expect_true(all(c("fisher_n_modules", "signed_kme_mwu_n_modules") %in% names(recurrence(res, "example_kme"))))
})

test_that("unnamed lists, tiny universes and runs without significant genes complete", {
    tab <- GOFER:::.read_score_file(file.path(ext, "example_deseq2.tsv"))
    res <- enrich_scores(list(tab, tab), ann, go_file = go, retain_go_obo = TRUE, cpus = 1,
                         output_dir = file.path(tempdir(), "test_unnamed"))
    expect_true(all(c("table_1|mwu|up", "table_2|mwu|up") %in% results(res)))
    small <- enrich_scores(tab[1:30, ], ann, go_file = go, retain_go_obo = TRUE, cpus = 1,
                           methods = "mwu,fisher,pr", pr_resamples = 2000, pr_max_resamples = 20000, pr_tail_size = 50,
                           output_dir = file.path(tempdir(), "test_small"))
    expect_true(all(summary(small)$n_tests == 0))
    none <- enrich_scores(tab, ann, go_file = go, retain_go_obo = TRUE, cpus = 1, methods = "fisher",
                          threshold = 1e-300, output_dir = file.path(tempdir(), "test_none"))
    expect_equal(nrow(results(none, "scores", "fisher", "enriched")), 0)
})

test_that("adaptive resampling extends only category sizes with few exceedances", {
    res <- enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                         methods = "pr", pr_resamples = 2000, pr_max_resamples = 20000, pr_tail_size = 50,
                         cpus = 1, output_dir = file.path(tempdir(), "test_adaptive"))
    pr <- results(res, "example_deseq2", "pr")
    expect_true(all(pr$null_resamples %in% c(2000, 20000)))
    extended <- pr$null_resamples == 20000
    expect_true(any(extended))
    expect_true(all(pr$null_exceedances[!extended] >= 10))
    emp <- pr$p_value_method == "empirical"
    expect_equal(pr$p_value[emp], (pr$null_exceedances[emp] + 1) / (pr$null_resamples[emp] + 1))
})

test_that("pr_seed NULL follows set.seed() and the seed is recorded", {
    run <- function() enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                                    methods = "pr", pr_seed = NULL, pr_resamples = 2000, pr_max_resamples = 2000,
                                    pr_tail_size = 50, cpus = 1, output_dir = file.path(tempdir(), "test_seed"))
    set.seed(42); a <- run()
    set.seed(42); b <- run()
    expect_equal(settings(a)$pr_seed, settings(b)$pr_seed)
    expect_equal(results(a, "example_deseq2", "pr"), results(b, "example_deseq2", "pr"))
    expect_true(any(grepl(paste0("^pr_seed +=  *", settings(a)$pr_seed, "$"),
                          readLines(file.path(tempdir(), "test_seed", "settings.txt")))))
})
