ext <- system.file("extdata", package = "GOFER")
go  <- file.path(ext, "example_go.obo")
ann <- file.path(ext, "example_annotations.tsv")

run <- function(dir, ...) {
    enrich_scores(file.path(ext, "example_deseq2.tsv"), ann, go_file = go, retain_go_obo = TRUE,
                  cpus = 1, report_fdr = 0.1, output_dir = file.path(tempdir(), dir), ...)
}

test_that("Wang similarity of a small graph matches the definition", {
    ## root <- a <- b (is_a), root <- c (is_a)
    hier <- list(level = c(0, 1, 2, 1),
                 edges = data.frame(child = c(2L, 3L, 4L), parent = c(1L, 2L, 1L), weight = 0.8),
                 ancestors = Matrix::sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4), j = c(1, 1, 2, 1, 2, 3, 1, 4),
                                                  x = 1, dims = c(4, 4)))
    sim <- GOFER:::.wang_similarity(c(2L, 3L, 4L), hier)
    ## S_a = {a 1, root 0.8}; S_b = {b 1, a 0.8, root 0.64}; S_c = {c 1, root 0.8}
    expect_equal(sim[1, 2], ((1 + 0.8) + (0.8 + 0.64)) / (1.8 + 2.44))
    expect_equal(sim[1, 3], (0.8 + 0.8) / (1.8 + 1.8))
    expect_equal(diag(sim), c(1, 1, 1))
})

test_that("information content similarity uses the most informative common ancestor", {
    anc <- Matrix::sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4), j = c(1, 1, 2, 1, 2, 3, 1, 4), x = 1, dims = c(4, 4))
    ic <- c(0, 1, 2, 1.5)
    sim <- GOFER:::.ic_similarity(c(2L, 3L, 4L), list(ancestors = anc), ic, "lin")
    expect_equal(sim[1, 2], 2 * 1 / (1 + 2))
    expect_equal(sim[1, 3], 0)
    rel <- GOFER:::.ic_similarity(c(2L, 3L, 4L), list(ancestors = anc), ic, "rel")
    expect_equal(rel[1, 2], 2 * 1 / (1 + 2) * (1 - exp(-1)))
})

test_that("term reduction assigns every significant term and writes tables", {
    res <- run("test_reduce", reduce_methods = "wang,ic,slim")
    down <- results(res, "example_deseq2", "mwu", "down")
    sig <- down$term_id[down$fdr <= 0.1]
    expect_gt(length(sig), 3)
    for (m in c("wang", "ic")) {
        mem <- theme_members(res, "example_deseq2", "mwu", "down", m)
        th  <- themes(res, "example_deseq2", "mwu", "down", m)
        expect_setequal(mem$term_id, sig)
        expect_false(anyDuplicated(mem$term_id) > 0)
        expect_equal(sum(th$n_terms), length(sig))
        expect_true(all(mem$similarity_to_medoid >= 0.4 - 1e-12))
    }
    expect_setequal(unique(theme_members(res, "example_deseq2", "mwu", "down", "slim")$term_id), sig)
    expect_true(file.exists(file.path(tempdir(), "test_reduce", "THEMES", "WANG",
                                      "example_deseq2_mwu_down_BP_wang_themes.tsv")))
    one <- reduce_terms(res, "wang", similarity_threshold = 0.99)
    expect_gte(nrow(themes(one, "example_deseq2", "mwu", "down", "wang")),
               nrow(themes(res, "example_deseq2", "mwu", "down", "wang")))
    expect_error(reduce_terms(res, "slim", slim_subset = "nope"), "Subsets in the file")
})

test_that("figures are written", {
    res <- run("test_plots", plots = "dots,dendrogram", plot_dpi = 72)
    fig <- file.path(tempdir(), "test_plots", "FIGURES")
    expect_true(file.exists(file.path(fig, "DOTS", "example_deseq2_BP_dots.png")))
    expect_true(file.exists(file.path(fig, "DENDROGRAM", "example_deseq2_mwu_BP_dendrogram.png")))
    skip_if_not_installed("Rgraphviz")
    f <- file.path(tempdir(), "graph.png")
    plot_go_graph(res, "example_deseq2", "mwu", file = f, dpi = 72)
    expect_gt(file.size(f), 1000)
    expect_error(plot_dots(res, "nope"), "No result table matches")
})

test_that("colour options are checked and used", {
    expect_error(run("test_bad_colour", color_up = "notacolour"), "is not a colour")
    expect_error(run("test_bad_palette", dots_colors = "notapalette"), "dots_colors")
    res <- run("test_colours", plots = "dots", plot_dpi = 72, color_up = "darkorange", dots_colors = "white,black")
    cols <- GOFER:::.plot_colours(settings(res))
    expect_equal(unname(cols$direction["up"]), "darkorange")
    expect_equal(cols$dots, c("white", "black"))
    f <- file.path(tempdir(), "dend_col.png")
    plot_dendrogram(res, "example_deseq2", "mwu", file = f, dpi = 72, colors = list(color_down = "#00FF00"))
    expect_true(file.exists(f))
})

test_that("theme figures are written", {
    res <- run("test_theme_plots", reduce_methods = "ic,slim", plots = "theme_dots", plot_dpi = 72)
    fig <- file.path(tempdir(), "test_theme_plots", "FIGURES", "THEME_DOTS")
    expect_true(file.exists(file.path(fig, "example_deseq2_BP_ic_theme_dots.png")))
    expect_true(file.exists(file.path(fig, "example_deseq2_BP_slim_theme_dots.png")))
    rows <- GOFER:::.theme_rows(res, results(res)[1:2], "ic")
    expect_true(all(rows$n_terms >= 1))
    expect_error(plot_theme_dots(res, "example_deseq2", method = "wang"), "No themes from method")
    skip_if_not_installed("Rgraphviz")
    f <- file.path(tempdir(), "theme_graph.png")
    plot_theme_graph(res, "example_deseq2", "mwu", method = "slim", file = f, dpi = 72)
    expect_gt(file.size(f), 1000)
})
