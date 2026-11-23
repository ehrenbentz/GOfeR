test_that("MWU p-values match wilcox.test with the normal approximation", {
    set.seed(1)
    N <- 400
    values <- round(rnorm(N), 1)
    memb <- Matrix::sparseMatrix(i = rep(1:3, c(20, 50, 7)),
                                 j = c(sample(N, 20), sample(N, 50), sample(N, 7)),
                                 x = 1, dims = c(3, N))
    res <- GOFER:::.mwu_tests(memb, values, c("greater", "less", "two.sided"))
    for (k in 1:3) {
        inside <- as.logical(memb[k, ])
        for (alt in c("greater", "less", "two.sided")) {
            w <- suppressWarnings(stats::wilcox.test(values[inside], values[!inside], alternative = alt,
                                                     exact = FALSE, correct = TRUE))
            expect_equal(res[[alt]][k], w$p.value, tolerance = 1e-10)
        }
    }
})

test_that("Fisher p-values match fisher.test", {
    set.seed(2)
    N <- 300
    sig <- runif(N) < 0.15
    memb <- Matrix::sparseMatrix(i = rep(1:2, c(30, 12)), j = c(sample(N, 30), sample(N, 12)),
                                 x = 1, dims = c(2, N))
    res <- GOFER:::.fisher_tests(memb, sig)
    for (k in 1:2) {
        inside <- as.logical(memb[k, ])
        ft <- stats::fisher.test(table(factor(inside, c(TRUE, FALSE)), factor(sig, c(TRUE, FALSE))),
                                 alternative = "greater")
        expect_equal(res$pval[k], ft$p.value, tolerance = 1e-10)
    }
})

test_that("average precision of random sets matches the direct computation", {
    N <- 500L
    sizes <- c(3L, 10L, 40L)
    col_of <- match(seq_len(max(sizes)), sizes)
    col_of[is.na(col_of)] <- 0L
    out <- GOFER:::.ap_null_chunk(5L, N, as.integer(col_of), length(sizes), 11)
    sets <- GOFER:::.random_sets(5L, N, max(sizes), 11)
    expect_true(all(apply(sets, 1, function(x) !anyDuplicated(x) && all(x >= 1 & x <= N))))
    for (b in 1:5) {
        for (j in seq_along(sizes)) {
            r <- sort(sets[b, seq_len(sizes[j])])
            expect_equal(out[b, j], mean(seq_along(r) / r), tolerance = 1e-12)
        }
    }
    expect_identical(GOFER:::.random_sets(3L, N, 10L, 7), GOFER:::.random_sets(3L, N, 10L, 7))
})

test_that("random ranks are uniform", {
    sets <- GOFER:::.random_sets(20000L, 10L, 1L, 3)
    expect_gt(stats::chisq.test(table(sets[, 1]))$p.value, 1e-4)
})

test_that("exceedance counts and kept tail values are correct", {
    set.seed(5)
    null_ap <- matrix(runif(200 * 3), 200, 3)
    ap_obs <- c(0.5, 0.99, 0.1, 0.7)
    size_col <- c(1L, 1L, 2L, 3L)
    counts <- numeric(4)
    top <- matrix(-Inf, 20, 3)
    GOFER:::.pr_count_chunk(null_ap, ap_obs, size_col, counts, top)
    expect_equal(counts, vapply(1:4, function(i) sum(null_ap[, size_col[i]] >= ap_obs[i]), numeric(1)))
    expect_equal(top[, 2], utils::tail(sort(null_ap[, 2]), 20))
})

test_that("average precision of categories is computed on the ranking", {
    values <- c(5, 4, 3, 2, 1)
    memb <- Matrix::sparseMatrix(i = c(1, 1, 2), j = c(1, 3, 5), x = 1, dims = c(2, 5))
    ap <- GOFER:::.average_precision(memb, values)
    expect_equal(ap, c((1 / 1 + 2 / 3) / 2, 1 / 5))
})
