## Enrichment tests on a design (categories x universe genes).

## Mann-Whitney U for all categories at once: rank sums from one sparse
## product with the global ranks; normal approximation with tie and
## continuity correction.
.mwu_tests <- function(memb, values, alternatives) {
    N  <- length(values)
    r  <- rank(values)
    n1 <- Matrix::rowSums(memb)
    n2 <- N - n1
    rank_sum   <- as.vector(memb %*% r)
    delta_rank <- rank_sum / n1 - (sum(r) - rank_sum) / n2
    u1     <- rank_sum - n1 * (n1 + 1) / 2
    mean_u <- n1 * n2 / 2
    counts <- tabulate(match(values, unique(values)))
    tie    <- sum(as.numeric(counts)^3 - counts)
    sd_u   <- sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))))
    bad    <- (n1 < 2) | (n2 < 2) | !(sd_u > 0)
    out <- list(delta_rank = delta_rank)
    for (alt in alternatives) {
        p <- switch(alt,
            greater   = stats::pnorm((u1 - mean_u - 0.5) / sd_u, lower.tail = FALSE),
            less      = stats::pnorm((u1 - mean_u + 0.5) / sd_u),
            two.sided = pmin(1, 2 * stats::pnorm((abs(u1 - mean_u) - 0.5) / sd_u, lower.tail = FALSE)))
        p[bad] <- NA
        out[[alt]] <- p
    }
    out
}

## Fisher's exact test (one-sided, over-representation) from the
## hypergeometric upper tail.
.fisher_tests <- function(memb, significant) {
    N     <- length(significant)
    n_sig <- sum(significant)
    size  <- Matrix::rowSums(memb)
    a <- round(as.vector(memb %*% as.numeric(significant)))
    b <- size - a
    c <- n_sig - a
    d <- N - n_sig - b
    list(sig_genes = a,
         odds_ratio = (a * d) / (b * c),
         enrichment_ratio = (a / size) / (n_sig / N),
         pval = stats::phyper(a - 1, n_sig, N - n_sig, size, lower.tail = FALSE))
}

## Average precision of each category in the ranking by decreasing value.
.average_precision <- function(memb, values) {
    ord <- order(-values)
    tm  <- methods::as(Matrix::t(memb[, ord, drop = FALSE]), "CsparseMatrix")
    counts <- diff(tm@p)
    ranks  <- tm@i + 1
    within <- sequence(counts)
    grp    <- rep(seq_along(counts), counts)
    out <- rep(NA_real_, length(counts))
    nz  <- counts > 0
    out[nz] <- as.vector(rowsum(within / ranks, grp, reorder = TRUE)) / counts[nz]
    out
}
