## Precision-recall: p-values of average precision from random gene sets.
## The null distribution depends only on the universe size and the category
## size, so one null serves every family with the same universe size.

## Runs on one worker: its chunks of random sets, one at a time. Each chunk
## has its own seed, so results do not depend on the number of workers.
.pr_worker <- function(chunk_ids, n_universe, col_of, n_sizes, chunk_reps, seed,
                       ap_obs, size_col, n_keep) {
    counts <- numeric(length(ap_obs))
    top    <- matrix(-Inf, n_keep, n_sizes)
    for (cid in chunk_ids) {
        null_ap <- .ap_null_chunk(chunk_reps[cid], n_universe, col_of, n_sizes, seed + cid)
        .pr_count_chunk(null_ap, ap_obs, size_col, counts, top)
    }
    list(counts = counts, top = top)
}

## Random sets for the category sizes in `sizes_needed` (indices into
## size_values): exceedance counts for the observed values of those sizes
## and the n_keep largest null values per size.
.pr_stage <- function(B, sizes_needed, size_values, size_col, ap_obs, n_universe, seed,
                      n_keep, bp, chunk) {
    stage_values <- size_values[sizes_needed]
    obs <- which(size_col %in% sizes_needed)
    col_of <- match(seq_len(max(stage_values)), stage_values)
    col_of[is.na(col_of)] <- 0L
    n_chunks   <- ceiling(B / chunk)
    chunk_reps <- as.integer(pmin(chunk, B - (seq_len(n_chunks) - 1) * chunk))
    n_workers  <- max(1L, BiocParallel::bpnworkers(bp))
    groups     <- unname(split(seq_len(n_chunks), (seq_len(n_chunks) - 1) %% n_workers))
    parts <- BiocParallel::bplapply(groups, .pr_worker, n_universe = as.integer(n_universe),
                                    col_of = as.integer(col_of), n_sizes = length(stage_values),
                                    chunk_reps = chunk_reps, seed = seed,
                                    ap_obs = as.numeric(ap_obs[obs]),
                                    size_col = as.integer(match(size_col[obs], sizes_needed)),
                                    n_keep = as.integer(n_keep), BPPARAM = bp)
    counts <- Reduce(`+`, lapply(parts, `[[`, "counts"))
    top <- parts[[1]]$top
    for (p in parts[-1]) {
        for (k in seq_along(stage_values)) top[, k] <- utils::tail(sort(c(top[, k], p$top[, k])), n_keep)
    }
    list(obs = obs, counts = counts, top = top)
}

## Adaptive resampling: pr_resamples random sets for every category size;
## then, for the sizes that have an observed value with fewer than
## pr_min_exceedances random sets at or above it, ten times as many random
## sets in total, repeated until every observed value reaches
## pr_min_exceedances or the size reaches pr_max_resamples. Each stage and
## chunk has its own seed, so results do not depend on the number of
## workers.
.pr_null_summary <- function(ap_obs, sizes, n_universe, s, bp, chunk = 1000L) {
    n_keep      <- s$pr_tail_size + 1L
    size_values <- sort(unique(sizes))
    size_col    <- match(sizes, size_values)
    counts      <- numeric(length(ap_obs))
    top         <- matrix(-Inf, n_keep, length(size_values))
    B_size      <- rep(0, length(size_values))
    stage       <- 0L
    repeat {
        if (stage == 0L) {
            needed <- seq_along(size_values)
            target <- s$pr_resamples
        } else {
            low    <- counts < s$pr_min_exceedances & B_size[size_col] < s$pr_max_resamples
            needed <- sort(unique(size_col[low]))
            if (!length(needed)) break
            target <- min(B_size[needed[1]] * 10, s$pr_max_resamples)
        }
        extra <- target - B_size[needed[1]]
        st <- .pr_stage(extra, needed, size_values, size_col, ap_obs, n_universe,
                        seed = s$pr_seed * 1000003 + stage * 1e7, n_keep = n_keep, bp = bp, chunk = chunk)
        counts[st$obs] <- counts[st$obs] + st$counts
        for (k in seq_along(needed)) {
            j <- needed[k]
            top[, j] <- utils::tail(sort(c(top[, j], st$top[, k])), n_keep)
        }
        B_size[needed] <- target
        if (stage > 0L) {
            .log("  adaptive resampling: ", length(needed), " category sizes extended to ",
                 format(target, big.mark = ",", scientific = FALSE), " random sets")
        }
        stage <- stage + 1L
    }
    list(counts = counts, top = top, size_col = size_col, n_resample = B_size[size_col])
}

## Maximum-likelihood fit of a generalized Pareto distribution to
## exceedances y > 0; starting values from probability-weighted moments,
## moved into the feasible region 1 + xi * max(y) / sigma > 0 if needed.
.gpd_fit <- function(y) {
    n  <- length(y)
    ys <- sort(y)
    a0 <- mean(ys)
    a1 <- mean(ys * (n - seq_len(n)) / (n - 1))
    sigma0 <- 2 * a0 * a1 / (a0 - 2 * a1)
    xi0    <- -(a0 / (a0 - 2 * a1) - 2)
    if (!is.finite(sigma0) || sigma0 <= 0) {
        sigma0 <- a0
        xi0 <- 0
    }
    xi0 <- max(xi0, -0.99 * sigma0 / max(ys))
    nll <- function(par) {
        sigma <- exp(par[1])
        xi <- par[2]
        if (abs(xi) < 1e-8) return(n * log(sigma) + sum(ys) / sigma)
        z <- 1 + xi * ys / sigma
        if (any(z <= 0)) return(1e300)
        n * log(sigma) + (1 / xi + 1) * sum(log(z))
    }
    o <- stats::optim(c(log(sigma0), xi0), nll, method = "Nelder-Mead",
                      control = list(maxit = 2000, reltol = 1e-12))
    list(sigma = exp(o$par[1]), xi = o$par[2])
}

.gpd_survival <- function(z, fit) {
    if (abs(fit$xi) < 1e-8) return(exp(-z / fit$sigma))
    base <- 1 + fit$xi * z / fit$sigma
    ifelse(base > 0, base^(-1 / fit$xi), 0)
}

## p = (M + 1) / (B + 1), M = random sets at or above the observed value.
## When M < pr_min_exceedances, a generalized Pareto distribution is fitted
## to the pr_tail_size largest null values of that category size above a
## threshold midway between the pr_tail_size-th and the next largest value
## (Knijnenburg et al. 2009, Bioinformatics 25:i161), and
## p = (pr_tail_size / B) * survival above the threshold. An observed value
## beyond the upper end of the fitted distribution keeps (M + 1) / (B + 1).
.pr_pvalues <- function(ap_obs, summ, s) {
    B <- summ$n_resample
    M <- summ$counts
    n_tail <- s$pr_tail_size
    p <- (M + 1) / (B + 1)
    method <- rep("empirical", length(ap_obs))
    tail_cats <- which(M < s$pr_min_exceedances)
    for (k in unique(summ$size_col[tail_cats])) {
        idx <- tail_cats[summ$size_col[tail_cats] == k]
        top <- summ$top[, k]
        nk  <- length(top)
        thr <- (top[nk - n_tail] + top[nk - n_tail + 1]) / 2
        fit <- .gpd_fit(top[(nk - n_tail + 1):nk] - thr)
        sv  <- .gpd_survival(ap_obs[idx] - thr, fit)
        inside <- sv > 0
        p[idx[inside]] <- (n_tail / B[idx[inside]]) * sv[inside]
        method[idx[inside]]  <- "gpd_tail"
        method[idx[!inside]] <- "empirical_bound_beyond_gpd_endpoint"
    }
    list(p = p, method = method, exceedances = M, resamples = B)
}

## Precision-recall p-values for a list of families (each with ap, sizes,
## n_universe). One null per distinct universe size.
.pr_families <- function(families, s, bp) {
    if (methods::is(bp, "SnowParam") && BiocParallel::bpisup(bp)) {
        .log("Workers: ", BiocParallel::bpnworkers(bp), " (SnowParam), already running")
    } else if (!methods::is(bp, "SerialParam") && !BiocParallel::bpisup(bp)) {
        t0 <- proc.time()[["elapsed"]]
        BiocParallel::bpstart(bp)
        .log("Workers: ", BiocParallel::bpnworkers(bp), " (", class(bp)[1], "), started in ",
             round(proc.time()[["elapsed"]] - t0, 1), " s")
    }
    n_u <- vapply(families, `[[`, numeric(1), "n_universe")
    out <- vector("list", length(families))
    names(out) <- names(families)
    for (N in sort(unique(n_u))) {
        t0   <- proc.time()[["elapsed"]]
        keys <- names(families)[n_u == N]
        ap    <- unlist(lapply(families[keys], `[[`, "ap"), use.names = FALSE)
        sizes <- unlist(lapply(families[keys], `[[`, "sizes"), use.names = FALSE)
        ok    <- !is.na(ap)
        full  <- list(p = rep(NA_real_, length(ap)), method = rep(NA_character_, length(ap)),
                      exceedances = rep(NA_real_, length(ap)), resamples = rep(NA_real_, length(ap)))
        if (any(ok)) {
            summ <- .pr_null_summary(ap[ok], sizes[ok], N, s, bp)
            pv   <- .pr_pvalues(ap[ok], summ, s)
            full$p[ok] <- pv$p
            full$method[ok] <- pv$method
            full$exceedances[ok] <- pv$exceedances
            full$resamples[ok] <- pv$resamples
        }
        ends <- cumsum(vapply(families[keys], function(f) length(f$ap), numeric(1)))
        starts <- c(1, utils::head(ends, -1) + 1)
        for (i in seq_along(keys)) {
            idx <- seq_len(ends[i] - starts[i] + 1) + starts[i] - 1
            out[[keys[i]]] <- list(p = full$p[idx], method = full$method[idx],
                                   exceedances = full$exceedances[idx], resamples = full$resamples[idx])
        }
        .log("Precision-recall null, universe ", N, ": ", length(keys), " families, ",
             length(unique(sizes)), " category sizes, ",
             round(proc.time()[["elapsed"]] - t0, 1), " s")
    }
    out
}
