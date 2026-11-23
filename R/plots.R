## Figures. Every figure is written as PNG; without a file name it is drawn
## on the current graphics device.

.open_png <- function(file, width, height, dpi) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    if (.Platform$OS.type != "windows" && isTRUE(capabilities("cairo"))) {
        grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = "cairo")
    } else {
        grDevices::png(file, width = width, height = height, units = "in", res = dpi)
    }
}

## Result tables selected by the leading parts of their names; "*" matches
## any value of a part.
.match_keys <- function(x, parts) {
    keys <- names(x@results)
    if (!length(parts)) return(keys)
    split_keys <- strsplit(keys, "|", fixed = TRUE)
    ok <- vapply(split_keys, function(k) {
        length(k) >= length(parts) && all(parts == "*" | parts == k[seq_along(parts)])
    }, logical(1))
    if (!any(ok))
        stop("No result table matches '", paste(parts, collapse = "|"), "'. Available: ",
             paste(keys, collapse = ", "), call. = FALSE)
    keys[ok]
}

## Significant terms of the selected result tables in long format.
.significant_terms <- function(x, keys, fdr) {
    out <- lapply(keys, function(k) {
        df <- x@results[[k]]
        df <- df[!is.na(df$fdr) & df$fdr <= fdr, c("term_id", "term_name", "n_genes", "fdr"), drop = FALSE]
        if (!nrow(df)) return(NULL)
        cbind(family = k, df)
    })
    out <- do.call(rbind, out)
    if (is.null(out)) return(data.frame(family = character(0), term_id = character(0),
                                        term_name = character(0), n_genes = integer(0), fdr = numeric(0)))
    out
}

## Short labels of result tables: the parts of their names that differ.
.family_labels <- function(keys) {
    parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    varying <- apply(parts, 2, function(col) length(unique(col)) > 1)
    if (!any(varying)) varying[ncol(parts)] <- TRUE
    stats::setNames(apply(parts[, varying, drop = FALSE], 1, paste, collapse = " "), keys)
}

.direction_of <- function(keys, mode) {
    if (mode != "scores") return(rep("enriched", length(keys)))
    vapply(strsplit(keys, "|", fixed = TRUE), function(k) k[length(k)], character(1))
}

## Colours of a run's settings, checked.
.plot_colours <- function(s) {
    one <- function(key) {
        v <- s[[key]]
        ok <- tryCatch({ grDevices::col2rgb(v); TRUE }, error = function(e) FALSE)
        if (!ok) stop("Option '", key, "': '", v, "' is not a colour", call. = FALSE)
        v
    }
    dots <- s$dots_colors
    pal <- grDevices::hcl.pals()[match(tolower(dots), tolower(grDevices::hcl.pals()))]
    if (length(dots) == 1 && !is.na(pal)) {
        dots <- grDevices::hcl.colors(9, pal, rev = TRUE)
    } else {
        bad <- dots[!vapply(dots, function(v) !inherits(try(grDevices::col2rgb(v), silent = TRUE), "try-error"), logical(1))]
        if (length(dots) < 2 || length(bad))
            stop("Option 'dots_colors': give a palette name from grDevices::hcl.pals() or two or more colours",
                 call. = FALSE)
    }
    list(direction = c(up = one("color_up"), down = one("color_down"),
                       enriched = one("color_enriched"), two_sided = one("color_enriched")),
         edge = c(is_a = one("color_is_a"), part_of = one("color_part_of")),
         dots = dots)
}

.shorten <- function(text, width) {
    long <- nchar(text) > width
    text[long] <- paste0(substr(text[long], 1, width - 3), "...")
    text
}

.dots_plot <- function(x, keys, top_terms, max_terms, fdr, file, dpi, title = NULL,
                       cols = .plot_colours(x@settings)) {
    d <- .significant_terms(x, keys, fdr)
    if (!nrow(d)) return(invisible(NULL))
    .dots_draw(d, keys, top_terms, max_terms, file, dpi, title, cols, "Genes")
}

## Dot plot of rows (term_id, term_name) against result tables (family);
## dot size from n_genes, colour from -log10 fdr.
.dots_draw <- function(d, keys, top_terms, max_terms, file, dpi, title, cols, size_name) {
    n_top <- top_terms
    repeat {
        sel <- unique(unlist(lapply(split(d, d$family), function(f)
            utils::head(f$term_id[order(f$term_id == "other", f$fdr)], n_top))))
        if (length(sel) <= max_terms || n_top == 1) break
        n_top <- n_top - 1
    }
    sel <- utils::head(sel, max_terms)
    d <- d[d$term_id %in% sel, , drop = FALSE]
    labels <- .family_labels(keys)
    d$family_label <- factor(labels[d$family], levels = unique(labels))
    best <- do.call(rbind, lapply(split(d, d$term_id), function(t) {
        i <- which.min(t$fdr)
        data.frame(term_id = t$term_id[1], term_name = t$term_name[1],
                   first = min(as.integer(t$family_label)), fdr = t$fdr[i])
    }))
    best <- best[order(best$term_id == "other", best$first, best$fdr, best$term_id), ]
    best$label <- .shorten(best$term_name, 60)
    dup <- duplicated(best$label) | duplicated(best$label, fromLast = TRUE)
    best$label[dup] <- paste0(best$label[dup], " (", best$term_id[dup], ")")
    d$term_label <- factor(best$label[match(d$term_id, best$term_id)], levels = rev(best$label))
    d$score <- -log10(pmax(d$fdr, 1e-300))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$family_label, y = .data$term_label,
                                         size = .data$n_genes, colour = .data$score)) +
        ggplot2::geom_point() +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_size_area(max_size = 4.5, name = size_name) +
        ggplot2::scale_colour_gradientn(colours = cols$dots, name = expression(-log[10] ~ FDR)) +
        ggplot2::labs(x = NULL, y = NULL, title = title) +
        ggplot2::theme_bw(base_size = 8) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
                       axis.text.y = ggplot2::element_text(colour = "black"),
                       panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "grey90"),
                       plot.title = ggplot2::element_text(size = 9),
                       legend.key.height = ggplot2::unit(0.14, "in"),
                       legend.key.width = ggplot2::unit(0.12, "in"))
    if (!is.null(file)) {
        n_fam <- nlevels(d$family_label)
        width  <- min(49, 0.5 + 0.062 * max(nchar(best$label)) + max(1.2, 0.24 * n_fam) + 1.1)
        height <- min(49, 0.5 + 0.15 * nrow(best) + 0.05 * max(nchar(levels(d$family_label))) +
                          if (is.null(title)) 0 else 0.25)
        height <- max(height, 2.2)
        .open_png(file, width, height, dpi)
        print(p)
        grDevices::dev.off()
    } else {
        print(p)
    }
    invisible(p)
}

## Significant terms clustered by shared genes: average linkage on
## 1 - shared genes / genes of the larger category. Label colour: direction
## of the test; type face: significance.
.dendrogram_plot <- function(x, keys, top_terms, fdr, file, dpi, title = NULL,
                             cols = .plot_colours(x@settings)) {
    d <- .significant_terms(x, keys, fdr)
    if (!nrow(d)) return(invisible(NULL))
    d$direction <- .direction_of(d$family, x@info$mode)
    d <- d[order(d$fdr, d$term_id), ]
    d <- d[!duplicated(d$term_id), ]
    table <- strsplit(keys[1], "|", fixed = TRUE)[[1]][1]
    all_sets <- .gene_sets_by_table(x)[[table]]
    d <- utils::head(d[d$term_id %in% names(all_sets), , drop = FALSE], top_terms)
    if (nrow(d) < 3) return(invisible(NULL))
    sets <- all_sets[d$term_id]
    genes <- unique(unlist(sets))
    memb <- Matrix::sparseMatrix(i = rep(seq_along(sets), lengths(sets)),
                                 j = match(unlist(sets), genes), x = 1,
                                 dims = c(length(sets), length(genes)))
    shared <- as.matrix(Matrix::tcrossprod(memb))
    sizes <- diag(shared)
    dist <- 1 - shared / outer(sizes, sizes, pmax)
    hc <- stats::hclust(stats::as.dist(dist), method = "average")
    ord <- hc$order
    labels <- paste0(.shorten(d$term_name, 62), "  (", d$n_genes, ")")
    font <- ifelse(d$fdr <= 0.001, 2, ifelse(d$fdr <= 0.01, 1, 3))
    cex  <- ifelse(d$fdr <= 0.001, 0.80, ifelse(d$fdr <= 0.01, 0.74, 0.68))
    col  <- unname(cols$direction[d$direction])
    n <- nrow(d)
    width <- 3.2 + 0.068 * max(nchar(labels))
    height <- 1.0 + 0.17 * n
    if (!is.null(file)) .open_png(file, width, height, dpi)
    old <- graphics::par(mai = c(0.55, 0.15, if (is.null(title)) 0.1 else 0.3, width - 2.6), xpd = NA)
    on.exit({
        graphics::par(old)
        if (!is.null(file)) grDevices::dev.off()
    })
    graphics::plot(stats::as.dendrogram(hc), horiz = TRUE, leaflab = "none", axes = FALSE,
                   xlab = "", ylab = "", yaxs = "i", ylim = c(0.4, n + 0.6))
    graphics::axis(1, at = pretty(c(0, max(hc$height))), cex.axis = 0.6, tcl = -0.25, mgp = c(1.2, 0.3, 0), lwd = 0.6)
    graphics::mtext("1 - shared genes / genes of the larger term", side = 1, line = 1.3, cex = 0.6)
    label_x <- graphics::grconvertX(graphics::grconvertX(0, "user", "inches") + 0.05, "inches", "user")
    graphics::text(label_x, seq_len(n), labels[ord], pos = 4, offset = 0, col = col[ord],
                   font = font[ord], cex = cex[ord])
    if (!is.null(title)) graphics::mtext(title, side = 3, line = 0.4, cex = 0.75, adj = 0)
    used <- intersect(names(cols$direction), unique(d$direction))
    legend_text <- c(if (x@info$mode == "scores" && any(used %in% c("up", "down"))) used,
                     "FDR <= 0.001", "FDR <= 0.01", paste0("FDR <= ", fdr))
    n_dir <- length(legend_text) - 3
    graphics::legend(x = graphics::grconvertX(1, "nfc", "user"), y = graphics::grconvertY(0, "nfc", "user"),
                     xjust = 1, yjust = 0, legend = legend_text, bty = "n", cex = 0.55, horiz = TRUE,
                     text.col = c(unname(cols$direction[legend_text[seq_len(n_dir)]]), rep("black", 3)),
                     text.font = c(rep(1, n_dir), 2, 1, 3))
    invisible(hc)
}

.wrap_label <- function(text, width) {
    vapply(text, function(t) paste(strwrap(t, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
}

## Fill colour of a significant box: from a light tint of the direction
## colour (FDR at the cutoff) to the full colour (highest significance).
.fill_ramp <- function(score, direction, colours) {
    top <- max(score, 3)
    vapply(seq_along(score), function(i) {
        base <- colours[[direction[i]]]
        light <- grDevices::colorRampPalette(c("#FFFFFF", base))(7)[2]
        f <- grDevices::colorRamp(c(light, base))(min(1, score[i] / top))
        grDevices::rgb(f[1], f[2], f[3], maxColorValue = 255)
    }, character(1))
}

## The most significant terms with all their ancestors as a rooted graph.
## Graphviz (dot, through Rgraphviz) computes node positions and edge
## routes; boxes, edges and labels are then drawn with base graphics, so
## every edge runs from the bottom edge of the parent box to the top edge of
## the child box. Filled boxes are significant (colour by direction, depth
## by -log10 FDR); is_a edges solid grey, part_of edges dashed blue.
.graph_plot <- function(x, keys, top_terms, fdr, file, dpi, title = NULL,
                        edge_labels = TRUE, term_ids = TRUE, cols = .plot_colours(x@settings)) {
    if (!requireNamespace("Rgraphviz", quietly = TRUE) || !requireNamespace("graph", quietly = TRUE))
        stop("The GO graph needs the Bioconductor package Rgraphviz: BiocManager::install(\"Rgraphviz\")",
             call. = FALSE)
    d <- .significant_terms(x, keys, fdr)
    if (!nrow(d)) return(invisible(NULL))
    d$direction <- .direction_of(d$family, x@info$mode)
    d <- d[order(d$fdr, d$term_id), ]
    d <- d[!duplicated(d$term_id), ]
    d <- utils::head(d, top_terms)
    .graph_draw(x, d, file, dpi, title, edge_labels, term_ids, cols)
}

## GO graph of the terms in d (term_id, fdr, direction; optional extra: one
## more label line) with all their ancestors.
.graph_draw <- function(x, d, file, dpi, title, edge_labels, term_ids, cols) {
    ont <- .result_ontology(x)
    obo <- ont$obo; hier <- ont$hier
    idx <- match(d$term_id, obo$current)
    nodes <- which(Matrix::colSums(hier$ancestors[idx, , drop = FALSE]) > 0)
    e <- hier$edges[hier$edges$child %in% nodes & hier$edges$parent %in% nodes, , drop = FALSE]
    ids <- obo$current[nodes]
    g <- graph::graphNEL(nodes = ids, edgemode = "directed")
    if (nrow(e)) g <- graph::addEdge(obo$current[e$parent], obo$current[e$child], g)

    sig <- match(ids, d$term_id)
    name <- .wrap_label(unname(obo$term_name[ids]), 24)
    label <- name
    if (term_ids) label <- paste0(ids, "\n", label)
    stats_line <- paste0("FDR ", formatC(d$fdr[sig], format = "g", digits = 2))
    if (!is.null(d$extra)) stats_line <- paste0(d$extra[sig], "\n", stats_line)
    label <- ifelse(is.na(sig), label, paste0(label, "\n", stats_line))
    n_lines <- lengths(strsplit(label, "\n", fixed = TRUE))
    box_w <- 2.0
    box_h <- 0.15 * n_lines + 0.12
    node_attrs <- list(label = stats::setNames(rep("", length(ids)), ids),
                       shape = stats::setNames(rep("box", length(ids)), ids),
                       fixedsize = stats::setNames(rep(TRUE, length(ids)), ids),
                       width = stats::setNames(rep(box_w, length(ids)), ids),
                       height = stats::setNames(box_h, ids))
    edge_names <- paste(obo$current[e$parent], obo$current[e$child], sep = "~")
    edge_type <- ifelse(e$weight < 0.7, "part_of", "is_a")
    edge_attrs <- list(arrowhead = stats::setNames(rep("none", nrow(e)), edge_names))
    if (edge_labels && nrow(e)) edge_attrs$label <- stats::setNames(edge_type, edge_names)
    g <- Rgraphviz::layoutGraph(g, nodeAttrs = node_attrs, edgeAttrs = edge_attrs,
                                attrs = list(graph = list(rankdir = "TB", nodesep = "0.25",
                                                          ranksep = if (edge_labels) "0.55" else "0.4"),
                                             edge = list(fontsize = "8")))
    ni <- graph::nodeRenderInfo(g)
    nx <- ni$nodeX[ids]; ny <- ni$nodeY[ids]
    hw <- (ni$lWidth[ids] + ni$rWidth[ids]) / 2
    hh <- ni$height[ids] / 2
    bb <- graph::graphRenderInfo(g, "bbox")
    xr <- c(bb[1, 1], bb[2, 1]); yr <- c(bb[1, 2], bb[2, 2])

    top_margin <- if (is.null(title)) 0.15 else 0.45
    width  <- diff(xr) / 72 + 0.3
    height <- diff(yr) / 72 + 0.15 + top_margin
    scale <- min(1, 45 / max(width, height))
    if (!is.null(file)) .open_png(file, width * scale, height * scale, dpi)
    old <- graphics::par(mai = c(0.15, 0.15, top_margin, 0.15) * scale, ps = 9 * scale, xpd = NA)
    on.exit({
        graphics::par(old)
        if (!is.null(file)) grDevices::dev.off()
    })
    graphics::plot.new()
    graphics::plot.window(xlim = xr, ylim = yr, xaxs = "i", yaxs = "i", asp = 1)

    edge_col <- cols$edge
    splines <- graph::edgeRenderInfo(g)$splines
    for (k in seq_len(nrow(e))) {
        pts <- do.call(rbind, lapply(splines[[edge_names[k]]], Rgraphviz::bezierPoints))
        p_id <- obo$current[e$parent[k]]; c_id <- obo$current[e$child[k]]
        if (pts[1, 2] < pts[nrow(pts), 2]) pts <- pts[nrow(pts):1, , drop = FALSE]
        pts <- rbind(c(pts[1, 1], ny[p_id] - hh[p_id]), pts,
                     c(pts[nrow(pts), 1], ny[c_id] + hh[c_id]))
        graphics::lines(pts[, 1], pts[, 2], col = edge_col[edge_type[k]], lwd = 1.6 * scale,
                        lty = if (edge_type[k] == "part_of") 2 else 1)
        if (edge_labels) {
            mid <- pts[ceiling(nrow(pts) / 2), ]
            tw <- graphics::strwidth(edge_type[k], cex = 0.8) / 2 + 2
            th <- graphics::strheight(edge_type[k], cex = 0.8) / 2 + 1.5
            graphics::rect(mid[1] - tw, mid[2] - th, mid[1] + tw, mid[2] + th, col = "white", border = NA)
            graphics::text(mid[1], mid[2], edge_type[k], cex = 0.8, font = 3, col = edge_col[edge_type[k]])
        }
    }

    fill <- rep("#FFFFFF", length(ids))
    score <- -log10(pmax(d$fdr[sig], 1e-300))
    fill[!is.na(sig)] <- .fill_ramp(score[!is.na(sig)], d$direction[sig[!is.na(sig)]], cols$direction)
    lum <- colSums(grDevices::col2rgb(fill) * c(0.299, 0.587, 0.114)) / 255
    text_col <- ifelse(lum < 0.55, "#FFFFFF", "#000000")
    graphics::rect(nx - hw, ny - hh, nx + hw, ny + hh, col = fill,
                   border = ifelse(is.na(sig), "#888888", "#333333"), lwd = 0.9 * scale)
    graphics::text(nx, ny, label, col = text_col, cex = 1)
    if (!is.null(title)) graphics::mtext(title, side = 3, line = 0.6, cex = 1.1, adj = 0)
    invisible(g)
}

## Themes of the selected result tables as rows keyed by the theme's GO
## term (common ancestor for wang and ic, slim term for slim; "other" for
## slim terms without one). Themes of one table with the same key are
## combined (member terms added, smallest FDR kept).
.theme_rows <- function(x, keys, method) {
    th <- x@info$themes
    if (is.null(th)) stop("No themes in this result; run reduce_terms() first", call. = FALSE)
    out <- lapply(keys, function(k) {
        t <- th[[paste(k, method, sep = "|")]]$themes
        if (is.null(t) || !nrow(t)) return(NULL)
        key <- if ("common_ancestor_id" %in% names(t)) t$common_ancestor_id else t$theme_id
        a <- stats::aggregate(list(n_terms = t$n_terms), list(term_id = key), sum)
        a$fdr <- stats::aggregate(t$min_fdr, list(key), min)$x[match(a$term_id, sort(unique(key)))]
        a$term_name <- t$theme_name[match(a$term_id, key)]
        cbind(family = k, a)
    })
    out <- do.call(rbind, out)
    if (is.null(out)) return(NULL)
    out
}

.theme_method <- function(x, method) {
    available <- unique(sub("^.*\\|", "", names(x@info$themes)))
    if (!length(available)) stop("No themes in this result; run reduce_terms() first", call. = FALSE)
    if (is.null(method)) return(available[1])
    if (!method %in% available)
        stop("No themes from method '", method, "'. Available: ", paste(available, collapse = ", "), call. = FALSE)
    method
}

.theme_dots_plot <- function(x, keys, method, top_themes, max_themes, file, dpi, title = NULL,
                             cols = .plot_colours(x@settings)) {
    d <- .theme_rows(x, keys, method)
    if (is.null(d)) return(invisible(NULL))
    d$n_genes <- d$n_terms
    .dots_draw(d, keys, top_themes, max_themes, file, dpi, title, cols, "Terms")
}

.theme_graph_plot <- function(x, keys, method, top_themes, file, dpi, title = NULL,
                              edge_labels = TRUE, term_ids = TRUE, cols = .plot_colours(x@settings)) {
    if (!requireNamespace("Rgraphviz", quietly = TRUE) || !requireNamespace("graph", quietly = TRUE))
        stop("The GO graph needs the Bioconductor package Rgraphviz: BiocManager::install(\"Rgraphviz\")",
             call. = FALSE)
    d <- .theme_rows(x, keys, method)
    if (is.null(d)) return(invisible(NULL))
    d <- d[d$term_id != "other", , drop = FALSE]
    if (!nrow(d)) return(invisible(NULL))
    d$direction <- .direction_of(d$family, x@info$mode)
    d <- d[order(d$fdr, -d$n_terms, d$term_id), ]
    d <- d[!duplicated(d$term_id), ]
    d <- utils::head(d, top_themes)
    d$extra <- paste0(d$n_terms, if (method == "slim") " terms mapped" else " terms in theme")
    .graph_draw(x, d, file, dpi, title, edge_labels, term_ids, cols)
}

## All figures requested by the option plots, written to FIGURES/ in the
## output directory.
.plot_object <- function(x) {
    s <- x@settings
    dir <- file.path(s$output_dir, "FIGURES")
    keys <- names(x@results)
    parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
    n_files <- 0
    if (x@info$mode == "scores") {
        groups_dots <- split(keys, parts[, 1])
        groups_detail <- split(keys, paste(parts[, 1], parts[, 2], sep = "|"))
    } else {
        groups_dots <- split(keys, paste(parts[, 1], parts[, 3], sep = "|"))
        groups_detail <- as.list(stats::setNames(keys, keys))
    }
    stem <- function(name) .family_stem(name, s)
    for (name in names(groups_dots)) {
        if ("dots" %in% s$plots) {
            f <- file.path(dir, "DOTS", paste0(stem(name), "_dots.png"))
            p <- .dots_plot(x, groups_dots[[name]], s$dots_top_terms, s$dots_max_terms, s$report_fdr,
                            f, s$plot_dpi, title = gsub("|", " ", name, fixed = TRUE))
            n_files <- n_files + !is.null(p)
        }
    }
    for (name in names(groups_detail)) {
        title <- gsub("|", " ", name, fixed = TRUE)
        if ("dendrogram" %in% s$plots) {
            f <- file.path(dir, "DENDROGRAM", paste0(stem(name), "_dendrogram.png"))
            p <- .dendrogram_plot(x, groups_detail[[name]], s$dendrogram_top_terms, s$report_fdr, f, s$plot_dpi, title)
            n_files <- n_files + !is.null(p)
        }
        if ("graph" %in% s$plots) {
            f <- file.path(dir, "GRAPH", paste0(stem(name), "_graph.png"))
            p <- .graph_plot(x, groups_detail[[name]], s$graph_top_terms, s$report_fdr, f, s$plot_dpi, title,
                             s$graph_edge_labels, s$graph_term_ids)
            n_files <- n_files + !is.null(p)
        }
    }
    theme_methods <- if (is.null(x@info$themes)) character(0) else unique(sub("^.*\\|", "", names(x@info$themes)))
    for (m in theme_methods) {
        for (name in names(groups_dots)) {
            if ("theme_dots" %in% s$plots) {
                f <- file.path(dir, "THEME_DOTS", paste0(stem(name), "_", m, "_theme_dots.png"))
                p <- .theme_dots_plot(x, groups_dots[[name]], m, s$dots_top_terms, s$dots_max_terms, f, s$plot_dpi,
                                      title = paste(gsub("|", " ", name, fixed = TRUE), m, "themes"))
                n_files <- n_files + !is.null(p)
            }
        }
        for (name in names(groups_detail)) {
            if ("theme_graph" %in% s$plots) {
                f <- file.path(dir, "THEME_GRAPH", paste0(stem(name), "_", m, "_theme_graph.png"))
                p <- .theme_graph_plot(x, groups_detail[[name]], m, s$graph_top_terms, f, s$plot_dpi,
                                       title = paste(gsub("|", " ", name, fixed = TRUE), m, "themes"),
                                       edge_labels = s$graph_edge_labels, term_ids = s$graph_term_ids)
                n_files <- n_files + !is.null(p)
            }
        }
    }
    .log("Figures (", paste(s$plots, collapse = ", "), "): ", n_files, " files in '", dir, "'")
    invisible(NULL)
}

#' Figures of enrichment results
#'
#' \describe{
#'   \item{plot_dots}{Dot plot of the most significant terms of several
#'     result tables side by side (for example the methods and directions
#'     of one score table, or the modules of a network). Dot size: genes in
#'     the term; colour: -log10 FDR.}
#'   \item{plot_dendrogram}{Significant terms clustered by the genes they
#'     share (average linkage on 1 - shared genes / genes of the larger
#'     term). Label colour: direction (red up, blue down); type face:
#'     significance (bold <= 0.001, plain <= 0.01, italic above).}
#'   \item{plot_go_graph}{The most significant terms with all their
#'     ancestors as a rooted GO graph (layout by Graphviz, package
#'     Rgraphviz). Each box shows the GO ID, the term name and, for
#'     significant terms, the FDR; fill colour shows direction and
#'     significance. Solid grey lines are is_a relations, dashed blue lines
#'     part_of relations; the relation can be printed on each line.}
#' }
#' The same figures are written during a run when the option \code{plots}
#' is set.
#'
#' @param x A \code{GOFERResult}.
#' @param ... Leading parts of result table names, as in
#'   \code{\link{results}}; \code{"*"} matches any value, for example
#'   \code{"AOB", "*", "signed_kme_mwu"} for all modules of a network.
#' @param top_terms Number of terms shown (per result table for
#'   \code{plot_dots}). Default: the options \code{dots_top_terms},
#'   \code{dendrogram_top_terms}, \code{graph_top_terms}.
#' @param fdr Significance cutoff. Default: \code{report_fdr} of the run.
#' @param file PNG file to write. NULL: draw on the current device.
#' @param dpi Resolution of the PNG file.
#' @param colors Named list of colour options for this figure only, for
#'   example \code{list(color_up = "darkorange", dots_colors = "inferno")};
#'   see \code{\link{write_config}} for the colour options. Default: the
#'   colours of the run.
#' @param edge_labels \code{plot_go_graph}: print the relation (is_a,
#'   part_of) on every edge. Default: option \code{graph_edge_labels}.
#' @param term_ids \code{plot_go_graph}: print the GO ID in every box.
#'   Default: option \code{graph_term_ids}.
#' @return Invisibly, the ggplot object (\code{plot_dots}), the clustering
#'   (\code{plot_dendrogram}) or the laid-out graph (\code{plot_go_graph});
#'   NULL when there are no significant terms.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      cpus = 1, use_cache = FALSE, report_fdr = 0.1,
#'                      output_dir = file.path(tempdir(), "example_plots"))
#' plot_dots(res, "example_deseq2")
#' plot_dendrogram(res, "example_deseq2", "mwu")
#' @name gofer-plots
NULL

#' @rdname gofer-plots
#' @export
plot_dots <- function(x, ..., top_terms = NULL, fdr = NULL, file = NULL, dpi = NULL, colors = list()) {
    s <- .update_settings(x, colors)@settings
    .dots_plot(x, .match_keys(x, c(...)), if (is.null(top_terms)) s$dots_top_terms else top_terms,
               s$dots_max_terms, if (is.null(fdr)) s$report_fdr else fdr, file,
               if (is.null(dpi)) s$plot_dpi else dpi, cols = .plot_colours(s))
}

#' @rdname gofer-plots
#' @export
plot_dendrogram <- function(x, ..., top_terms = NULL, fdr = NULL, file = NULL, dpi = NULL, colors = list()) {
    s <- .update_settings(x, colors)@settings
    .dendrogram_plot(x, .match_keys(x, c(...)), if (is.null(top_terms)) s$dendrogram_top_terms else top_terms,
                     if (is.null(fdr)) s$report_fdr else fdr, file, if (is.null(dpi)) s$plot_dpi else dpi,
                     cols = .plot_colours(s))
}

#' @rdname gofer-plots
#' @export
plot_go_graph <- function(x, ..., top_terms = NULL, fdr = NULL, file = NULL, dpi = NULL,
                          edge_labels = NULL, term_ids = NULL, colors = list()) {
    s <- .update_settings(x, colors)@settings
    .graph_plot(x, .match_keys(x, c(...)), if (is.null(top_terms)) s$graph_top_terms else top_terms,
                if (is.null(fdr)) s$report_fdr else fdr, file, if (is.null(dpi)) s$plot_dpi else dpi,
                edge_labels = if (is.null(edge_labels)) s$graph_edge_labels else edge_labels,
                term_ids = if (is.null(term_ids)) s$graph_term_ids else term_ids, cols = .plot_colours(s))
}

#' Figures of themes
#'
#' Figures of the themes from term reduction (see \code{\link{reduce_terms}}).
#' \describe{
#'   \item{plot_theme_dots}{One row per theme, result tables as columns.
#'     Dot size: number of member terms; colour: -log10 of the smallest FDR
#'     of the members. A theme is identified by its GO term (the common
#'     ancestor of its members for wang and ic, the slim term for slim), so
#'     themes with the same term line up across result tables.}
#'   \item{plot_theme_graph}{The GO graph with one box per theme, placed at
#'     the theme's GO term with all its ancestors; the box gives the number of
#'     member terms and the smallest FDR. Slim terms without a slim ancestor
#'     ("other") are not shown.}
#' }
#'
#' @inheritParams gofer-plots
#' @param method Reduction method whose themes are shown (\code{"wang"},
#'   \code{"ic"} or \code{"slim"}). Default: the first method of the run.
#' @param top_themes Number of themes shown (per result table for
#'   \code{plot_theme_dots}). Default: the options \code{dots_top_terms}
#'   and \code{graph_top_terms}.
#' @return Invisibly, the ggplot object or the laid-out graph; NULL when
#'   there are no themes.
#' @examples
#' ext <- system.file("extdata", package = "GOFER")
#' res <- enrich_scores(file.path(ext, "example_deseq2.tsv"),
#'                      file.path(ext, "example_annotations.tsv"),
#'                      go_file = file.path(ext, "example_go.obo"), retain_go_obo = TRUE,
#'                      cpus = 1, use_cache = FALSE, report_fdr = 0.1,
#'                      output_dir = file.path(tempdir(), "example_theme_plots"))
#' plot_theme_dots(res, "example_deseq2")
#' @name gofer-theme-plots
NULL

#' @rdname gofer-theme-plots
#' @export
plot_theme_dots <- function(x, ..., method = NULL, top_themes = NULL, file = NULL, dpi = NULL,
                            colors = list()) {
    s <- .update_settings(x, colors)@settings
    .theme_dots_plot(x, .match_keys(x, c(...)), .theme_method(x, method),
                     if (is.null(top_themes)) s$dots_top_terms else top_themes, s$dots_max_terms,
                     file, if (is.null(dpi)) s$plot_dpi else dpi, cols = .plot_colours(s))
}

#' @rdname gofer-theme-plots
#' @export
plot_theme_graph <- function(x, ..., method = NULL, top_themes = NULL, file = NULL, dpi = NULL,
                             edge_labels = NULL, term_ids = NULL, colors = list()) {
    s <- .update_settings(x, colors)@settings
    .theme_graph_plot(x, .match_keys(x, c(...)), .theme_method(x, method),
                      if (is.null(top_themes)) s$graph_top_terms else top_themes,
                      file, if (is.null(dpi)) s$plot_dpi else dpi,
                      edge_labels = if (is.null(edge_labels)) s$graph_edge_labels else edge_labels,
                      term_ids = if (is.null(term_ids)) s$graph_term_ids else term_ids,
                      cols = .plot_colours(s))
}
