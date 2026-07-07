# Visualization Improvements — Milestone 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a read × position methylation heatmap view to `plot_methylation()` in which reads are drawn one-per-row and ordered by methylation-pattern clustering, to surface epiallele / allele-specific structure.

**Architecture:** A new `view = "heatmap"` mode swaps the packed read panel for a one-read-per-row tile panel that bypasses `pack_reads()`. Reads are pivoted into a read × position matrix (`$sites` is already long-format with shared genomic coordinates), clustered with base `stats::hclust`, and row-ordered by the dendrogram. The heatmap x-axis defaults to genomic coordinates (aligned with the gene/smooth/delta panels); an opt-in compact `index` mode evenly spaces CpG columns. All other panels and the patchwork assembly are unchanged.

**Tech Stack:** R, ggplot2, patchwork, base `stats` (`hclust`, `dist`), testthat (edition 3), devtools, roxygen2.

## Global Constraints

- R package; tests use `testthat` edition 3 via `devtools::test()`.
- All new behaviour argument-gated; default `view = "packed"` output must be **byte-for-byte unchanged** (existing snapshots must still pass).
- No new dependencies — clustering uses base `stats::hclust`/`dist`; rendering uses already-imported `ggplot2`/`patchwork`/`scales`/`ggnewscale`.
- Internal helpers prefixed `.`, reached in tests via `ggmethylation:::`.
- Exported functions use roxygen2 (`markdown = TRUE`); run `devtools::document()` after doc edits.
- **Assumes Milestone 1 is merged** — reuses `theme_ggmethylation()`, `.PROB_GRADIENT`, `.classify_calls()`, and the `call_threshold` argument. Independent of Milestone 2.
- Out of scope for v1 (document as such): dendrogram strip, variant/CIGAR/SA overlays on the heatmap panel, and heatmap for `multi_methylation_data` (falls back to packed with a message).

---

### Task 1: Read × position matrix builder

**Files:**
- Create: `R/read_matrix.R`
- Create: `tests/testthat/test-read_matrix.R`

**Interfaces:**
- Consumes: `data$sites` (`read_name`, `position`, `mod_prob`).
- Produces: `.build_read_matrix(sites)` → numeric matrix, rownames = read_name, colnames = **numerically sorted** genomic positions, `NA` where a read has no site at a position; `0×0` matrix on empty input.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-read_matrix.R
test_that("build_read_matrix collapses shared positions across reads", {
  sites <- data.frame(
    read_name = c("r1", "r1", "r2", "r2"),
    position  = c(100L, 200L, 100L, 300L),
    mod_prob  = c(0.9, 0.1, 0.2, 0.8),
    stringsAsFactors = FALSE
  )
  m <- ggmethylation:::.build_read_matrix(sites)
  expect_equal(sort(rownames(m)), c("r1", "r2"))
  expect_equal(colnames(m), c("100", "200", "300"))   # numeric order
  expect_equal(m["r1", "100"], 0.9)
  expect_true(is.na(m["r2", "200"]))                   # r2 has no site at 200
})

test_that("build_read_matrix orders columns numerically not lexically", {
  sites <- data.frame(
    read_name = c("r1", "r1"),
    position  = c(900L, 1000L),
    mod_prob  = c(0.5, 0.6),
    stringsAsFactors = FALSE
  )
  m <- ggmethylation:::.build_read_matrix(sites)
  expect_equal(colnames(m), c("900", "1000"))
})

test_that("build_read_matrix returns 0x0 on empty input", {
  m <- ggmethylation:::.build_read_matrix(
    data.frame(read_name = character(0), position = integer(0),
               mod_prob = numeric(0), stringsAsFactors = FALSE)
  )
  expect_equal(dim(m), c(0L, 0L))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-read_matrix.R")'`
Expected: FAIL — `.build_read_matrix` not found.

- [ ] **Step 3: Implement**

```r
# R/read_matrix.R

# Pivot long-format sites into a reads x positions matrix.
# tapply over two factors yields character-sorted column names, so columns are
# reordered by numeric genomic position. NA marks positions a read does not cover.
.build_read_matrix <- function(sites) {
  if (is.null(sites) || nrow(sites) == 0L) {
    return(matrix(numeric(0L), nrow = 0L, ncol = 0L))
  }
  m <- tapply(sites$mod_prob,
              list(sites$read_name, sites$position),
              FUN = mean)
  m <- matrix(as.double(m), nrow = nrow(m), ncol = ncol(m),
              dimnames = dimnames(m))
  ord <- order(as.numeric(colnames(m)))
  m[, ord, drop = FALSE]
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-read_matrix.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/read_matrix.R tests/testthat/test-read_matrix.R
git commit -m "feat: add read x position matrix builder"
```

---

### Task 2: Methylation-pattern read ordering

**Files:**
- Create: `R/cluster_reads.R`
- Create: `tests/testthat/test-cluster_reads.R`

**Interfaces:**
- Consumes: matrix from `.build_read_matrix()` (Task 1); `.classify_calls()` is not needed here — binarization is inline against `threshold`.
- Produces: `.cluster_read_order(mat, metric = "euclidean", method = "ward.D2", threshold = 0.5)` → character vector: a permutation of all `rownames(mat)`. Reads sharing no positions with any other (all-`NaN` distances) are appended at the end in original order; `≤ 2` clusterable reads returns input order.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-cluster_reads.R
test_that("cluster_read_order returns a permutation of all reads", {
  set.seed(1)
  mat <- matrix(runif(20), nrow = 4,
                dimnames = list(paste0("r", 1:4), paste0("p", 1:5)))
  ord <- ggmethylation:::.cluster_read_order(mat)
  expect_setequal(ord, rownames(mat))
})

test_that("cluster_read_order groups identical patterns adjacently", {
  mat <- rbind(
    a1 = c(0.9, 0.9, 0.1, 0.1),
    b1 = c(0.1, 0.1, 0.9, 0.9),
    a2 = c(0.9, 0.9, 0.1, 0.1),
    b2 = c(0.1, 0.1, 0.9, 0.9)
  )
  colnames(mat) <- paste0("p", 1:4)
  ord <- ggmethylation:::.cluster_read_order(mat)
  # a1/a2 adjacent and b1/b2 adjacent
  pos <- match(c("a1", "a2", "b1", "b2"), ord)
  expect_equal(abs(pos[1] - pos[2]), 1)
  expect_equal(abs(pos[3] - pos[4]), 1)
})

test_that("cluster_read_order appends non-overlapping reads at the end", {
  mat <- rbind(
    r1 = c(0.9, 0.1, NA,  NA),
    r2 = c(0.8, 0.2, NA,  NA),
    r3 = c(NA,  NA,  0.5, 0.5)   # shares no positions with r1/r2
  )
  colnames(mat) <- paste0("p", 1:4)
  ord <- ggmethylation:::.cluster_read_order(mat)
  expect_equal(ord[length(ord)], "r3")
})

test_that("cluster_read_order binary metric returns a permutation", {
  mat <- rbind(r1 = c(0.9, 0.1), r2 = c(0.1, 0.9), r3 = c(0.9, 0.9))
  colnames(mat) <- c("p1", "p2")
  ord <- ggmethylation:::.cluster_read_order(mat, metric = "binary", threshold = 0.5)
  expect_setequal(ord, rownames(mat))
})

test_that("cluster_read_order returns input order for <= 2 reads", {
  mat <- rbind(r1 = c(0.1, 0.2), r2 = c(0.3, 0.4))
  colnames(mat) <- c("p1", "p2")
  expect_equal(ggmethylation:::.cluster_read_order(mat), c("r1", "r2"))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-cluster_reads.R")'`
Expected: FAIL — `.cluster_read_order` not found.

- [ ] **Step 3: Implement**

Base `dist()` excludes `NA` pairwise; reads sharing no positions yield `NaN`.
Drop all-`NaN` reads to the end, impute remaining `NaN` to the max finite
distance, then cluster.

```r
# R/cluster_reads.R

# Order reads by methylation-pattern similarity via hierarchical clustering.
# Robust to NA (uncovered positions) and to reads that share no positions.
.cluster_read_order <- function(mat, metric = "euclidean",
                                method = "ward.D2", threshold = 0.5) {
  rn <- rownames(mat)
  if (is.null(rn) || length(rn) <= 2L) return(rn)

  x <- mat
  if (identical(metric, "binary")) {
    x[] <- ifelse(is.na(x), NA_real_, as.double(x >= threshold))
    d <- stats::dist(x, method = "binary")
  } else {
    d <- stats::dist(x, method = "euclidean")
  }

  dm  <- as.matrix(d)
  bad <- rn[apply(dm, 1L, function(r) all(is.na(r) | is.nan(r)))]
  ok  <- setdiff(rn, bad)
  if (length(ok) <= 2L) return(c(ok, bad))

  d_ok <- stats::as.dist(dm[ok, ok, drop = FALSE])
  finite_max <- max(d_ok[is.finite(d_ok)], na.rm = TRUE)
  d_ok[is.na(d_ok) | is.nan(d_ok)] <- finite_max

  hc <- stats::hclust(d_ok, method = method)
  c(ok[hc$order], bad)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-cluster_reads.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/cluster_reads.R tests/testthat/test-cluster_reads.R
git commit -m "feat: add methylation-pattern read clustering order"
```

---

### Task 3: Heatmap panel builder (genomic mode)

**Files:**
- Create: `R/build_heatmap_panel.R`
- Create: `tests/testthat/test-build_heatmap_panel.R`

**Interfaces:**
- Consumes: a `methylation_data` whose `$reads` already has a one-read-per-lane
  `lane` column; `theme_ggmethylation()`, `.PROB_GRADIENT`, `.classify_calls()`
  (M1).
- Produces: `build_heatmap_panel(data, region_start, region_end, colour_low, colour_high, group_colours, call_mode = "continuous", call_threshold = 0.5, heatmap_x = "genomic", separator_lanes = numeric(0), show_x_axis = FALSE)` → a ggplot.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-build_heatmap_panel.R
# Reuse make_test_data pattern (copy the minimal constructor here).
make_test_data <- function(reads_df, sites_df, region_start = 1000L, region_end = 2000L) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = region_start, end = region_end))
  structure(list(reads = reads_df, sites = sites_df, region = gr,
                 mod_code = "m", group_tag = NULL,
                 cigar_features = data.frame()),
            class = "methylation_data")
}

test_that("build_heatmap_panel returns a ggplot with tile + background layers", {
  reads <- data.frame(
    read_name = c("r1", "r2"), start = c(1000L, 1100L), end = c(1900L, 2000L),
    strand = "+", lane = c(1L, 2L), stringsAsFactors = FALSE
  )
  sites <- data.frame(
    read_name = rep(c("r1", "r2"), each = 3),
    position = rep(c(1200L, 1400L, 1600L), 2),
    mod_prob = c(0.9, 0.1, 0.8, 0.2, 0.7, 0.3),
    mod_code = "m", stringsAsFactors = FALSE
  )
  md <- make_test_data(reads, sites)
  p <- ggmethylation:::build_heatmap_panel(
    md, 1000L, 2000L, colour_low = "#BDBDBD", colour_high = "#C62828",
    group_colours = NULL
  )
  expect_s3_class(p, "ggplot")
  # at least a background segment layer and a tile layer
  expect_gte(length(p$layers), 2L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-build_heatmap_panel.R")'`
Expected: FAIL — `build_heatmap_panel` not found.

- [ ] **Step 3: Implement (genomic mode; index handled in Task 4)**

```r
# R/build_heatmap_panel.R

#' Build a one-read-per-row methylation heatmap panel
#'
#' @keywords internal
build_heatmap_panel <- function(data, region_start, region_end,
                                colour_low, colour_high, group_colours,
                                call_mode = "continuous", call_threshold = 0.5,
                                heatmap_x = "genomic",
                                separator_lanes = numeric(0),
                                show_x_axis = FALSE) {
  reads <- data$reads
  sites <- merge(
    data$sites,
    reads[, c("read_name", "lane"), drop = FALSE],
    by = "read_name"
  )

  # x mapping: genomic uses true position; index remaps to column rank (Task 4)
  sites$.x <- sites$position
  tile_w <- (region_end - region_start) * 0.006

  half_height <- 0.45
  p <- ggplot2::ggplot()

  # Faint full-extent row background (group fill when grouped, else grey)
  if (!is.null(data$group_tag) && "group" %in% names(reads)) {
    p <- p +
      ggplot2::geom_rect(
        data = reads,
        ggplot2::aes(xmin = .data$start, xmax = .data$end,
                     ymin = .data$lane - half_height,
                     ymax = .data$lane + half_height, fill = .data$group),
        alpha = 0.25, colour = NA
      )
    if (!is.null(group_colours)) {
      p <- p + ggplot2::scale_fill_manual(values = group_colours,
                                          na.value = "grey50", name = "Group")
    }
    p <- p + ggnewscale::new_scale_fill()
  } else {
    p <- p +
      ggplot2::geom_rect(
        data = reads,
        ggplot2::aes(xmin = .data$start, xmax = .data$end,
                     ymin = .data$lane - half_height,
                     ymax = .data$lane + half_height),
        fill = "grey90", colour = NA
      )
  }

  # CpG tiles
  if (identical(call_mode, "binary") && nrow(sites) > 0L) {
    sites$.call <- .classify_calls(sites$mod_prob, call_threshold, NULL)
    p <- p +
      ggplot2::geom_tile(
        data = sites,
        ggplot2::aes(x = .data$.x, y = .data$lane, fill = .data$.call),
        width = tile_w, height = 2 * half_height
      ) +
      ggplot2::scale_fill_manual(
        values = c(unmethylated = colour_low, methylated = colour_high),
        name = "Call"
      )
  } else if (nrow(sites) > 0L) {
    p <- p +
      ggplot2::geom_tile(
        data = sites,
        ggplot2::aes(x = .data$.x, y = .data$lane, fill = .data$mod_prob),
        width = tile_w, height = 2 * half_height
      ) +
      ggplot2::scale_fill_gradient(
        low = colour_low, high = colour_high, limits = c(0, 1),
        name = "Modification\nprobability"
      )
  }

  if (length(separator_lanes) > 0L) {
    p <- p + ggplot2::geom_hline(
      yintercept = separator_lanes, linetype = "dashed",
      colour = "grey40", linewidth = 0.4
    )
  }

  p <- p +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(xlim = c(region_start, region_end)) +
    theme_ggmethylation() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL)

  if (!show_x_axis) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }
  p
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::document(); devtools::test_file("tests/testthat/test-build_heatmap_panel.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/build_heatmap_panel.R tests/testthat/test-build_heatmap_panel.R man/
git commit -m "feat: add genomic-coordinate heatmap panel builder"
```

---

### Task 4: Compact index x-axis mode

**Files:**
- Modify: `R/build_heatmap_panel.R`
- Modify: `tests/testthat/test-build_heatmap_panel.R`

**Interfaces:**
- Produces: when `heatmap_x = "index"`, each site's x becomes its column rank
  (1..K over sorted unique positions), tiles drawn contiguous (`width = 1`).

- [ ] **Step 1: Write the failing test**

```r
test_that("index mode remaps x to monotonic column ranks", {
  reads <- data.frame(read_name = "r1", start = 1000L, end = 2000L,
                      strand = "+", lane = 1L, stringsAsFactors = FALSE)
  sites <- data.frame(read_name = "r1",
                      position = c(1600L, 1200L, 1400L),
                      mod_prob = c(0.7, 0.9, 0.1), mod_code = "m",
                      stringsAsFactors = FALSE)
  md <- make_test_data(reads, sites)
  # Access the internal remap helper directly
  idx <- ggmethylation:::.position_to_index(sites$position)
  # 1200 -> 1, 1400 -> 2, 1600 -> 3 regardless of input order
  expect_equal(idx[sites$position == 1200L], 1L)
  expect_equal(idx[sites$position == 1600L], 3L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-build_heatmap_panel.R")'`
Expected: FAIL — `.position_to_index` not found.

- [ ] **Step 3: Implement the remap and wire it in**

Add the helper and branch the x mapping / coord in `build_heatmap_panel()`:

```r
# Map genomic positions to 1-based column ranks over sorted unique positions.
.position_to_index <- function(positions) {
  levels_sorted <- sort(unique(positions))
  match(positions, levels_sorted)
}
```

In `build_heatmap_panel()`, replace the x-mapping / width / coord blocks with:

```r
  if (identical(heatmap_x, "index")) {
    sites$.x <- .position_to_index(sites$position)
    tile_w   <- 1
    x_limits <- c(0.5, max(sites$.x, 1) + 0.5)
    # background row spans full index width (start/end lose meaning in index space)
    reads$.xmin <- 0.5
    reads$.xmax <- max(sites$.x, 1) + 0.5
  } else {
    sites$.x <- sites$position
    tile_w   <- (region_end - region_start) * 0.006
    x_limits <- c(region_start, region_end)
    reads$.xmin <- reads$start
    reads$.xmax <- reads$end
  }
```

Update the two `geom_rect` background layers to use `.xmin`/`.xmax`, the
`geom_tile` `width = tile_w`, and the final `coord_cartesian(xlim = x_limits)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-build_heatmap_panel.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/build_heatmap_panel.R tests/testthat/test-build_heatmap_panel.R
git commit -m "feat: add compact index x-axis mode to heatmap panel"
```

---

### Task 5: Wire `view = "heatmap"` into `plot_methylation()`

**Files:**
- Modify: `R/plot_methylation.R` (signature; row-order + lane assignment; panel swap; multi-sample guard)
- Create: `tests/testthat/test-plot_methylation_heatmap.R`

**Interfaces:**
- Consumes: `.build_read_matrix()`, `.cluster_read_order()`, `build_heatmap_panel()`, `.ordered_plot_groups()`, `.match_plot_group()`.
- Produces: `plot_methylation(view = c("packed","heatmap"), cluster = TRUE, cluster_metric = c("euclidean","binary"), cluster_method = "ward.D2", heatmap_x = c("genomic","index"))`.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-plot_methylation_heatmap.R
# (reuse make_test_data constructor as in test-build_heatmap_panel.R)

test_that("plot_methylation heatmap view renders (ungrouped, clustered)", {
  reads <- data.frame(
    read_name = paste0("r", 1:4), start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- do.call(rbind, lapply(1:4, function(i) data.frame(
    read_name = paste0("r", i),
    position = seq(1100, 1900, length.out = 6),
    mod_prob = if (i <= 2) c(.9,.9,.1,.1,.9,.9) else c(.1,.1,.9,.9,.1,.1),
    mod_code = "m", stringsAsFactors = FALSE
  )))
  md <- make_test_data(reads, sites)
  p <- ggmethylation::plot_methylation(md, view = "heatmap",
                                       show_supplementary = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("plot_methylation packed view is the default", {
  reads <- data.frame(
    read_name = "r1", start = 1000L, end = 2000L, strand = "+",
    lane = 0L, mean_mod_prob = 0.5, clip_side = NA_character_,
    sa_chrom = NA_character_, stringsAsFactors = FALSE
  )
  sites <- data.frame(read_name = "r1", position = seq(1100, 1900, length.out = 5),
                      mod_prob = c(.1,.3,.5,.7,.9), mod_code = "m",
                      stringsAsFactors = FALSE)
  md <- make_test_data(reads, sites)
  expect_no_error(ggmethylation::plot_methylation(md, show_supplementary = FALSE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-plot_methylation_heatmap.R")'`
Expected: FAIL — unused argument `view`.

- [ ] **Step 3: Add args + `match.arg`**

Add to the `plot_methylation()` signature:
`view = c("packed", "heatmap")`, `cluster = TRUE`,
`cluster_metric = c("euclidean", "binary")`, `cluster_method = "ward.D2"`,
`heatmap_x = c("genomic", "index")`. At the top of the body:

```r
  view          <- match.arg(view)
  cluster_metric <- match.arg(cluster_metric)
  heatmap_x     <- match.arg(heatmap_x)
```

- [ ] **Step 4: Row-order + lane assignment for heatmap view**

Replace the lane-packing block (section 5) so that in heatmap view each read
gets its own lane in cluster (or `sort_by`) order, per group:

```r
  data$reads$lane <- integer(nrow(data$reads))
  separator_lanes <- numeric(0)

  assign_rows <- function(idx) {
    # idx: row indices into data$reads for one group (or all). Returns lanes.
    rn <- data$reads$read_name[idx]
    if (view == "heatmap" && isTRUE(cluster)) {
      sub_sites <- data$sites[data$sites$read_name %in% rn, , drop = FALSE]
      mat   <- .build_read_matrix(sub_sites)
      order_rn <- .cluster_read_order(mat, cluster_metric, cluster_method,
                                      call_threshold)
      # reads with no sites are absent from mat: append them
      order_rn <- c(order_rn, setdiff(rn, order_rn))
    } else {
      order_rn <- rn  # already sorted by sort_by earlier
    }
    stats::setNames(match(rn, order_rn), rn)
  }

  if (!is.null(data$group_tag)) {
    groups_ordered <- .ordered_plot_groups(data$reads$group)
    lane_offset <- 0L
    for (grp in groups_ordered) {
      idx <- which(.match_plot_group(data$reads$group, grp))
      if (view == "heatmap") {
        data$reads$lane[idx] <- assign_rows(idx) + lane_offset
      } else {
        data$reads$lane[idx] <- pack_reads(data$reads[idx, ],
                                            clip_side = data$reads$clip_side[idx]) + lane_offset
      }
      lane_offset <- max(data$reads$lane[idx]) + 2L
      separator_lanes <- c(separator_lanes, lane_offset - 1L)
    }
    separator_lanes <- separator_lanes[-length(separator_lanes)]
  } else {
    if (view == "heatmap") {
      data$reads$lane <- assign_rows(seq_len(nrow(data$reads)))
    } else {
      data$reads$lane <- pack_reads(data$reads, clip_side = data$reads$clip_side)
    }
  }
```

- [ ] **Step 5: Swap the top panel + guards**

Replace the `build_read_panel(...)` call for `p_top` with a branch:

```r
  if (view == "heatmap") {
    if (identical(heatmap_x, "index")) {
      message("heatmap_x = 'index': gene/smooth/delta panels use genomic ",
              "coordinates and will not align horizontally with the heatmap.")
    }
    p_top <- build_heatmap_panel(
      data, region_start, region_end,
      colour_low = colour_low, colour_high = colour_high,
      group_colours = group_colours,
      call_mode = call_mode, call_threshold = call_threshold,
      heatmap_x = heatmap_x, separator_lanes = separator_lanes,
      show_x_axis = FALSE
    )
  } else {
    p_top <- build_read_panel(...)  # unchanged existing call
  }
```

In `.plot_multi_methylation()`, guard heatmap at the top:

```r
  if (identical(view, "heatmap")) {
    message("view = 'heatmap' is not supported for multi-sample objects; ",
            "using packed view.")
    view <- "packed"
  }
```

(Thread `view`/`cluster`/`cluster_metric`/`cluster_method`/`heatmap_x` into the
`.plot_multi_methylation()` signature so the arg exists to guard.)

- [ ] **Step 6: Run tests**

Run: `Rscript -e 'devtools::test_file("tests/testthat/test-plot_methylation_heatmap.R")'`
Expected: PASS. Then the full suite to prove packed view is unchanged:

Run: `Rscript -e 'devtools::test()'`
Expected: PASS (all existing snapshots/tests still green).

- [ ] **Step 7: Commit**

```bash
git add R/plot_methylation.R tests/testthat/test-plot_methylation_heatmap.R
git commit -m "feat: add heatmap view with methylation-pattern clustering (view='heatmap')"
```

---

### Task 6: Documentation, NEWS, vignette note

**Files:**
- Modify: `R/plot_methylation.R` (roxygen)
- Modify: `NEWS.md`
- Modify: `vignettes/ggmethylation.Rmd`

- [ ] **Step 1: Roxygen for new args**

Document `@param view`, `@param cluster`, `@param cluster_metric`,
`@param cluster_method`, `@param heatmap_x` on `plot_methylation()`. Note that
`cluster*`/`heatmap_x` apply only to `view = "heatmap"`, that grouped data
clusters within each group, that `heatmap_x = "index"` breaks x-alignment with
genomic panels, and that heatmap view is unsupported for multi-sample objects.

Run: `Rscript -e 'devtools::document()'`

- [ ] **Step 2: NEWS bullet**

```markdown
* New `view = "heatmap"` renders reads one-per-row ordered by
  methylation-pattern clustering (`cluster`, `cluster_metric`,
  `cluster_method`), with genomic (default) or compact `index` x-axis.
```

- [ ] **Step 3: Vignette paragraph**

Add a short "Heatmap view" section to `vignettes/ggmethylation.Rmd` showing
`plot_methylation(md, view = "heatmap")` and `... cluster_metric = "binary"`.

- [ ] **Step 4: Check + commit**

Run: `Rscript -e 'devtools::check()'`
Expected: no new NOTEs/WARNINGs (no new dependencies).

```bash
git add R/plot_methylation.R NEWS.md vignettes/ggmethylation.Rmd man/
git commit -m "docs: document heatmap view and clustering options"
```

---

## Self-Review

- **Spec coverage:** matrix → Task 1; clustering (euclidean default + binary) → Task 2; genomic heatmap → Task 3; index mode → Task 4; view wiring + within-group clustering + multi-sample guard → Task 5; docs → Task 6.
- **Type consistency:** `.build_read_matrix` (rownames = read_name) → `.cluster_read_order` (returns read_name permutation) → `match(rn, order_rn)` lane assignment in Task 5 → `build_heatmap_panel` consumes `$reads$lane`. `.position_to_index` produced/consumed within Tasks 3–4.
- **Risk flags:** numeric column ordering (Task 1); NaN/no-overlap distances (Task 2); index-mode panel misalignment (message + docs, Tasks 4–5); reads with zero sites appended after clustering (Task 5 `assign_rows`); multi-sample fallback (Task 5).
- **Ordering:** Tasks 1→2→(3→4)→5→6. Assumes Milestone 1 merged (`theme_ggmethylation`, `.PROB_GRADIENT`, `.classify_calls`, `call_threshold`); independent of Milestone 2.
