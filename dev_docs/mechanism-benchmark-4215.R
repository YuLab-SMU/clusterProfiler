#!/usr/bin/env Rscript
# Offline exploratory benchmark; run from the package root. No package load/network.
args <- commandArgs(trailingOnly = TRUE)
self_test <- "--self-test" %in% args
freeze <- "--freeze" %in% args
args <- args[!args %in% c("--self-test", "--freeze")]
stopifnot(length(args) <= 1L)
root <- if (length(args)) args[[1L]] else "../adjudicate-interpret/results/audit-benchmark-4215"
input <- file.path(root, "inputs")
outdir <- file.path(root, "outputs")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                         quote = FALSE, na = "NA", fileEncoding = "UTF-8")
manifest <- function(paths, base) data.frame(
    file = substring(paths, nchar(base) + 2L), bytes = file.info(paths)$size,
    md5 = unname(tools::md5sum(paths)), stringsAsFactors = FALSE)
if (freeze) {
    if (dir.exists(input)) stop("Frozen inputs already exist; refusing to overwrite.")
    original <- dirname(root)
    rel <- c("all11_de.rds", "p3/p3_new_de.rds", "collectri_net.tsv",
             "p6_discover_verify.tsv", "p7_module_all17.tsv")
    paths <- file.path(original, rel)
    stopifnot(all(file.exists(paths)))
    profiles <- c(readRDS(paths[1L]), readRDS(paths[2L]))
    profiles <- lapply(profiles, function(z) z[c("label", "driver", "signed")])
    stopifnot(length(profiles) == 17L)
    dir.create(input, recursive = TRUE)
    saveRDS(profiles, file.path(input, "profiles.rds"), compress = "xz", version = 2)
    stopifnot(all(file.copy(paths[3:5], input)))
    write_tsv(manifest(paths, original), file.path(input, "original-manifest.tsv"))
    # Study assignments are manual, evidenced by generation scripts, not inferred from row IDs.
    studies <- c("GSE319526", "unknown", "GSE246363", "GSE286890", "GSE311462",
                 "GSE275316", "GSE286515", "GSE285498", "GSE348245", "GSE315572",
                 "GSE276150", rep("GSE245918", 3), rep("GSE348821", 2), "GSE351288")
    evidence <- c("36_all11_de.R:92-99", "36_all11_de.R:101-108; accession not recorded",
                  "36_all11_de.R:110-119", "36_all11_de.R:121-129", "36_all11_de.R:131-141",
                  "36_all11_de.R:143-152", "36_all11_de.R:154-165", "36_all11_de.R:167-176",
                  "36_all11_de.R:178-188", "36_all11_de.R:190-202", "36_all11_de.R:204-213",
                  rep("50_p3_new_datasets.R:74-91", 3), rep("50_p3_new_datasets.R:93-107", 2),
                  "50_p3_new_datasets.R:109-122")
    label <- vapply(profiles, `[[`, character(1), "label")
    driver <- vapply(profiles, function(z) z$driver[1L], character(1))
    driver[driver == "GR"] <- "NR3C1"
    meta <- data.frame(dataset = label, expected_label = driver,
        declared_drivers = vapply(profiles, function(z) paste(z$driver, collapse = ","), character(1)),
        design = c(rep("activation", 2), "KO", rep("activation", 14)),
        study_group = studies, study_evidence = evidence,
        label_evidence = c(rep("36_all11_de.R:finish calls; first driver only; GR=NR3C1", 11),
                           rep("50_p3_new_datasets.R:spec drv", 6)))
    write_tsv(meta, file.path(input, "metadata.tsv"))
    p <- sort(list.files(input, full.names = TRUE))
    write_tsv(manifest(p, input), file.path(root, "input-manifest.tsv"))
}
expected <- read.delim(file.path(root, "input-manifest.tsv"), stringsAsFactors = FALSE)
actual <- manifest(file.path(input, expected$file), input)
stopifnot(isTRUE(all.equal(actual, expected, check.attributes = TRUE)))
profiles <- readRDS(file.path(input, "profiles.rds"))
meta <- read.delim(file.path(input, "metadata.tsv"), stringsAsFactors = FALSE)
net <- read.delim(file.path(input, "collectri_net.tsv"), stringsAsFactors = FALSE)
hist <- read.delim(file.path(input, "p6_discover_verify.tsv"), stringsAsFactors = FALSE)
modules <- read.delim(file.path(input, "p7_module_all17.tsv"), stringsAsFactors = FALSE)
source("R/mechanism_discover.R", local = TRUE)
stopifnot(length(profiles) == 17L, nrow(meta) == 17L, !anyDuplicated(meta$dataset),
          identical(meta$dataset, vapply(profiles, `[[`, character(1), "label")),
          setequal(meta$dataset, unique(hist$dataset)), setequal(meta$dataset, modules$dataset))

# Full candidate universe: never infer absent/unscorable from a truncated top20 table.
rows <- list(); rankings <- list()
for (i in seq_along(profiles)) {
    z <- profiles[[i]]; truth <- meta$expected_label[i]
    for (mode in c("score", "abs")) {
        d <- suppressMessages(mechanism_discover(z$signed, net,
            top = length(unique(net$source)), rank_by = mode, min_targets = 5L,
            keep_self_targets = FALSE))
        ntested <- attr(d, "n_factors_tested")
        stopifnot(nrow(d) == ntested, !anyDuplicated(d$source), all(is.finite(d$score)))
        r <- match(truth, d$source)
        status <- if (!truth %in% net$source) "missing_from_network" else
                  if (is.na(r)) "unscorable" else "scored"
        rows[[length(rows) + 1L]] <- data.frame(dataset = z$label, expected_label = truth,
            design = meta$design[i], study_group = meta$study_group[i], rank_by = mode,
            status = status, truncated = FALSE, expected_rank = r, n_tested = ntested,
            score = if (is.na(r)) NA_real_ else d$score[r],
            hit1 = !is.na(r) && r <= 1L, hit5 = !is.na(r) && r <= 5L,
            hit10 = !is.na(r) && r <= 10L)
        rankings[[length(rankings) + 1L]] <- cbind(dataset = z$label, rank_by = mode,
                                                as.data.frame(d))
    }
}
ranks <- do.call(rbind, rows)
full <- do.call(rbind, rankings)
metric <- function(d, group_type, group) data.frame(group_type = group_type, group = group,
    rank_by = d$rank_by[1L], n = nrow(d), scored = sum(d$status == "scored"),
    hit1 = sum(d$hit1), hit5 = sum(d$hit5), hit10 = sum(d$hit10),
    rate1 = mean(d$hit1), rate5 = mean(d$hit5), rate10 = mean(d$hit10))
summary <- list()
for (mode in c("score", "abs")) {
    d <- ranks[ranks$rank_by == mode, ]
    summary[[length(summary) + 1L]] <- metric(d, "all", "all17")
    for (kind in c("design", "study_group", "expected_label")) {
        for (g in unique(d[[kind]])) summary[[length(summary) + 1L]] <-
            metric(d[d[[kind]] == g, ], kind, g)
    }
    # Prespecified design rule, not selected using observed ranks.
}
selected <- ranks[(ranks$design == "KO" & ranks$rank_by == "abs") |
                  (ranks$design == "activation" & ranks$rank_by == "score"), ]
selected$rank_by <- "design_prespecified"
summary[[length(summary) + 1L]] <- metric(selected, "all", "all17")
for (kind in c("design", "study_group", "expected_label")) {
    for (g in unique(selected[[kind]])) summary[[length(summary) + 1L]] <-
        metric(selected[selected[[kind]] == g, ], kind, g)
}
summary <- do.call(rbind, summary)
# Equal weight per known study/driver, not per repeated contrast; unknown not a unique study.
macro <- do.call(rbind, lapply(c("score", "abs", "design_prespecified"), function(mode) {
    do.call(rbind, lapply(c("study_group", "expected_label"), function(kind) {
        s <- summary[summary$rank_by == mode & summary$group_type == kind, ]
        if (kind == "study_group") s <- s[s$group != "unknown", ]
        data.frame(rank_by = mode, grouping = kind, groups = nrow(s),
                   rate1 = mean(s$rate1), rate5 = mean(s$rate5), rate10 = mean(s$rate10),
                   unknown_excluded = if (kind == "study_group") 1L else 0L)
    }))
}))
# Frozen historical modules: direct top1 neighbours within historical top5 only.
# Exact expectation: keep top1, add a uniform (k-1)-subset of the SAME remaining four.
modules <- modules[match(meta$dataset, modules$dataset), ]
modrows <- lapply(seq_len(nrow(meta)), function(i) {
    h <- hist[hist$dataset == meta$dataset[i], ]; h <- h[order(h$rank), ]
    m <- strsplit(modules$module[i], ",", fixed = TRUE)[[1L]]
    current <- full[full$dataset == meta$dataset[i] & full$rank_by == "score", ]
    stopifnot(identical(head(current$source, 5L), h$candidate))
    targets <- function(tf) unique(net$target[net$source == tf])
    top_targets <- targets(h$candidate[1L])
    jaccard <- vapply(h$candidate, function(tf) {
        other <- targets(tf)
        length(intersect(top_targets, other)) / length(union(top_targets, other))
    }, numeric(1))
    # Historical threshold is fixed at 0.10; this audit does not select a threshold.
    rebuilt <- h$candidate[seq_len(nrow(h)) == 1L | (!is.na(jaccard) & jaccard >= 0.10)]
    stopifnot(identical(m, rebuilt))
    stopifnot(nrow(h) == 5L, !anyDuplicated(h$candidate), all(m %in% h$candidate),
              length(m) == modules$module_size[i], h$candidate[1L] == modules$top1[i])
    truth <- meta$expected_label[i]
    data.frame(dataset = meta$dataset[i], expected_label = truth,
        candidate_set = paste(h$candidate, collapse = ","), candidate_set_size = nrow(h),
        module = paste(m, collapse = ","), module_size = length(m),
        top1_match = h$candidate[1L] == truth, coverage = truth %in% m,
        label_in_candidate_set = truth %in% h$candidate,
        matched_size_random_expected = if (truth == h$candidate[1L]) 1 else
            as.numeric(truth %in% h$candidate[-1L]) * (length(m) - 1L)/(nrow(h) - 1L))
})
mod <- do.call(rbind, modrows)
module_summary <- data.frame(n = 17L, coverage = sum(mod$coverage),
    coverage_rate = mean(mod$coverage), top1_matches = sum(mod$top1_match),
    mean_size = mean(mod$module_size), random_expected_matches = sum(mod$matched_size_random_expected),
    random_expected_rate = mean(mod$matched_size_random_expected),
    label_in_top5 = sum(mod$label_in_candidate_set))
# Historical A/B calls are not rerun and cannot represent current fixed grading.
hist$expected_label_match <- hist$candidate == meta$expected_label[match(hist$dataset, meta$dataset)]
hist$historical_pass <- hist$verdict %in% c("A", "B")
historical_summary <- data.frame(candidate_rows = nrow(hist), historical_pass = sum(hist$historical_pass),
    expected_label_matches_among_pass = sum(hist$historical_pass & hist$expected_label_match))
# Regression guards lock denominators and the original unmodified signed-score result.
b <- ranks[ranks$rank_by == "score", ]
stopifnot(nrow(b) == 17L, all(b$status == "scored"),
          identical(as.integer(c(sum(b$hit1), sum(b$hit5), sum(b$hit10))), c(6L, 10L, 13L)),
          sum(mod$coverage) == 7L, sum(mod$top1_match) == 6L,
          abs(sum(mod$matched_size_random_expected) - 6.5) < 1e-12,
          historical_summary$historical_pass == 37L,
          historical_summary$expected_label_matches_among_pass == 7L)
if (self_test) {
    stopifnot(!any(ranks$truncated), all(mod$matched_size_random_expected >= 0),
        all(mod$matched_size_random_expected <= 1), sum(meta$study_group == "GSE245918") == 3L,
        sum(meta$study_group == "GSE348821") == 2L,
        nrow(selected) == 17L, sum(selected$design == "KO") == 1L)
    # Exhaustive size-matched baseline check, avoiding RNG/version dependence.
    for (i in seq_len(nrow(mod))) {
        candidates <- strsplit(mod$candidate_set[i], ",", fixed = TRUE)[[1L]]
        subsets <- combn(candidates[-1L], mod$module_size[i] - 1L, simplify = FALSE)
        p <- mean(vapply(subsets, function(s)
            mod$expected_label[i] %in% c(candidates[1L], s), logical(1)))
        stopifnot(abs(p - mod$matched_size_random_expected[i]) < 1e-12)
    }
    message("Self-tests passed: full denominators, study groups, baseline, exact random subsets.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
artifacts <- list(ranks = ranks, full_rankings = full, metrics = summary, grouped_macro = macro,
    modules = mod, module_summary = module_summary, historical_verification = hist,
    historical_summary = historical_summary)
for (nm in names(artifacts)) write_tsv(artifacts[[nm]], file.path(outdir, paste0(nm, ".tsv")))
write_tsv(data.frame(file = "R/mechanism_discover.R", md5 = unname(tools::md5sum("R/mechanism_discover.R")),
    baseline_commit = "978eabaee502350e12e5b6a3988c54091ac79fe3"), file.path(outdir, "implementation.tsv"))
p <- sort(list.files(outdir, full.names = TRUE)); p <- p[basename(p) != "output-manifest.tsv"]
write_tsv(manifest(p, outdir), file.path(outdir, "output-manifest.tsv"))
print(summary[summary$group_type %in% c("all", "design"), ], row.names = FALSE)
print(module_summary, row.names = FALSE)
print(historical_summary, row.names = FALSE)
cat("Outputs:", normalizePath(outdir), "\n")
