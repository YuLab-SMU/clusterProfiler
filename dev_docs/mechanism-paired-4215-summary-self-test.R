#!/usr/bin/env Rscript
# Offline integration self-check against a complete snapshot; writes only a new temp directory.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "../adjudicate-interpret/results/mechanism-paired-4215-paced"
source("dev_docs/mechanism-paired-4215-summary.R")
out <- tempfile("paired-analysis-self-test-")
summary <- paired_summarize(root, out)
e <- new.env(); sys.source(file.path(root, "replay-runner.R"), e)
e$paired_check(out, "manifest.tsv")
read_tab <- function(name) read.delim(file.path(out, paste0(name, ".tsv")), stringsAsFactors = FALSE)
p <- read_tab("profile-endpoints")
r <- read.delim(file.path(root, "replay", "rows.tsv"), stringsAsFactors = FALSE)
m <- read.delim(file.path(root, "inputs", "metadata.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(summary) == 4L, all(summary$n_profiles == nrow(m)), nrow(p) == 4L * nrow(m))
for (i in seq_len(nrow(p))) {
    z <- r[r$profile == p$profile[i] & r$code == p$code[i] & r$cap == p$cap[i], ]
    support <- is.finite(z$es) & is.finite(z$padj) & z$padj < .05 &
        ((z$direction %in% "up" & z$es < 0) | (z$direction %in% "down" & z$es > 0))
    conflict <- is.finite(z$es) & is.finite(z$padj) & z$padj < .05 &
        ((z$direction %in% "up" & z$es > 0) | (z$direction %in% "down" & z$es < 0))
    stopifnot(p$n_selected_support[i] == sum(support, na.rm = TRUE),
        p$n_selected_conflict[i] == sum(conflict, na.rm = TRUE), p$n_pathways[i] == nrow(z))
}
stopifnot(all(is.na(summary$native_outcomes[summary$code == "baseline"])))
pairs <- read_tab("paired-profile-endpoints")
stopifnot(nrow(pairs) == 4L * nrow(m))
raw <- read_tab("unique-raw-requests")
stopifnot(!anyDuplicated(raw$key), nrow(raw) == length(list.files(file.path(root, "raw"))))
cat("Summary self-test passed: external-output manifest, all-profile denominators, native NA and selected endpoint reconciliation.\n")
