#!/usr/bin/env Rscript
# Tiny synthetic, OFFLINE fixture; never writes the observational result root.
source("dev_docs/mechanism-paired-4215.R")
stopifnot(requireNamespace("clusterProfiler"), requireNamespace("gpsadb"))
root <- tempfile("paired-selftest-")
# Leave fixture artifacts in the reported directory for inspection; no destructive cleanup.
origin <- tempfile(fileext = ".rds"); saveRDS("synthetic fixture, not observations", origin)
genes <- LETTERS[1:8]
ego <- data.frame(ID = paste0("P", 1:8), Description = paste0("path", 1:8),
    geneID = vapply(1:8, function(i) paste(c(genes[i], "X", "Y", "Z", "W"), collapse = "/"), ""),
    pathway_direction = "up")
profiles <- list(list(label = "tiny", signed = setNames(rep(1, 12), c(genes, "X", "Y", "Z", "W")), ego = ego),
    list(label = "empty", signed = c(A = 1), ego = NULL),
    list(label = "lookup failure", signed = setNames(rep(1, 12), c(genes, "X", "Y", "Z", "W")), ego = ego[1,,drop=FALSE]))
meta <- data.frame(dataset = c("tiny", "empty", "lookup failure"), expected_label = c("FIXTURE", "FIXTURE", "FAIL"))
paired_snapshot(root, profiles, meta, origin)
# Restore namespace function hashing comparison: no scoped source may alter package functions.
original <- paired_hash(body(get("pathway_evidence", asNamespace("clusterProfiler"))))
network_calls <- 0L
fixture <- function(kind, args) {
    network_calls <<- network_calls + 1L
    if (kind == "lookup") {
        if (args$perturbed_gene == "FAIL") return(list(ok = FALSE, error = "fixture lookup transport failure"))
        return(list(ok = TRUE, value = list(data.frame(dataset_id = paste0("DS", 1:4),
            cell_line = paste0("CELL", 1:4), n_samples = 4:1, in_reference = TRUE))))
    }
    i <- match(intersect(args$genes, genes), genes)
    if (i == 5L) return(list(ok = FALSE, error = "fixture test service error"))
    value <- if (i == 3L) data.frame(es = NA_real_, padj = NA_real_) else
        data.frame(es = if (i == 4L) 1 else -1, padj = if (i == 2L) 0.8 else 0.001)
    list(ok = TRUE, value = list(value))
}
state <- paired_requester(root, live = fixture)
collected <- paired_run(root, state)
stopifnot(nrow(collected$coverage) == 12L, length(unique(collected$coverage$code)) == 2L,
    length(unique(collected$coverage$cap)) == 2L,
    all(collected$coverage$status[collected$coverage$profile == "empty"] == "no-enrichment"),
    all(c("supports", "conflicts", "not-computable", "not-significant", "failed") %in% collected$ledger$status),
    all(is.na(collected$coverage$n_ledger[collected$coverage$code == "baseline"])),
    all(collected$coverage$ledger_available[collected$coverage$code == "fixed"]),
    any(collected$rows$reason[collected$rows$code == "fixed"] == "lookup-failed", na.rm = TRUE),
    all(collected$coverage$observed_lookup_errors[collected$coverage$profile == "lookup failure"] > 0),
    sum(!collected$traces$cache_hit) == network_calls,
    any(collected$traces$cache_hit),
    length(list.files(file.path(root, "raw"))) == network_calls,
    identical(original, paired_hash(body(get("pathway_evidence", asNamespace("clusterProfiler"))))))
paired_check(root, "frozen-manifest.tsv")
replayed <- paired_run(root, paired_requester(root))
stopifnot(identical(collected$rows, replayed$rows), identical(collected$ledger, replayed$ledger),
          identical(collected$coverage, replayed$coverage), all(replayed$traces$cache_hit))
# Same canonical key across implementations and caps is reused, not recollected.
stopifnot(any(duplicated(replayed$traces$key)))
# A missing key must fail outside the historical implementation's catch-and-NULL.
e <- paired_source(root, "baseline", paired_requester(root)$request)
missing_state <- paired_requester(root)
missing_env <- paired_source(root, "baseline", missing_state$request)
missing <- tryCatch(paired_run(root, {
    s <- paired_requester(root)
    original_request <- s$request
    s$request <- function(kind, args) {
        if (kind == "lookup") args$perturbed_gene <- "UNCACHED"
        original_request(kind, args)
    }
    s
}), error = identity)
stopifnot(inherits(missing, "error"), grepl("CACHE MISS", conditionMessage(missing)))
# Enforce budget even if source catches the wrapper's error.
other <- tempfile("paired-budget-")
paired_snapshot(other, profiles, meta, origin)
budget <- tryCatch(paired_run(other, paired_requester(other, fixture, max_new = 0L)), error = identity)
stopifnot(inherits(budget, "error"), grepl("budget exhausted", conditionMessage(budget)))
trans <- paired_transitions(replayed$rows)
stopifnot(all(c("code", "cap") %in% trans$contrast), all(c("level_from", "verdict_to", "reason_to") %in% names(trans)))
# Malformed raw records must escape historical catch-and-NULL and never receive a seal.
key <- collected$traces$key[1L]
for (tamper in c("key", "unreadable", "kind", "args", "ok", "missing-value", "missing-error")) {
    badroot <- tempfile(paste0("paired-invalid-", tamper, "-")); dir.create(badroot)
    stopifnot(all(file.copy(list.files(root, full.names = TRUE), badroot, recursive = TRUE)))
    path <- file.path(badroot, "raw", paste0(key, ".rds"))
    rec <- readRDS(path)
    if (tamper == "key") rec$key <- "wrong"
    if (tamper == "kind") rec$kind <- "wrong"
    if (tamper == "args") rec$args <- list(wrong = TRUE)
    if (tamper == "ok") rec$ok <- c(TRUE, FALSE)
    if (tamper == "missing-value") rec$value <- NULL
    if (tamper == "missing-error") { rec$ok <- FALSE; rec$error <- NULL }
    if (tamper == "unreadable") writeLines("not an RDS file", path) else saveRDS(rec, path)
    failed <- tryCatch(paired_main(c("finalize", badroot)), error = identity)
    stopifnot(inherits(failed, "error"), grepl("CACHE INVALID", conditionMessage(failed)),
              !file.exists(file.path(badroot, "raw-manifest.tsv")))
}
stopifnot("inputs/plan.md" %in% read.delim(file.path(root, "frozen-manifest.tsv"))$file,
    identical(readLines(file.path(root, "inputs", "plan.md")), readLines("dev_docs/mechanism-paired-4215-plan.md")))
# Exercise the HTTP-attempt pacing primitive without network access or real sleeping.
clock <- 0; sleeps <- numeric(); starts <- numeric()
pace <- paired_pacer(2, now = function() clock, sleep = function(s) {
    sleeps <<- c(sleeps, s); clock <<- clock + s
})
for (elapsed in c(0, 0.3, 0.4, 3, 0)) {
    clock <- clock + elapsed; pace(); starts <- c(starts, clock)
}
stopifnot(all(diff(starts) >= 2), identical(starts, c(0, 2, 4, 7, 9)),
          length(sleeps) == 3L)
# The captured runner must finalize/replay from a fresh process outside the checkout.
oldwd <- getwd(); setwd(tempdir())
for (mode in c("finalize", "replay")) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
        c(shQuote(file.path(root, "replay-runner.R")), mode, shQuote(root)))
    stopifnot(status == 0L)
}
setwd(oldwd)
stopifnot(file.exists(file.path(root, "replay", "full-results.rds")))
paired_check(root, "raw-manifest.tsv")
cat("PASS: all four conditions; fresh-process/fresh-cwd archive replay; full raw reuse; offline identical rows/ledger/coverage; cache-miss and budget guards;\n",
    "NC/NS/conflict/test failure and lookup failure distinct; empty profiles kept; production namespace unchanged.\n",
    "Synthetic fixture root: ", root, "\n", sep = "")
