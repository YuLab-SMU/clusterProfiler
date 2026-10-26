suppressMessages({devtools::load_all(".", quiet=TRUE)})
options(clusterProfiler_cache_dir="../biomedical-knowledge-mining-book/cache",
        clusterProfiler.progress="none")
truth <- readRDS("/tmp/cases/truth.rds")
sym <- mechanism_known_symbols("Homo sapiens")
set.seed(42)
pick <- sample(names(truth), min(30, length(truth)))
res <- list()
for (x in pick) {
  id <- truth[[x]]
  d <- tryCatch(gpsadb::gpsadb_de(id, limit=500), error=function(e) NULL)
  if (is.null(d) || is.null(d$data)) next
  dd <- as.data.frame(d$data)
  dd <- dd[!is.na(dd$padj) & dd$padj < 0.05, ]
  gs <- unique(dd$gene); gs <- gs[!is.na(gs) & nzchar(gs)]
  if (length(gs) < 30) next
  ev <- tryCatch(suppressWarnings(evidence(gs, upstream="tf")), error=function(e) NULL)
  if (is.null(ev)) next
  r <- ev$regulators
  r <- r[mechanism_significant(r), ]
  r <- r[order(r$padj_within_source), ]
  regs <- unique(r$regulator)
  rank <- match(mechanism_alias_regulator(x), regs)
  res[[x]] <- data.frame(
    truth = x, dataset = id, n_input = length(gs), n_sig_entities = length(regs),
    rank = if (is.na(rank)) NA_integer_ else rank,
    found = !is.na(rank), stringsAsFactors = FALSE)
  cat(sprintf("%-10s %-8s 输入%4d 显著实体%4d 真值名次 %s\n", x, id, length(gs), length(regs),
      if (is.na(rank)) "未出现" else rank))
}
out <- do.call(rbind, res)
saveRDS(out, "/tmp/cases/eval.rds")
cat("\n=== 汇总 ===\n")
cat("测试数据集:", nrow(out), "\n")
cat("真值被找回(出现在显著实体中):", sum(out$found), sprintf("(%.0f%%)\n", 100*mean(out$found)))
for (k in c(1,3,5,10,20)) cat(sprintf("  recall@%d: %.0f%%\n", k, 100*mean(!is.na(out$rank) & out$rank<=k)))
cat("MRR:", sprintf("%.3f", mean(ifelse(is.na(out$rank), 0, 1/out$rank))), "\n")
cat("真值名次中位数:", median(out$rank, na.rm=TRUE), "\n")
