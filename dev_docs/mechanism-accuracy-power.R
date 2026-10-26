suppressMessages({devtools::load_all(".", quiet=TRUE)})
options(clusterProfiler_cache_dir="../biomedical-knowledge-mining-book/cache",
        clusterProfiler.progress="none")
eval <- readRDS("/tmp/cases/eval.rds")
g <- mechanism_sets("tf"); m <- g@gsid2gene
sym <- mechanism_known_symbols("Homo sapiens"); N <- length(sym)

## 真值靶集在输入中的观测重叠 vs 背景期望
rows <- list()
for (i in seq_len(nrow(eval))) {
  x <- eval$truth[i]; id <- eval$dataset[i]
  d <- as.data.frame(gpsadb::gpsadb_de(id, limit=500)$data)
  gs <- unique(d$gene[!is.na(d$padj) & d$padj < 0.05])
  reg <- mechanism_alias_regulator(x)
  terms <- unique(m$gsid[mechanism_regulator(m$gsid, symbols=sym) == reg])
  setg <- unique(m$gene[m$gsid %in% terms])
  setg <- intersect(setg, sym)
  if (!length(setg)) next
  obs <- length(intersect(setg, gs))
  exp <- length(gs) * length(setg) / N
  rows[[length(rows)+1]] <- data.frame(truth=x, n_input=length(gs), set_size=length(setg),
      observed=obs, expected=exp, ratio=obs/exp, rank=eval$rank[i], stringsAsFactors=FALSE)
}
r <- do.call(rbind, rows)
cat("=== 真值靶集的富集情况（观测/期望）===\n")
cat("n =", nrow(r), "\n")
cat("观测 > 期望 的比例:", sprintf("%.0f%%", 100*mean(r$ratio > 1)), "\n")
cat("比值中位数:", sprintf("%.2f", median(r$ratio)), "\n")
cat("\n按是否找回分组:\n")
r$found <- !is.na(r$rank)
for (f in c(TRUE, FALSE)) {
  s <- r[r$found==f,]
  if (nrow(s)) cat(sprintf("  %s n=%2d  比值中位 %.2f  观测中位 %.0f  期望中位 %.0f\n",
      if(f) "找回" else "未找回", nrow(s), median(s$ratio), median(s$observed), median(s$expected)))
}
cat("\n=== 结论 ===\n")
cat("若比值中位数 < 1，说明扰动后的差异基因并不富含该因子的直接靶标\n")
saveRDS(r, "/tmp/cases/ratio.rds")
