## 共用 fixture。testthat 只跨文件共享 helper-*.R，所以这些定义不能留在
## 任何一个 test-*.R 里：另一个测试文件会看不到它们。

## ---- fixtures -------------------------------------------------------------

write_gmt <- function(path, sets) {
    lines <- vapply(names(sets), function(nm) {
        paste(c(nm, "https://example.invalid/" , sets[[nm]]), collapse = "\t")
    }, "")
    writeLines(lines, path)
    path
}

## A tiny collection shaped like C3:TFT: one GTRD set, one motif set naming a
## factor, and one unnamed motif set.
tft_fixture <- function() {
    list(
        FOXO3_TARGET_GENES = paste0("g", 1:8),
        TGANTCA_AP1_C      = paste0("g", 1:6),
        AACTTT_UNKNOWN     = paste0("g", 1:5)
    )
}

## Not withr::local_tempdir(): its cleanup runs when *this* function returns,
## so the directory is gone by the time the test uses it. Make it here and let
## the caller decide when it dies.
evidence_fixture <- function() {
    dir <- tempfile("mechanism")
    dir.create(dir)
    write_gmt(file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt"), tft_fixture())
    writeLines(c("release=2026.2.Hs", "accessed_date=2026-10-25"),
               file.path(dir, "c3.tft.v2026.2.Hs.symbols.gmt.meta"))
    dir
}

use_fixture_cache <- function(env = parent.frame()) {
    dir <- evidence_fixture()
    withr::defer(unlink(dir, recursive = TRUE), envir = env)
    old <- options(clusterProfiler_cache_dir = dir)
    withr::defer(options(old), envir = env)
    dir
}
