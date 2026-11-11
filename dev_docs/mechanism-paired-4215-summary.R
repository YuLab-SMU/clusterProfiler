#!/usr/bin/env Rscript
# Offline descriptive endpoints; never reconstruct a missing historical native ledger.
paired_summarize <- function(root, out = file.path(root, "analysis")) {
    source(file.path(root, "replay-runner.R"), local = TRUE)
    paired_check(root, "frozen-manifest.tsv")
    paired_check(root, "raw-manifest.tsv")
    paired_check(root, "replay/manifest.tsv")
    if (dir.exists(out)) stop("Analysis output exists; refusing overwrite: ", out)
    read_tab <- function(name) read.delim(file.path(root, "replay", paste0(name, ".tsv")),
        stringsAsFactors = FALSE, check.names = FALSE)
    coverage <- read_tab("coverage"); rows <- read_tab("rows"); traces <- read_tab("traces")
    trans <- read_tab("transitions")
    meta <- read.delim(file.path(root, "inputs", "metadata.tsv"), stringsAsFactors = FALSE)
    stopifnot(nrow(coverage) == 4L * nrow(meta), !anyDuplicated(coverage[c("profile", "code", "cap")]))
    selected_state <- function(z) {
        finite <- is.finite(z$es) & is.finite(z$padj) & z$es != 0
        significant <- finite & z$padj < 0.05
        opposite <- (z$direction %in% "up" & z$es < 0) | (z$direction %in% "down" & z$es > 0)
        known <- z$direction %in% c("up", "down")
        list(finite = finite, support = significant & known & opposite,
             conflict = significant & known & !opposite)
    }
    extra <- lapply(seq_len(nrow(coverage)), function(i) {
        a <- coverage[i, ]
        z <- rows[rows$profile == a$profile & rows$code == a$code & rows$cap == a$cap, , drop = FALSE]
        st <- selected_state(z)
        data.frame(n_A = sum(z$level %in% "A"), n_B = sum(z$level %in% "B"),
            n_C = sum(z$level %in% "C"), any_AB = any(z$level %in% c("A", "B")),
            n_selected_finite = sum(st$finite), n_selected_support = sum(st$support, na.rm = TRUE),
            n_selected_conflict = sum(st$conflict, na.rm = TRUE),
            any_selected_support = any(st$support, na.rm = TRUE),
            any_selected_conflict = any(st$conflict, na.rm = TRUE))
    })
    profiles <- cbind(coverage, do.call(rbind, extra))
    # Unique raw requests, independent of whether either implementation retained a ledger.
    raw_paths <- sort(list.files(file.path(root, "raw"), full.names = TRUE))
    raw <- do.call(rbind, lapply(raw_paths, function(path) {
        x <- readRDS(path)
        d <- if (isTRUE(x$ok)) as.data.frame(x$value[[1]]) else data.frame()
        finite <- x$kind == "test" && nrow(d) > 0L && all(c("es", "padj") %in% names(d)) &&
            is.finite(d$es[1]) && is.finite(d$padj[1]) && d$es[1] != 0
        data.frame(key = x$key, kind = x$kind, ok = isTRUE(x$ok),
            response_status = if (!isTRUE(x$ok)) "failed" else if (x$kind == "lookup") {
                if (nrow(d)) "lookup-rows" else "lookup-empty"
            } else if (finite) "finite-test-statistics" else "noncomputable-test-statistics",
            n_response_rows = nrow(d), stringsAsFactors = FALSE)
    }))
    tr <- merge(traces, raw, by = c("key", "kind", "ok"), all.x = TRUE, sort = FALSE)
    stopifnot(!anyNA(tr$response_status))
    summary <- lapply(seq_len(nrow(paired_conditions)), function(i) {
        a <- paired_conditions[i, ]
        z <- profiles[profiles$code == a$code & profiles$cap == a$cap, , drop = FALSE]
        q <- tr[tr$code == a$code & tr$cap == a$cap, , drop = FALSE]
        q <- q[!duplicated(q$key), , drop = FALSE]
        # Native ledger counts are absent, not zero, for the historical implementation.
        total <- function(field) if (all(is.na(z[[field]]))) NA_integer_ else sum(z[[field]], na.rm = TRUE)
        data.frame(code = a$code, cap = a$cap, n_profiles = nrow(z),
            n_graded = sum(z$status == "graded"), n_no_enrichment = sum(z$status == "no-enrichment"),
            n_implementation_error = sum(z$status == "implementation-error"),
            n_pathways = total("n_pathways"), n_A = total("n_A"), n_B = total("n_B"), n_C = total("n_C"),
            profiles_any_AB = sum(z$any_AB), profiles_selected_support = sum(z$any_selected_support),
            profiles_selected_conflict = sum(z$any_selected_conflict),
            n_selected_finite = total("n_selected_finite"), n_selected_support = total("n_selected_support"),
            n_selected_conflict = total("n_selected_conflict"),
            native_ledger_available = all(z$ledger_available), native_outcomes = total("n_ledger"),
            native_support = total("n_support"), native_conflict = total("n_conflict"),
            native_NS = total("n_NS"), native_NC = total("n_NC"), native_failed = total("n_failed"),
            profiles_native_support = if (all(is.na(z$n_support))) NA_integer_ else sum(z$n_support > 0, na.rm = TRUE),
            profiles_native_conflict = if (all(is.na(z$n_conflict))) NA_integer_ else sum(z$n_conflict > 0, na.rm = TRUE),
            profiles_native_finite = if (all(is.na(z$n_finite))) NA_integer_ else sum(z$n_finite > 0, na.rm = TRUE),
            unique_lookup_requests = sum(q$kind == "lookup"), unique_test_requests = sum(q$kind == "test"),
            unique_request_errors = sum(!q$ok),
            unique_finite_test_requests = sum(q$response_status == "finite-test-statistics"),
            unique_noncomputable_test_requests = sum(q$response_status == "noncomputable-test-statistics"))
    })
    summary <- do.call(rbind, summary)
    equal <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    pathway_pairs <- lapply(split(trans, interaction(trans$contrast, trans$held_value, drop = TRUE)), function(z) {
        common <- !is.na(z$level_from) & !is.na(z$level_to)
        data.frame(contrast = z$contrast[1], held_value = z$held_value[1], from = z$from[1], to = z$to[1],
            common_pathways = sum(common), removed_pathways = sum(!is.na(z$level_from) & is.na(z$level_to)),
            added_pathways = sum(is.na(z$level_from) & !is.na(z$level_to)),
            changed_grade = sum(common & !equal(z$level_from, z$level_to)),
            changed_verdict = sum(common & !equal(z$verdict_from, z$verdict_to)),
            changed_reason = sum(common & !equal(z$reason_from, z$reason_to)))
    })
    pathway_pairs <- do.call(rbind, pathway_pairs); rownames(pathway_pairs) <- NULL
    contrasts <- list(c("code", "8-4", "baseline", "fixed"), c("code", "5-2", "baseline", "fixed"),
                      c("cap", "baseline", "8-4", "5-2"), c("cap", "fixed", "8-4", "5-2"))
    fields <- c("profile", "expected_driver", "status", "any_AB", "any_selected_support", "any_selected_conflict",
                "n_pathways", "n_A", "n_B", "n_C", "n_support", "n_conflict", "n_NC", "n_failed")
    profile_pairs <- do.call(rbind, lapply(contrasts, function(p) {
        variable <- p[1]; held <- if (variable == "code") "cap" else "code"
        z <- profiles[profiles[[held]] == p[2], ]
        a <- z[z[[variable]] == p[3], fields]; b <- z[z[[variable]] == p[4], fields]
        pair <- merge(a, b, by = c("profile", "expected_driver"), all = TRUE, sort = TRUE,
                      suffixes = c("_from", "_to"))
        stopifnot(nrow(pair) == nrow(meta))
        cbind(contrast = p[1], held_value = p[2], from = p[3], to = p[4], pair)
    }))
    selected_pairs <- do.call(rbind, lapply(contrasts, function(p) {
        variable <- p[1]; held <- if (variable == "code") "cap" else "code"
        z <- rows[rows[[held]] == p[2], ]
        fields <- c("profile", "expected_driver", "pathway", "driver", "dataset", "cell_line", "es", "padj",
                    "direction", "support", "level", "verdict", "reason")
        a <- z[z[[variable]] == p[3], fields]; b <- z[z[[variable]] == p[4], fields]
        pair <- merge(a, b, by = fields[1:3], all = TRUE, sort = TRUE, suffixes = c("_from", "_to"))
        common <- !is.na(pair$level_from) & !is.na(pair$level_to)
        pair$common <- common
        for (field in c("driver", "dataset", "cell_line", "es", "padj", "direction", "support"))
            pair[[paste0(field, "_changed")]] <- common & !equal(pair[[paste0(field, "_from")]], pair[[paste0(field, "_to")]])
        cbind(contrast = p[1], held_value = p[2], from = p[3], to = p[4], pair)
    }))
    dir.create(out, recursive = TRUE)
    tables <- list("profile-endpoints" = profiles, "condition-summary" = summary,
                   "unique-raw-requests" = raw, "paired-pathway-summary" = pathway_pairs,
                   "paired-profile-endpoints" = profile_pairs, "paired-selected-transitions" = selected_pairs)
    for (name in names(tables)) paired_tsv(tables[[name]], file.path(out, paste0(name, ".tsv")))
    writeLines(c("Selected finite and unique finite-test counts require finite ES/padj and nonzero ES; native n_finite in replay coverage only checks finite values.",
        "Selected-row support/conflict uses only the returned row's finite ES, adjusted p < 0.05 and known direction.",
        "It is NOT a reconstruction of the historical full native ledger or all tested datasets.",
        "No selected significant record does not mean no usable test; use raw-query coverage separately.",
        "Unique raw-query counts are not pathway outcomes or independent biological replicates.",
        "Native historical ledger fields remain NA, not zero. All 17 profiles remain in each denominator.",
        "Code contrasts include dataset lookup/selection policies; cap contrasts hold code fixed.",
        "No endpoint here establishes improved accuracy or an upstream causal driver."), file.path(out, "endpoint-definitions.txt"))
    paired_tsv(paired_manifest(out, sort(list.files(out, full.names = TRUE))), file.path(out, "manifest.tsv"))
    invisible(summary)
}
if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    root <- if (length(args)) args[1] else "../adjudicate-interpret/results/mechanism-paired-4215"
    out <- if (length(args) > 1L) args[2] else file.path(root, "analysis")
    print(paired_summarize(root, out))
}
