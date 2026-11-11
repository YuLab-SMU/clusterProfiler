#!/usr/bin/env Rscript
# Paired, expected-driver-only experiment. Run from clusterProfiler with R_LIBS_USER=~/Rlib.
# freeze -> collect (LIVE, explicit) -> finalize -> replay (strictly offline).
# This compares whole source implementations, INCLUDING dataset-selection policies.
# Cap contrasts hold code fixed. Selected expected-driver preconditions are not accuracy.
paired_refs <- c(baseline = "978eabaee502350e12e5b6a3988c54091ac79fe3", fixed = "fe7c060")
paired_conditions <- expand.grid(code = names(paired_refs), cap = c("8-4", "5-2"),
                                 stringsAsFactors = FALSE)
paired_tsv <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
    quote = TRUE, na = "NA", fileEncoding = "UTF-8")
paired_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
paired_manifest <- function(root, paths) data.frame(file = substring(paths, nchar(root) + 2L),
    bytes = unname(file.info(paths)$size), md5 = unname(tools::md5sum(paths)),
    sha256 = unname(vapply(paths, digest::digest, "", algo = "sha256", file = TRUE)), row.names = NULL)
paired_check <- function(root, name) {
    m <- read.delim(file.path(root, name), check.names = FALSE)
    paths <- file.path(root, m$file)
    if (!all(file.exists(paths)) || !isTRUE(all.equal(paired_manifest(root, paths), m, check.attributes = FALSE)))
        stop("Checksum verification failed: ", name)
    invisible(m)
}
paired_once <- function(x, path) {
    if (file.exists(path)) stop("Refusing to overwrite immutable file: ", path)
    pending <- tempfile(".pending-", tmpdir = dirname(path))
    saveRDS(x, pending, version = 2)
    if (!file.rename(pending, path)) stop("Cannot atomically finalize: ", path)
}
paired_git <- function(ref, file) {
    x <- system2("git", c("show", shQuote(paste0(ref, ":", file))), stdout = TRUE)
    if (!is.null(attr(x, "status"))) stop("Cannot read pinned source: ", ref, ":", file)
    x
}
paired_transform <- function(expr) {
    if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return(expr)
    if (is.call(expr) && identical(expr[[1]], as.name("::")) &&
        identical(expr[[2]], as.name("gpsadb"))) {
        n <- as.character(expr[[3]])
        if (!n %in% c("gpsadb_search_datasets", "gpsadb_test_gene_set"))
            stop("Unexpected GPSAdb namespace call: ", n)
        return(as.name(paste0(".paired_", n)))
    }
    # Preserve missing arguments in function formals and indexing calls.
    was_pairlist <- is.pairlist(expr)
    for (i in seq_along(expr)) if (!identical(expr[[i]], quote(expr = )))
        expr[i] <- list(paired_transform(expr[[i]]))
    if (was_pairlist) as.pairlist(expr) else expr
}
paired_source <- function(root, code, request) {
    e <- new.env(parent = asNamespace("clusterProfiler"))
    deps <- readRDS(file.path(root, "inputs", "helpers.rds"))
    for (n in names(deps)) {
        f <- deps[[n]]
        if (is.function(f)) environment(f) <- e
        assign(n, f, e)
    }
    e$.paired_gpsadb_search_datasets <- function(...) request("lookup", list(...))
    e$.paired_gpsadb_test_gene_set <- function(dataset_id, genes)
        request("test", list(dataset_id = dataset_id, genes = genes))
    original <- parse(file.path(root, "sources", paste0(code, ".R")), keep.source = FALSE)
    eval(paired_transform(original), e)
    e
}
paired_snapshot <- function(root, profiles, metadata, original_paths) {
    if (dir.exists(root)) stop("Result root exists; freeze must never overwrite: ", root)
    dir.create(file.path(root, "inputs"), recursive = TRUE)
    dir.create(file.path(root, "sources")); dir.create(file.path(root, "raw"))
    stopifnot(length(profiles) == nrow(metadata),
        identical(vapply(profiles, `[[`, "", "label"), metadata$dataset))
    # Keep full original profile objects, including ego S4 where present and signed FC.
    paired_once(profiles, file.path(root, "inputs", "profiles.rds"))
    paired_tsv(metadata, file.path(root, "inputs", "metadata.tsv"))
    paired_tsv(paired_manifest(dirname(original_paths[1]), original_paths),
               file.path(root, "inputs", "origin-manifest.tsv"))
    refs <- vapply(paired_refs, function(r) {
        z <- system2("git", c("rev-parse", r), stdout = TRUE)
        stopifnot(is.null(attr(z, "status")), length(z) == 1L); z
    }, "")
    helpers <- list(); ns <- asNamespace("clusterProfiler")
    for (code in names(refs)) {
        txt <- paired_git(refs[[code]], "R/pathway_evidence.R")
        writeLines(txt, file.path(root, "sources", paste0(code, ".R")), useBytes = TRUE)
        expr <- parse(text = txt, keep.source = FALSE)
        transformed <- paired_transform(expr)
        writeLines(vapply(transformed, function(z) paste(deparse(z, width.cutoff = 120L), collapse = "\n"), ""),
                   file.path(root, "sources", paste0(code, "-injected.R")), useBytes = TRUE)
        e <- new.env(parent = ns); eval(expr, e)
        queue <- ls(e, all.names = TRUE); local <- queue; visited <- character()
        while (length(queue)) {
            n <- queue[1]; queue <- queue[-1]
            if (n %in% visited) next
            visited <- c(visited, n)
            f <- if (n %in% local) get(n, e) else get(n, ns)
            if (!is.function(f)) next
            globals <- codetools::findGlobals(f, merge = FALSE)$functions
            more <- globals[!globals %in% c(local, names(helpers)) &
                vapply(globals, exists, TRUE, envir = ns, inherits = FALSE)]
            for (g in more) {
                val <- get(g, ns)
                if (is.function(val)) {
                    # Store bodies/formals only: no serialized namespace or process state.
                    environment(val) <- baseenv(); helpers[[g]] <- val
                    queue <- c(queue, g)
                }
            }
        }
    }
    stopifnot(file.copy("dev_docs/mechanism-paired-4215-plan.md", file.path(root, "inputs", "plan.md")))
    stopifnot(file.copy("dev_docs/mechanism-paired-4215.R", file.path(root, "replay-runner.R")))
    stopifnot(file.copy("dev_docs/mechanism-paired-4215-self-test.R", file.path(root, "fixture-self-test.R")))
    paired_once(list(R = R.version, session = utils::sessionInfo()), file.path(root, "inputs", "R-environment.rds"))
    paired_once(helpers, file.path(root, "inputs", "helpers.rds"))
    paired_tsv(data.frame(helper = names(helpers), sha256 = vapply(helpers, paired_hash, "")),
               file.path(root, "inputs", "helpers.tsv"))
    pkgs <- unique(c("clusterProfiler", "gpsadb", "aisdk.mcp", "httr2", "jsonlite", "digest",
                     "DOSE", "AnnotationDbi", "org.Hs.eg.db", loadedNamespaces()))
    versions <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
    paired_once(versions, file.path(root, "inputs", "package-versions.rds"))
    paired_once(list(refs = refs, created_utc = format(Sys.time(), tz = "UTC"),
        conditions = paired_conditions, scope = "All17 expected drivers only; no rank tuning; no accuracy claim",
        policy = "Same full raw cache; code includes lookup/selection policy changes; cap contrast same code",
        transform = "Only AST gpsadb::gpsadb_search_datasets and gpsadb::gpsadb_test_gene_set become .paired_*; all other source unchanged",
        helpers = "Recursively captured clusterProfiler namespace function dependencies; bodies/checksums frozen; external package versions checked",
        rate_policy = "At least 2 seconds between starts of ALL HTTP attempts, including MCP initialization/notifications and proxy retry; no HTTP status retries; prior roots never reused",
        request_keys = "SHA256 of kind and exact named arguments, including original genes order/multiplicity; no canonicalization",
        disclosure = "Baseline native ledger/lookup unavailable (NA), never reconstructed as absent evidence; compare selected rows for code effects",
        R_version = R.version.string,
        raw = "Complete returned gpsadb object + complete tool payload + full HTTP request/response bodies/headers/errors; failures cached once"),
        file.path(root, "inputs", "design.rds"))
    paths <- sort(list.files(root, recursive = TRUE, full.names = TRUE))
    paired_tsv(paired_manifest(root, paths), file.path(root, "frozen-manifest.tsv"))
}
paired_freeze <- function(root) {
    paths <- c("../adjudicate-interpret/results/all11_de.rds",
               "../adjudicate-interpret/results/p3/p3_new_de.rds",
               "../adjudicate-interpret/results/audit-benchmark-4215/inputs/metadata.tsv")
    stopifnot(all(file.exists(paths)))
    profiles <- c(readRDS(paths[1]), readRDS(paths[2]))
    meta <- read.delim(paths[3], stringsAsFactors = FALSE)
    stopifnot(length(profiles) == 17L, nrow(meta) == 17L, !anyDuplicated(meta$dataset))
    paired_snapshot(root, profiles, meta, paths)
}
# Rate-limit HTTP attempts, not logical queries; the clock is monotonic within a process.
paired_pacer <- function(interval = 2, now = function() unname(proc.time()["elapsed"]), sleep = Sys.sleep) {
    last <- -Inf
    function() {
        wait <- max(0, interval - (now() - last))
        if (wait > 0) sleep(wait)
        last <<- now()
        wait
    }
}
# No gpsadb namespace modifications: original result adapters cloned into one local environment.
# Custom MCP transport disables the SDK's implicit three-attempt HTTP retry.
paired_live <- function() {
    state <- new.env(parent = emptyenv()); state$proxy <- FALSE; state$client <- NULL
    state$wire <- list(); state$payload <- NULL
    pace <- paired_pacer(2)
    url <- get("gpsadb_url", asNamespace("gpsadb"))()
    transport <- function(message) {
        routes <- if (state$proxy) "proxy" else c("direct", "proxy")
        for (route in routes) {
            started <- format(Sys.time(), tz = "UTC")
            body <- jsonlite::toJSON(message, auto_unbox = TRUE, null = "null")
            req <- httr2::request(url)
            req <- httr2::req_headers(req, Accept = "application/json, text/event-stream",
                `Content-Type` = "application/json", `MCP-Protocol-Version` = "2025-06-18")
            if (!is.null(state$sid)) req <- httr2::req_headers(req, `Mcp-Session-Id` = state$sid)
            req <- httr2::req_body_raw(req, body, type = "application/json")
            req <- httr2::req_timeout(req, 25)
            req <- httr2::req_options(req, proxy = if (route == "proxy") "http://127.0.0.1:7897" else "",
                                      noproxy = if (route == "proxy") "" else "*")
            req <- httr2::req_error(req, is_error = function(resp) FALSE)
            waited <- pace()
            started <- format(Sys.time(), tz = "UTC")
            got <- tryCatch(httr2::req_perform(req), error = function(e) e)
            record <- list(url = url, route = route, pacing_wait_seconds = waited,
                min_http_interval_seconds = 2, started_utc = started,
                finished_utc = format(Sys.time(), tz = "UTC"), request = message,
                request_body = body, session_id = state$sid)
            if (inherits(got, "error")) {
                record$error <- conditionMessage(got)
                state$wire[[length(state$wire) + 1L]] <- record
                if (route == "direct") { state$proxy <- TRUE; next }
                stop(got)
            }
            record$status <- httr2::resp_status(got)
            record$headers <- httr2::resp_headers(got)
            record$body <- httr2::resp_body_raw(got)
            state$wire[[length(state$wire) + 1L]] <- record
            sid <- httr2::resp_header(got, "Mcp-Session-Id")
            if (!is.null(sid)) state$sid <- sid
            if (record$status >= 400L) stop("MCP HTTP ", record$status, ": ", rawToChar(record$body))
            text <- rawToChar(record$body)
            if (!nzchar(text)) return(NULL)
            if (grepl("text/event-stream", httr2::resp_header(got, "content-type"), fixed = TRUE))
                return(get("mcp_sse_payload", asNamespace("aisdk.mcp"))(text))
            return(jsonlite::fromJSON(text, simplifyVector = FALSE))
        }
    }
    adapt <- new.env(parent = asNamespace("gpsadb"))
    adapt$gpsadb_client <- function() {
        if (is.null(state$client)) state$client <- aisdk.mcp::create_mcp_http_client(
            url, timeout = 25, transport = transport)
        state$client
    }
    adapt$gpsadb_call <- function(client, tool, arguments = list()) {
        out <- client$call_tool_data(tool, arguments)
        state$payload <- out
        attr(out, "provenance") <- list(tool = tool, arguments = arguments,
            queried_at = format(Sys.time(), tz = "UTC"), url = url,
            gpsadb_version = as.character(utils::packageVersion("gpsadb")))
        out
    }
    for (n in c("gpsadb_search_datasets", "gpsadb_test_gene_set")) {
        f <- get(n, asNamespace("gpsadb")); environment(f) <- adapt; assign(n, f, adapt)
    }
    function(kind, args) {
        state$wire <- list(); state$payload <- NULL
        f <- adapt[[if (kind == "lookup") "gpsadb_search_datasets" else "gpsadb_test_gene_set"]]
        val <- tryCatch(list(ok = TRUE, value = do.call(f, args)), error = function(e) {
            state$client <- NULL; state$sid <- NULL
            list(ok = FALSE, error = conditionMessage(e), error_class = class(e))
        })
        c(val, list(payload = state$payload, wire = state$wire))
    }
}
paired_requester <- function(root, live = NULL, max_new = 350L, seconds = 1800) {
    state <- new.env(parent = emptyenv()); state$trace <- list(); state$guard <- character()
    state$n_new <- 0L; state$started <- Sys.time()
    state$request <- function(kind, args) {
        key <- paired_hash(list(kind = kind, args = args))
        path <- file.path(root, "raw", paste0(key, ".rds"))
        hit <- file.exists(path)
        if (!hit) {
            if (is.null(live)) {
                state$guard <- c(state$guard, paste("CACHE MISS", kind, key))
                stop(tail(state$guard, 1L))
            }
            if (state$n_new >= max_new || as.numeric(difftime(Sys.time(), state$started, units = "secs")) > seconds) {
                state$guard <- c(state$guard, "Collection budget exhausted; resume collect; no request attempted")
                stop(tail(state$guard, 1L))
            }
            state$n_new <- state$n_new + 1L
            rec <- c(list(key = key, kind = kind, args = args,
                queried_utc = format(Sys.time(), tz = "UTC")), live(kind, args))
            paired_once(rec, path)
        }
        rec <- tryCatch({
            x <- readRDS(path)
            stopifnot(is.list(x), !is.data.frame(x), !anyDuplicated(names(x)),
                all(c("key", "kind", "args", "ok") %in% names(x)),
                identical(x$key, key), identical(x$args, args), identical(x$kind, kind),
                is.logical(x$ok), length(x$ok) == 1L, !is.na(x$ok))
            if (x$ok) stopifnot("value" %in% names(x)) else
                stopifnot(is.character(x$error), length(x$error) == 1L, !is.na(x$error), nzchar(x$error))
            x
        }, error = function(e) {
            message <- paste("CACHE INVALID", kind, key, conditionMessage(e))
            state$guard <- c(state$guard, message)
            stop(message)
        })
        event <- data.frame(kind = kind, key = key, cache_hit = hit, ok = rec$ok,
            error = if (isTRUE(rec$ok)) NA_character_ else rec$error)
        state$trace[[length(state$trace) + 1L]] <- event
        if (!is.null(live)) {
            pathlog <- file.path(root, "collection-events.tsv")
            write.table(cbind(utc = format(Sys.time(), tz = "UTC"), context = state$context, event),
                pathlog, sep = "\t", quote = TRUE, row.names = FALSE, append = file.exists(pathlog),
                col.names = !file.exists(pathlog))
        }
        if (!isTRUE(rec$ok)) stop(rec$error)
        rec$value
    }
    state
}
paired_bind <- function(xs) if (length(xs)) do.call(rbind, xs) else data.frame()
paired_run <- function(root, request_state) {
    profiles <- readRDS(file.path(root, "inputs", "profiles.rds"))
    meta <- read.delim(file.path(root, "inputs", "metadata.tsv"), stringsAsFactors = FALSE)
    results <- list(); traces <- list(); rows <- list(); ledgers <- list(); coverage <- list(); lookups <- list()
    for (k in seq_len(nrow(paired_conditions))) {
        cond <- paired_conditions[k, ]; caps <- as.integer(strsplit(cond$cap, "-", fixed = TRUE)[[1]])
        e <- paired_source(root, cond$code, request_state$request)
        for (i in seq_along(profiles)) {
            z <- profiles[[i]]; warnings <- character(); start <- length(request_state$trace)
            request_state$context <- paste(cond$code, cond$cap, i, z$label, sep = ":")
            pe <- tryCatch(withCallingHandlers(e$pathway_evidence(
                if (is.null(z$ego)) data.frame() else z$ego,
                drivers = meta$expected_label[i], signed = z$signed,
                max_pathways = caps[1], max_datasets_per_driver = caps[2]),
                warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
                error = function(err) err)
            # Both implementations catch service errors internally; infrastructure guards must escape.
            if (length(request_state$guard)) stop(paste(unique(request_state$guard), collapse = "; "))
            id <- paste(cond$code, cond$cap, i, sep = ":")
            results[[id]] <- list(profile = z$label, code = cond$code, cap = cond$cap,
                result = pe, warnings = warnings)
            add <- function(d) cbind(profile = z$label, expected_driver = meta$expected_label[i],
                                     code = cond$code, cap = cond$cap, d, stringsAsFactors = FALSE)
            if (length(request_state$trace) > start) {
                trace <- paired_bind(request_state$trace[seq.int(start + 1L, length(request_state$trace))])
                traces[[id]] <- add(trace)
            } else trace <- data.frame(kind = character(), ok = logical())
            err <- inherits(pe, "error")
            ledger <- if (err) NULL else attr(pe, "perturbation_tests", exact = TRUE)
            lookup <- if (err) NULL else attr(pe, "dataset_lookup", exact = TRUE)
            if (!is.null(ledger) && nrow(ledger)) ledgers[[id]] <- add(ledger)
            if (!is.null(lookup) && nrow(lookup)) lookups[[id]] <- add(lookup)
            if (!err && nrow(pe)) rows[[id]] <- add(as.data.frame(pe))
            counts <- function(status) if (is.null(ledger)) NA_integer_ else sum(ledger$status %in% status)
            coverage[[id]] <- add(data.frame(status = if (err) "implementation-error" else
                if (!nrow(pe)) "no-enrichment" else "graded", n_pathways = if (err) NA_integer_ else nrow(pe),
                candidate_found = if (is.null(lookup)) NA else any(lookup$status == "available"),
                ledger_available = !is.null(ledger), lookup_available = !is.null(lookup),
                n_ledger = if (is.null(ledger)) NA_integer_ else nrow(ledger),
                n_finite = if (is.null(ledger)) NA_integer_ else sum(is.finite(ledger$es) & is.finite(ledger$padj)),
                n_support = counts("supports"), n_conflict = counts("conflicts"),
                n_NS = counts("not-significant"), n_NC = counts("not-computable"), n_failed = counts("failed"),
                observed_lookup_errors = sum(trace$kind == "lookup" & !trace$ok),
                observed_test_errors = sum(trace$kind == "test" & !trace$ok),
                error = if (err) conditionMessage(pe) else NA_character_,
                warnings = paste(warnings, collapse = " | ")))
        }
    }
    list(results = results, rows = paired_bind(rows), coverage = paired_bind(coverage),
         traces = paired_bind(traces), ledger = paired_bind(ledgers), lookup = paired_bind(lookups))
}
paired_transitions <- function(rows) {
    if (!nrow(rows)) return(data.frame())
    out <- list()
    contrasts <- list(c("code", "8-4", "baseline", "fixed"), c("code", "5-2", "baseline", "fixed"),
                      c("cap", "baseline", "8-4", "5-2"), c("cap", "fixed", "8-4", "5-2"))
    for (p in contrasts) {
        variable <- p[1]; held <- if (variable == "code") "cap" else "code"
        d <- rows[rows[[held]] == p[2], ]
        fields <- c("profile", "expected_driver", "pathway", "level", "verdict", "reason")
        a <- d[d[[variable]] == p[3], fields]; b <- d[d[[variable]] == p[4], fields]
        z <- merge(a, b, by = fields[1:3], all = TRUE, suffixes = c("_from", "_to"))
        out[[length(out) + 1L]] <- cbind(contrast = variable, held_value = p[2], from = p[3], to = p[4], z)
    }
    paired_bind(out)
}
paired_main <- function(args = commandArgs(trailingOnly = TRUE)) {
    mode <- if (length(args)) args[1] else "help"
    root <- if (length(args) >= 2) args[2] else "../adjudicate-interpret/results/mechanism-paired-4215"
    root <- sub("/$", "", root)
    if (mode == "help") {
        cat("Rscript dev_docs/mechanism-paired-4215.R freeze|collect|finalize|replay [root] [max-new=350] [seconds=1800]\n",
            "freeze is local only; collect explicitly performs LIVE sequential requests; resumable, failures never retried.\n",
            "finalize verifies complete offline union and seals raw cache; replay requires seal and verifies all checksums.\n")
        return(invisible(NULL))
    }
    if (mode == "freeze") return(paired_freeze(root))
    stopifnot(mode %in% c("collect", "finalize", "replay"))
    paired_check(root, "frozen-manifest.tsv")
    design <- readRDS(file.path(root, "inputs", "design.rds"))
    if (!identical(design$R_version, R.version.string)) stop("R version differs from frozen environment")
    versions <- readRDS(file.path(root, "inputs", "package-versions.rds"))
    now <- vapply(names(versions), function(p) as.character(utils::packageVersion(p)), "")
    if (!identical(versions, now)) stop("External package versions differ from frozen environment")
    sealed <- file.path(root, "raw-manifest.tsv")
    if (mode == "collect" && file.exists(sealed)) stop("Raw cache is sealed; collection forbidden")
    if (mode == "replay" && !file.exists(sealed)) stop("Finalize raw cache before replay")
    if (file.exists(sealed)) {
        m <- paired_check(root, "raw-manifest.tsv")
        raw_files <- sort(list.files(file.path(root, "raw")))
        if (!identical(raw_files, sort(basename(m$file[startsWith(m$file, "raw/")]))))
            stop("Unsealed extra raw cache files found")
    }
    state <- paired_requester(root, live = if (mode == "collect") paired_live() else NULL,
        max_new = if (length(args) >= 3) as.integer(args[3]) else 350L,
        seconds = if (length(args) >= 4) as.numeric(args[4]) else 1800)
    result <- paired_run(root, state)
    if (mode == "finalize") {
        if (file.exists(sealed)) stop("Raw cache is already sealed")
        paths <- sort(list.files(file.path(root, "raw"), full.names = TRUE))
        n_requests <- length(paths)
        if (file.exists(file.path(root, "collection-events.tsv")))
            paths <- c(paths, file.path(root, "collection-events.tsv"))
        paired_tsv(paired_manifest(root, paths), sealed)
        message("Finalized immutable raw cache: ", n_requests, " requests; retained service failures included")
    }
    if (mode == "collect") {
        message("Collection union traversal complete; new requests: ", state$n_new,
                "; run finalize then replay (no network).")
        return(invisible(result))
    }
    if (mode == "replay") {
        out <- file.path(root, "replay")
        if (dir.exists(out)) stop("Replay output exists; refusing overwrite")
        dir.create(out)
        paired_once(result$results, file.path(out, "full-results.rds"))
        for (n in c("rows", "coverage", "traces", "ledger", "lookup"))
            paired_tsv(result[[n]], file.path(out, paste0(n, ".tsv")))
        paired_tsv(paired_transitions(result$rows), file.path(out, "transitions.tsv"))
        paths <- sort(list.files(out, full.names = TRUE))
        paired_tsv(paired_manifest(root, paths), file.path(out, "manifest.tsv"))
        message("Offline replay complete: ", out)
    }
    invisible(result)
}
if (sys.nframe() == 0L) paired_main()
