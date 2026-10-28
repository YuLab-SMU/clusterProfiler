#' Evidence level rules for upstream gene sets
#'
#' The rules [evidence()] applies when it turns an enriched collection into a
#' candidate mechanism. Each row is a rule, evaluated in order, and the first
#' match wins -- so the order is part of the definition, not a display choice.
#'
#' The levels express how strong the claim is, and therefore how it may be
#' written:
#'
#' * **A** -- the collection came from a measured perturbation. Only perturbagen
#'   evidence reaches this level, and only when `design$treatment` reports one.
#' * **B** -- a mechanism is recorded: a ChIP-supported target set, a database
#'   annotation, or a membership in a curated complex.
#' * **C** -- an association only: a predicted target repertoire, a
#'   co-expression pattern, a cell-state or tissue signature. These may be
#'   reported as correlation and never as regulation.
#'
#' @return A `data.frame` with columns `source`, `rule`, `match`,
#'   `default_level` and `note`.
#' @seealso [evidence()]
#' @examples
#' mechanism_xref()
#' @export
mechanism_xref <- function() {
    data.frame(
        source = c("tf", "tf", "tf", "mirna", "perturbagen", "cellstate", "tissue"),
        rule = c(
            "unnamed_motif", "gtrd", "motif",
            "predicted_target", "recorded_perturbation", "signature", "signature"
        ),
        ## `match` is a regular expression against the set identifier; NA matches
        ## everything and must therefore come last within its source block.
        match = c(
            "_UNKNOWN$",
            "_TARGET_GENES$",
            NA_character_,
            NA_character_,
            NA_character_,
            NA_character_,
            NA_character_
        ),
        default_level = c("C", "B", "B", "C", "B", "C", "C"),
        note = c(
            "motif without a named factor: a consensus sequence, not a regulator",
            "ChIP-supported target set (GTRD)",
            "legacy motif set naming a factor; still predicted, not measured",
            "predicted target repertoire (miRDB); no measured binding",
            "transcriptional signature of a known perturbation",
            "expression signature; may reflect cell composition, not a program",
            "tissue signature; same caveat as cellstate"
        ),
        stringsAsFactors = FALSE
    )
}

## Level of one set, given the source it came from and its identifier.
##
## Unnamed motif sets are checked first: 'AACTTT_UNKNOWN' is the single most
## significant set in a typical gene list, and it must not inherit the B grade
## that the generic motif rule would give it.
mechanism_set_level <- function(source, gsid, xref = mechanism_xref()) {
    vapply(as.character(gsid), function(one) mechanism_set_level1(source, one, xref), "")
}

mechanism_set_level1 <- function(source, gsid, xref) {
    ## Same prefix-first dispatch as mechanism_rule_for(), and for the same
    ## reason: the two must agree, or the level and the evidence kind describe
    ## different things about one row.
    lib <- mechanism_union_library(gsid)
    if (!is.na(lib)) {
        ur <- mechanism_union_rules()
        return(ur$level[match(lib, ur$prefix)])
    }
    if (identical(source, "tfunion")) {
        source <- "tf"
        gsid <- sub("^C3TFT:", "", gsid)
    }
    rules <- xref[xref$source == source, , drop = FALSE]
    if (nrow(rules) == 0 || is.na(gsid)) return(NA_character_)
    for (i in seq_len(nrow(rules))) {
        pattern <- rules$match[i]
        if (is.na(pattern) || grepl(pattern, gsid)) {
            return(rules$default_level[i])
        }
    }
    NA_character_
}

## The level actually reported: the rule's grade, demoted to C when the set
## does not name a regulator at all.
##
## The demotion is separate from mechanism_set_level() because the rule table
## matches raw set identifiers, and a name like 'TATA_01' only becomes visible
## as unnamed after mechanism_regulator() has collapsed it. Reporting B for a
## set that names no factor would let a consensus sequence be quoted as a
## mechanism, which is exactly the failure the level column exists to prevent.
mechanism_reported_level <- function(source, gsid, regulator = mechanism_regulator(gsid)) {
    level <- mechanism_set_level(source, gsid)
    ifelse(mechanism_unnamed(regulator), "C", level)
}

## Collapse a set identifier to the regulator it names.
##
## Two naming conventions coexist inside C3:TFT and each needs its own rule:
##   GTRD:   '<TF>_TARGET_GENES'        -> drop the suffix
##   motif:  '[motif_]FFF_Qn', 'FFF_C'  -> drop the motif and the quality tier
##
## A single `sub("_.*$")` is not enough: 'AP1_Q4_01' carries an index after
## the quality tier, and 'TGANTCA_AP1_C' carries a motif before the factor.
## shortcut: motif-only names that never contained a factor ('AACTTT_UNKNOWN')
## still come out as 'UNKNOWN'; they are marked named = FALSE downstream and
## excluded from the supported grade rather than being silently renamed.
## Collapse a set identifier to the entity it names.
##
## The rules differ by source and must not be shared: a TF set identifier is a
## factor name wrapped in a motif prefix and a quality tier, while a
## cell-state or tissue identifier IS the entity name. Applying the TF rules to
## 'CAR_HPX' turns a cardiac signature into the gene HPX, and to
## 'HAY_BONE_MARROW_MONOCYTE' drops the first word. Those sources are therefore
## passed through unchanged.
mechanism_entity <- function(gsid, source, species = "Homo sapiens") {
    if (source %in% c("cellstate", "tissue")) return(as.character(gsid))
    ## The union's added libraries name their factor as the first token of a
    ## free-text set id, which mechanism_regulator()'s motif-tier patterns would
    ## mangle -- 'HIF1A 21447827 ChIP-Seq MCF-7 Human' is not a TRANSFAC name.
    if (identical(source, "tfunion")) {
        return(mechanism_union_entity(gsid, species = species))
    }
    mechanism_regulator(gsid, symbols = mechanism_known_symbols(species))
}

mechanism_regulator <- function(gsid, symbols = NULL) {
    tf <- sub("_TARGET_GENES$", "", gsid)

    ## Strip a quality tier, most specific pattern first. A tier is only
    ## removed when the remainder still ends in a letter, so 'SRF_C' -> 'SRF'
    ## and 'AP1_Q4_01' -> 'AP1', while 'TATA_01' keeps its name instead of
    ## collapsing to '01' (which every downstream test would read as a factor).
    ##
    ## The loop mutates `tf` in place, so each pattern sees the result of the
    ## previous one -- applying all three to the original would let the first
    ## match win even when a later, more specific pattern applies.
    for (pattern in c("_[CQ][0-9A-Z]*_[0-9]+$",   # AP1_Q4_01
                      "_[CQ][0-9]+$",             # HNF3_Q6, AP1_Q4
                      "_[A-Z]$")) {               # SRF_C, WHN_B
        candidate <- sub(pattern, "", tf)
        ## `grepl("[A-Za-z]", ...)` accepted a single character, so 'WHN_B'
        ## collapsed to 'B' -- a quality tier read as a factor name. Two
        ## characters is the shortest real label (AR, ER, GR, PR).
        keep <- candidate != tf & grepl("[A-Za-z]{2}", candidate)
        tf[keep] <- candidate[keep]
    }

    ## A trailing index is not part of a factor name: GATA1_01 and STAT5A_02 are
    ## the first and second motif instances for one factor. Strip it only when
    ## the result is a gene symbol the annotation knows -- that is the check
    ## that separates 'GATA1_01' (an index) from 'NKX2_2' (a real symbol whose
    ## underscore matters) and from 'TATA_01' (a motif name). Without the check
    ## every indexed term became an entity of its own and 101 symbols were
    ## reported under a name GPSAdb cannot look up.
    known <- if (is.null(symbols)) NULL else symbols
    if (!is.null(known)) {
        stripped <- sub("_[0-9]+$", "", tf)
        keep <- stripped != tf & stripped %in% known & !tf %in% known
        tf[keep] <- stripped[keep]
    }

    ## Drop a leading motif. IUPAC covers R Y S W K M B D H V as well as ACGTUN,
    ## so 'TGACCTY_ERR1' and 'GCANCTGNY_MYOD' are handled too; matching only
    ## [ACGTUN] left the motif attached and produced labels naming no factor.
    ##
    ## A motif is only dropped when what follows is itself letter-containing:
    ## 'TATA_01' is a motif name with a position after it, not a motif followed
    ## by a factor, and stripping it would leave '01'.
    label <- sub("^[ACGTRYSWKMBDHVN]+_", "", tf)
    keep <- grepl("[A-Za-z]", label)

    ## An unnamed motif keeps its motif in the label. Collapsing 'AACTTT_UNKNOWN'
    ## and 'TGACATY_UNKNOWN' both to 'UNKNOWN' would report one regulator twice
    ## and make the better-supported of the two impossible to identify.
    unnamed <- !grepl("[A-Za-z]", label) | grepl("^UNKNOWN$", label)
    keep <- keep & !unnamed

    tf[keep] <- label[keep]
    tf
}

## TRUE when the collapsed label is not the name of a factor at all.
mechanism_unnamed <- function(regulator) {
    ## Four ways a label fails to name a factor:
    ##   "UNKNOWN"                  -- the motif had no annotated factor
    ##   "TGACATY_UNKNOWN"          -- same, with the motif still attached
    ##   "CCAWYNNGAAR_UNKNOWN"      -- a degenerate motif that looks like a name
    ##   "01", "TATA_01"            -- a name whose factor never existed
    ## The last one is why a label must contain a letter somewhere: a numeric
    ## remainder is a position or a tier index, never a regulator.
    out <- grepl("^UNKNOWN(_|$)", regulator) |
        grepl("(^|_)UNKNOWN$", regulator) |
        grepl("^[ACGTRYSWKMBDHVN]+_", regulator) |
        !grepl("[A-Za-z]", regulator)
    out[is.na(regulator)] <- FALSE
    out
}

## Every spelling the annotation recognises, not only the official symbol.
##
## Gene symbols appear in prose as protein names and legacy aliases far more
## often than as HGNC symbols: "p65" for RELA, "c-Jun" for JUN, "NF-kB" for
## NFKB1. A checker that knows only SYMBOL silently passes all of them.
##
## Cached alongside the symbol set because the ALIAS table is large and this is
## consulted for every token of every reply.
mechanism_known_names <- function(species = "Homo sapiens") {
    key <- paste0("clusterProfiler_names_", species)
    if (!is.null(.mechanism_cache[[key]])) return(.mechanism_cache[[key]])

    orgdb <- mechanism_orgdb(species)
    syms <- mechanism_known_symbols(species)
    aliases <- if (is.null(orgdb) || !length(syms)) character(0) else {
        a <- tryCatch(
            suppressMessages(AnnotationDbi::select(orgdb, keys = syms,
                                                   columns = "ALIAS",
                                                   keytype = "SYMBOL")),
            error = function(e) NULL)
        if (is.null(a) || !nrow(a)) character(0) else unique(a$ALIAS)
    }
    aliases <- aliases[!is.na(aliases) & nzchar(aliases)]
    ## Case matters here: several aliases differ from an English word only by
    ## case, and folding them would make "was" a gene.
    out <- unique(c(syms, aliases))
    assign(key, out, envir = .mechanism_cache)
    out
}

## TRUE when a label names a protein family or a complex rather than one gene.
##
## 'AP1' is the AP-1 complex (JUN/FOS heterodimers); 'NFAT' is a family of five
## genes; 'MEF2' is four. These are legitimate transcription-factor labels but
## they are not gene symbols, so no perturbation dataset can ever exist for
## them -- GPSAdb indexes gene symbols.
##
## They are reported as such rather than resolved to a member. Mapping AP1 to
## JUN would let a knockdown of JUN stand as evidence about the AP-1 complex,
## which is the same mistake as the NF1 collision: one name, two molecules.
## The organism package's ALIAS table is not a safe substitute either -- it maps
## STAT to SOAT1 (a sterol O-acyltransferase) and OCT1 to POU2F1 as well as to
## SLC22A1 (an organic-cation transporter).
mechanism_is_family_label <- function(regulator, species = "Homo sapiens") {
    if (!length(regulator)) return(logical(0))
    known <- mechanism_known_symbols(species)
    if (!length(known)) return(rep(FALSE, length(regulator)))

    ## TRANSFAC writes some symbols with an underscore where HGNC uses a
    ## hyphen: the collection says NKX2_5, the gene is NKX2-5. Testing only the
    ## literal spelling read those as family names, which they are not -- NKX2-5
    ## is one gene, and it arrives from GTRD, so it is a measurement.
    hyphenated <- gsub("_", "-", regulator, fixed = TRUE)

    ## A name the collections use for one gene is not a family, however little
    ## it looks like a symbol. TRANSFAC names a nuclear receptor by its protein
    ## -- GR, PXR, FXR, HIF2A -- none of which are HGNC symbols, so testing
    ## against the symbol list alone classified each as a family. The cost was
    ## not cosmetic: mechanism_perturb_candidates() skips family labels, so `GR`
    ## was never queried against GPSAdb and could not reach `supported` on a
    ## glucocorticoid experiment, where it is the one answer that is certain.
    ## A name that resolves to a gene is not a family, in either direction of
    ## the alias table. Checking only "which legacy names point here" caught GR
    ## (NR3C1 -> GR) but missed HIF1, which is itself the legacy side of
    ## HIF1 -> HIF1A. The consequence was concrete: on the hypoxia dataset the
    ## HIF1 motif sets reached padj 1.4e-4 and were still dropped from the
    ## perturbation candidates, so HIF1A -- the regulator the experiment is
    ## about -- could not be corroborated by the five HIF1A knockdown datasets
    ## GPSAdb holds.
    aliases <- mechanism_regulator_aliases()
    aliased <- vapply(regulator, function(r) {
        any(names(aliases)[aliases == r] %in% known) ||
            any(aliases[names(aliases) == r] %in% known)
    }, logical(1), USE.NAMES = FALSE)

    !(regulator %in% known | hyphenated %in% known | aliased)
}

## Gene symbols the organism package knows, used to decide whether a stripped
## label names a gene at all. Cached per call because annotation lookups are
## slow and this runs over every term.
mechanism_known_symbols <- function(species) {
    key <- paste0("clusterProfiler_symbols_", species)
    if (!is.null(.mechanism_cache[[key]])) return(.mechanism_cache[[key]])
    orgdb <- mechanism_orgdb(species)
    syms <- tryCatch(AnnotationDbi::keys(orgdb, keytype = "SYMBOL"),
                     error = function(e) character(0))
    ## The environment is locked once the package is loaded, so the cache is
    ## written into it rather than rebound.
    assign(key, syms, envir = .mechanism_cache)
    syms
}

.mechanism_cache <- new.env(parent = emptyenv())

## The OrgDb for a species, or NULL when the package is not installed.
mechanism_orgdb <- function(species) {
    pkg <- switch(species,
        "Homo sapiens" = "org.Hs.eg.db",
        "Mus musculus" = "org.Mm.eg.db",
        "Rattus norvegicus" = "org.Rn.eg.db",
        NULL)
    if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) return(NULL)
    getExportedValue(pkg, pkg)
}

## Labels that name two different things, so agreement between a motif set and
## a perturbation database cannot be read as two sources backing one entity.
##
## Verified case: TRANSFAC's `NF1` is Nuclear Factor I (motif
## NNTTGGCNNNNNNCCNNN, MSigDB TFT_LEGACY), while the HGNC symbol `NF1` is
## neurofibromin, a RAS-GAP tumour suppressor. GPSAdb's dataset D27128 knocks
## out neurofibromin; the C3:TFT terms NF1_Q6 / NF1_Q6_01 / TGCCAAR_NF1_Q6 are
## the Nuclear Factor I motif. Matched on the string, they made `NF1` the
## toolkit's headline verdict (`supported`, level A) for a gene list the two
## sources describe differently.
##
## shortcut: a hand-kept table, not a TRANSFAC-to-HGNC map. Building that map is
## the real fix -- roughly 610 of the 1114 C3:TFT terms are motif-style and any
## of them may use a family name that collides with an unrelated symbol. Extend
## this table when a run turns up another collision; until then a collided label
## is refused the cross-source verdict rather than silently granted it.
mechanism_symbol_collisions <- function() {
    c(NF1 = "TRANSFAC NF1 is Nuclear Factor I; HGNC NF1 is neurofibromin")
}

## TRUE when this label means different things to the collections and to the
## measured sources, so their agreement is a coincidence of spelling.
mechanism_is_collision <- function(regulator) {
    regulator %in% names(mechanism_symbol_collisions())
}

## Measured fact worth knowing before trusting a microRNA result.
##
## A microRNA represses its targets, so its set should be enriched among the
## DOWN-regulated genes and not among the up ones. On the airway list the
## opposite holds: 7 miRNAs were significant in the up half and 0 in the down
## half. miRDB sets are predicted seed matches, so the enrichment is describing
## which genes carry the motif rather than which microRNA is active, and no
## direction is available here to separate the two. That is why "mirna" is not
## in the default `upstream` and why its sentence says so.

## What kind of evidence a row is, as opposed to how strong it is.
##
## `level` answers "how much does this row support a mechanism" and `source`
## answers "which collection produced it". Neither says whether the row is a
## measurement or a prediction, and that distinction is what makes two sources
## agreeing meaningful: a motif set and a knockdown are different kinds of
## evidence, while a motif set and a miRDB repertoire are both predictions and
## their agreement is one observation counted twice.
##
## The mapping is mechanical from the xref rule:
##   recorded_perturbation -> measured   a real perturbation was performed
##   gtrd                  -> measured   ChIP: the factor was found bound
##   motif                 -> predicted  a consensus sequence in a promoter
##   predicted_target      -> predicted  a seed match in a 3' UTR
##   signature             -> observed   expression was measured, but of a
##                                       tissue or state, not of a perturbation
##   unnamed_motif         -> none       names no factor at all
mechanism_evidence_kind <- function(source, rule) {
    kind <- c(
        recorded_perturbation = "measured",
        chip_union            = "measured",
        perturb_signature     = "observed",
        curated_target        = "observed",
        union_motif           = "none",
        gtrd                  = "measured",
        motif                 = "predicted",
        predicted_target      = "predicted",
        signature             = "observed",
        unnamed_motif         = "none"
    )
    out <- unname(kind[rule])
    ## An unrecognised rule is reported as unknown rather than assumed
    ## measured: a wrong "measured" is the kind of error this column exists to
    ## prevent.
    out[is.na(out)] <- "unknown"
    out
}

## Which xref rule a set identifier matches.
##
## Walks the table in the same order as mechanism_set_level1(), so the rule and
## the level can never disagree about which row applied.
mechanism_rule_for <- function(source, gsid, xref = mechanism_xref()) {
    out <- rep(NA_character_, length(gsid))
    if (!length(gsid)) return(out)

    ## A union set's library prefix decides its rule, before any pattern is
    ## consulted: the prefix is a fact about where the set came from, while the
    ## patterns below are guesses from its name.
    lib <- mechanism_union_library(gsid)
    ## A mask, not a subsetting of gsid: the pattern loop below indexes `out` by
    ## positions in the full vector, so shortening gsid would misalign them and
    ## silently assign rules to the wrong sets.
    todo <- is.na(lib)
    if (!all(todo)) {
        ur <- mechanism_union_rules()
        out[!todo] <- ur$rule[match(lib[!todo], ur$prefix)]
    }

    ## A C3TFT:-prefixed set is a C3:TFT set wearing a prefix, so it takes the
    ## `tf` rules on its bare name. Without this the union source has no rows of
    ## its own in the table and every C3:TFT set in it grades as unknown --
    ## which reads as "no evidence", not "graded by a rule we did not apply".
    if (identical(source, "tfunion")) {
        source <- "tf"
        gsid <- sub("^C3TFT:", "", gsid)
    }

    rules <- xref[xref$source == source, , drop = FALSE]
    if (!nrow(rules) || !any(todo)) return(out)

    ## Vectorised over `gsid`; the first matching row wins, as in
    ## mechanism_set_level1(). Assigning to the first still-unresolved hit keeps
    ## the order without a per-element loop.
    gsid <- as.character(gsid)
    for (i in seq_len(nrow(rules))) {
        pattern <- rules$match[i]
        hit <- if (is.na(pattern)) todo else todo & grepl(pattern, gsid)
        out[hit] <- rules$rule[i]
        todo[hit] <- FALSE
        if (!any(todo)) break
    }
    out
}

## The evidence kind of a set identifier, using the same rule table as grading.
mechanism_kind_of <- function(source, gsid, xref = mechanism_xref()) {
    mechanism_evidence_kind(source, mechanism_rule_for(source, gsid, xref))
}

## Bring a regulator label to one spelling, so that evidence from different
## sources about the *same factor* is recognised as the same entity.
##
## The collections and the perturbation databases do not agree on symbol
## spelling: C3:TFT names the hypoxia factor 'HIF1' (the legacy symbol) while
## GPSAdb names it 'HIF1A'. Without this step a factor supported by both its
## target sets and a knockdown is two findings, and consensus() can never
## report the agreement that a mechanism needs.
##
## shortcut: a fixed alias table, not an OrgDb lookup. OrgDb would be more
## complete but needs a species package and a keytype round trip for every
## label; extend the table when a real analysis hits a miss. Aliases are
## applied after mechanism_regulator(), so only collapsed labels are rewritten.
mechanism_regulator_aliases <- function() {
    c(
        ## Legacy symbol -> current symbol. C3:TFT was built over several
        ## releases and mixes both spellings; GPSAdb uses the current one.
        HIF1  = "HIF1A",
        EPAS1 = "HIF2A",
        ## Nuclear-receptor gene symbol -> the name the collections use.
        NR3C1 = "GR",
        NR1I2 = "PXR",
        NR1H4 = "FXR"
    )
}

mechanism_alias_regulator <- function(regulator) {
    aliases <- mechanism_regulator_aliases()
    hit <- match(regulator, names(aliases))
    out <- regulator
    out[!is.na(hit)] <- aliases[hit[!is.na(hit)]]
    unname(out)
}

#' Names a regulator is known by in a perturbation database
#'
#' @description
#' The name a collection uses is not always the name a database uses, and the
#' difference silently costs the only evidence that can reach level A.
#'
#' @details
#' C3:TFT names a nuclear receptor by the motif it was built from -- `GR`,
#' `PXR`, `FXR` -- while GPSAdb indexes the same factor by its gene symbol,
#' `NR3C1`, `NR1I2`, `NR1H4`. `mechanism_alias_regulator()` normalises towards
#' the collection, which is right for the verdict table and wrong for the
#' lookup: querying GPSAdb for `GR` returns nothing, so the regulator can never
#' collect perturbation evidence and can never be reported as `supported`.
#'
#' Measured on the airway dexamethasone list, that is exactly what happened to
#' `GR`: the only regulator a glucocorticoid experiment is guaranteed to
#' involve came out `single_source` at level B, while GPSAdb held two datasets
#' for `NR3C1` that both move the list at padj 1e-22.
#'
#' @param regulator Regulator names as the collections spell them.
#' @return A list, one character vector per input, holding the name itself and
#'   every gene symbol that normalises to it.
#' @name mechanism-gpsadb-names
#' @keywords internal
mechanism_gpsadb_names <- function(regulator) {
    aliases <- mechanism_regulator_aliases()
    lapply(regulator, function(r) {
        ## Both directions: a collection may use either spelling, and asking
        ## under both costs one extra request only when they differ.
        extra <- c(names(aliases)[aliases == r], aliases[names(aliases) == r])
        unique(c(r, unname(extra)))
    })
}

## Grading for the transcription-factor union.
##
## The union's sets carry their library as a prefix, and each library's evidence
## is worth something different, so the prefix decides the rule rather than a
## pattern-matched row in mechanism_xref(). Folding them into that table would
## need six rows and would still get C3TFT: sets wrong, because those keep the
## motif-versus-GTRD distinction the plain `tf` source already makes.
mechanism_union_rules <- function() {
    data.frame(
        prefix = c("ENCODE", "ChEA", "TFP", "TRRUST"),
        rule   = c("chip_union", "chip_union", "perturb_signature", "curated_target"),
        level  = c("B", "B", "B", "C"),
        stringsAsFactors = FALSE
    )
}

## The library prefix of a union set, or NA when it is a C3:TFT set or a plain id.
mechanism_union_library <- function(gsid) {
    pre <- mechanism_union_prefix(gsid)
    rules <- mechanism_union_rules()
    out <- rep(NA_character_, length(gsid))
    hit <- !is.na(pre) & pre %in% rules$prefix
    out[hit] <- pre[hit]
    out
}

## Which way a perturbation moves the factor, in the direction the caller's own
## contrast is stated in.
##
## One place, because both the sign test and the concordance flag have to agree
## about it, and a second copy would drift. GPSAdb's `method` is a mixed
## vocabulary ("siRNA", "knockout") while `method_label` is normalised
## ("Knock_down", "Knock_out"), so the label is read first and the raw method
## second.
##
## "decrease" means the factor was removed, so its targets should FALL when it
## is; "increase" means it was supplied, so its targets should RISE. Anything
## unrecognised returns NA rather than a guess -- a wrong direction here would
## silently invert every concordance call it touches.
mechanism_method_sign <- function(method, method_label = NULL) {
    lab <- tolower(gsub("[^a-z]", "", tolower(paste(method_label %||% ""))))
    raw <- tolower(gsub("[^a-z]", "", tolower(paste(method %||% ""))))
    txt <- paste(lab, raw)
    out <- rep(NA_character_, length(txt))
    ## Removal first: "knockout" contains "knock" and "overexpression" does not.
    remove <- grepl("knockdown|knockout|sirna|shrna|crispr|deletion|inhibition", txt)
    add <- grepl("overexpression|overexpress|activation|induction|oe", txt)
    out[remove] <- "decrease"
    out[add & !remove] <- "increase"
    out
}

## The direction the responding genes should move if the regulator drives the
## programme the collection layer found.
##
## The collection layer says the programme moved `direction` in the caller's
## contrast. A factor that was removed cannot have pushed it that way: removing
## an activator makes its targets fall, so a programme seen rising under a
## knockdown is one the factor was restraining. Flipping the expectation for a
## removal is what makes the two observations comparable at all.
mechanism_expected_measured <- function(direction, method_sign) {
    out <- rep(NA_character_, length(direction))
    ok <- !is.na(direction) & !is.na(method_sign)
    flip <- function(d) ifelse(d == "up", "down", ifelse(d == "down", "up", NA))
    out[ok] <- ifelse(method_sign[ok] == "decrease",
                      flip(direction[ok]), direction[ok])
    out
}
