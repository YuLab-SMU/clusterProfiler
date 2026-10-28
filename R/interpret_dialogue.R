#' Discuss an evidence bundle with a model
#'
#' Open a conversation about an [evidence()] result. Unlike [interpret()], which
#' produces one finished narrative, this returns a *dialogue object* the caller
#' keeps and sends further questions to, so a mechanism hypothesis can be argued
#' with rather than announced.
#'
#' The conversation is bounded by what the evidence says, and the boundary is
#' enforced rather than requested. Three rules hold for every turn:
#'
#' * Only entities that [adjudicate()] already reports may be discussed. A
#'   question about a regulator the evidence does not contain is answered with
#'   that fact, not with model knowledge, because a plausible answer about an
#'   unexamined factor is indistinguishable from an invented one.
#' * Anything the caller adds -- a lab observation, an unpublished result -- is
#'   recorded as a **user assertion** with level `U`. It enters the conversation
#'   and the transcript, and it never changes a verdict: `CEBPB` stays
#'   `single_source` even after the user says their knockdown had no phenotype.
#'   The alternative would let a sentence in the chat raise the grade of a
#'   finding, which is exactly the claim the grades exist to protect.
#' * Every reply is checked against the admit list, so a model that names an
#'   unexamined regulator is caught the same way [interpret()] catches it.
#'
#' @param x An `evidenceBundle` from [evidence()].
#' @param context Optional free text describing the experiment.
#' @param model A model specification understood by `aisdk`, e.g.
#'   `"deepseek:deepseek-v4-flash"`.
#' @param assertions Character vector of things the caller knows that the data
#'   does not, e.g. `"our CEBPB knockdown showed no phenotype"`. Recorded as
#'   level `U` and never used to grade anything.
#' @param max_turns Stop accepting questions after this many. Counted in
#'   questions, not transcript rows: a question writes two rows and an assertion
#'   one, so counting rows made the limit mean different things depending on how
#'   the caller used the dialogue.
#' @param species Species whose annotation bounds the check, e.g.
#'   `"Homo sapiens"`. Taken from the bundle when `NULL`. Without a matching
#'   OrgDb the dialogue refuses to open rather than running a check that would
#'   pass everything.
#' @return An `evidenceDialogue`: a list with `session` (the live
#'   `aisdk::ChatSession`), `transcript` (a data frame of turns), `assertions`,
#'   and `admitted` (the entities a reply may name). Print it to see the
#'   conversation; call [ask()] to continue it.
#'
#'   The session is a reference object and is mutated in place, so assigning the
#'   result of [ask()] rebinds the variable but does not copy the conversation.
#'   Use [fork_dialogue()] to take a branch that can be abandoned.
#' @seealso [ask()], [evidence()], [adjudicate()]
#' @examples
#' \dontrun{
#' ev <- evidence(genes, upstream = c("tf", "perturbagen"),
#'                design = list(treatment = "knockdown"))
#' d  <- interpret_dialogue(ev, context = "airway smooth muscle, dex vs untreated")
#' d  <- ask(d, "Why is CEBPB ranked first?")
#' d  <- say(d, "Our CEBPB knockdown showed no phenotype in this line.")
#' }
#' @export
interpret_dialogue <- function(x,
                               context = NULL,
                               model = NULL,
                               assertions = character(0),
                               max_turns = 50L,
                               species = NULL) {
    if (!inherits(x, "evidenceBundle")) {
        rlang::abort(c(
            "`x` must be an evidenceBundle from evidence().",
            i = "Run evidence() first; a dialogue needs the verdicts to be bounded by."
        ))
    }

    ## The bundle already carries the verdicts adjudicate() computed for it; the
    ## dialogue reads them rather than recomputing, so what is discussed is
    ## exactly what print() and $writable report.
    cv <- x$writable %||% adjudicate(x)$verdicts
    if (is.null(cv) || !nrow(cv)) {
        rlang::abort(c(
            "There is nothing to discuss: this bundle carries no entity the sources graded.",
            i = "Check `x$diagnostics` -- a failed or empty run leaves nothing to talk about."
        ))
    }

    admitted <- cv[cv$verdict %in% c("supported", "single_source"), , drop = FALSE]
    if (!nrow(admitted)) {
        rlang::abort(c(
            "No entity in this bundle reached a supported or single-source verdict.",
            i = "An absence is a result, but it is not something to build a mechanism on."
        ))
    }

    ## The species decides which annotation the boundary check consults; taking
    ## it from the bundle keeps a mouse analysis from being checked against
    ## human symbols, which would pass silently.
    if (is.null(species)) species <- x$input$design$species %||% "Homo sapiens"
    if (!length(mechanism_known_names(species))) {
        rlang::abort(c(
            paste0("No annotation is available for species '", species, "'."),
            i = "Install the matching OrgDb, or a dialogue cannot bound its answers."
        ))
    }

    model <- .normalize_interpret_model(model)

    d <- structure(
        list(
            session     = aisdk::create_chat_session(
                model = model,
                system_prompt = .dialogue_system_prompt(x, admitted, context)
            ),
            transcript  = .dialogue_empty_transcript(),
            assertions  = .dialogue_assertions(assertions),
            admitted    = admitted,
            findings    = consensus_findings(cv),
            species     = species,
            context     = context,
            max_turns   = as.integer(max_turns)
        ),
        class = c("evidenceDialogue", "list")
    )
    d
}

#' Send a question to an evidence dialogue
#'
#' Continues the conversation opened by [interpret_dialogue()]. The reply is
#' checked against the admitted entities before it is returned, so a model that
#' reaches for a regulator the evidence does not contain is reported rather than
#' believed.
#'
#' @param dialogue An `evidenceDialogue` from [interpret_dialogue()].
#' @param question The question to ask.
#' @param ... Ignored, for method compatibility.
#' @return The updated `evidenceDialogue`, with the exchange appended to
#'   `transcript` and any out-of-bounds naming recorded in the last turn's
#'   `unsupported`.
#' @seealso [interpret_dialogue()], [say()]
#' @export
ask <- function(dialogue, question, ...) {
    if (!inherits(dialogue, "evidenceDialogue")) {
        rlang::abort("`dialogue` must come from interpret_dialogue().")
    }
    ## Count questions, not transcript rows: ask() writes two rows and say()
    ## one, so a row limit of 5 allowed three questions and fewer if the caller
    ## had stated anything.
    asked <- if (.nrow0(dialogue$transcript)) {
        sum(dialogue$transcript$role == "assistant")
    } else 0L
    if (asked >= dialogue$max_turns) {
        rlang::abort(c(
            paste0("This dialogue reached its ", dialogue$max_turns, "-turn limit."),
            i = "Raise `max_turns` in interpret_dialogue() if the conversation should continue."
        ))
    }
    if (!is.character(question) || length(question) != 1L || !nzchar(question)) {
        rlang::abort("`question` must be one non-empty string.")
    }

    ## A plain `d2 <- d` shares the R6 session but not the transcript, so the
    ## model's history and print()'s view drift apart with no warning: the user
    ## reads two turns while the model answers from four. Detect the mismatch
    ## and name the fix rather than producing an audit trail that is wrong.
    .dialogue_check_in_sync(dialogue)

    ## A failed call leaves the conversation unchanged, which is right, but the
    ## question the caller just typed must not vanish. It is recorded in the
    ## transcript *before* the call, and rolled back if the session went through
    ## -- so a failure leaves the question visible and a success does not
    ## duplicate it.
    pending <- .dialogue_add_turn(
        dialogue$transcript,
        data.frame(role = "user", text = question, kind = "question",
                   unsupported = NA_character_, stringsAsFactors = FALSE))

    reply <- tryCatch(
        dialogue$session$send(question),
        error = function(e) {
            ## The recorded question has to travel with the error: R copies on
            ## assignment, so a failure inside this function cannot update the
            ## caller's object. `transcript` is a list *field* rather than an
            ## attribute, because stop() drops attributes but keeps fields.
            cond <- structure(
                class = c("evidenceDialogueError", "error", "condition"),
                list(message = paste0("The model call failed: ", conditionMessage(e),
                                      "\nThe question was kept in `$transcript` on this condition; ",
                                      "the session is unchanged."),
                     call = sys.call(-1),
                     transcript = pending))
            stop(cond)
        }
    )
    text <- .dialogue_reply_text(reply)

    ## The same check interpret() runs, applied to a turn instead of a whole
    ## narrative: a regulator the evidence does not contain is reported, not
    ## silently accepted because it sounded competent.
    out_of_bounds <- .dialogue_out_of_bounds(text, dialogue$admitted,
                                             species = dialogue$species %||% "Homo sapiens")

    dialogue$transcript <- .dialogue_add_turn(
        pending,
        data.frame(role = "assistant", text = text, kind = "answer",
                   unsupported = if (length(out_of_bounds))
                       paste(out_of_bounds, collapse = ", ")
                   else NA_character_,
                   stringsAsFactors = FALSE))
    dialogue
}

#' Tell an evidence dialogue something the data does not contain
#'
#' Records a caller-supplied fact -- an unpublished observation, a lab result, a
#' constraint -- so later answers can take it into account.
#'
#' The fact is stored as a **user assertion** at level `U` and is deliberately
#' inert with respect to grading. It cannot raise a verdict, because a sentence
#' typed into a chat is not a measurement, and letting one move `CEBPB` from
#' `single_source` to `supported` would make the grade meaningless.
#'
#' @param dialogue An `evidenceDialogue`.
#' @param statement One thing the caller knows, in their own words.
#' @param ... Ignored.
#' @return The updated `evidenceDialogue`, with the assertion appended and sent
#'   to the model as context for later turns.
#' @seealso [interpret_dialogue()], [ask()]
#' @export
say <- function(dialogue, statement, ...) {
    if (!inherits(dialogue, "evidenceDialogue")) {
        rlang::abort("`dialogue` must come from interpret_dialogue().")
    }
    if (!is.character(statement) || length(statement) != 1L || !nzchar(statement)) {
        rlang::abort("`statement` must be one non-empty string.")
    }

    ## rbind, not c(): combining a data.frame with c() strips it to a list and
    ## the object stops being a table.
    dialogue$assertions <- rbind(
        .dialogue_assertions(dialogue$assertions$statement %||% character(0)),
        data.frame(statement = statement, level = "U", stringsAsFactors = FALSE)
    )
    dialogue$transcript <- .dialogue_add_turn(
        dialogue$transcript,
        data.frame(role = "user", text = statement, kind = "assertion",
                   unsupported = NA_character_, stringsAsFactors = FALSE))

    ## Told to the model as context, phrased so it is not mistaken for evidence
    ## the analysis produced.
    dialogue$session$append_message(
        "user",
        paste0("[caller-supplied context, unverified, level U -- treat as a stated ",
               "premise for discussion, never as support for a finding] ",
               statement)
    )
    dialogue
}

#' Fork a dialogue to explore a side question
#'
#' Copies a dialogue so the branch can be taken somewhere else without losing
#' the original.
#'
#' This exists because the underlying chat session is a reference object: in
#' `d <- ask(d, "q")` the variable is rebound but the *session* is mutated in
#' place, so two dialogue objects that came from the same
#' [interpret_dialogue()] share one conversation. Without a fork, asking a
#' speculative question and then wanting to go back is not possible -- the
#' speculation is already in the history the next reply sees.
#'
#' @param dialogue An `evidenceDialogue`.
#' @param name Optional label for the branch.
#' @return A new `evidenceDialogue` with its own session, whose history starts
#'   at the point the fork was taken.
#' @seealso [interpret_dialogue()], [ask()]
#' @export
fork_dialogue <- function(dialogue, name = NULL) {
    if (!inherits(dialogue, "evidenceDialogue")) {
        rlang::abort("`dialogue` must come from interpret_dialogue().")
    }
    out <- dialogue
    ## The session is an R6 object, so `out <- dialogue` shares it. Cloning is
    ## what actually gives the branch its own history; a copy that still pointed
    ## at the original session would defeat the purpose.
    ##
    ## Note: aisdk::session_fork_branch() returns a branch *identifier*, not a
    ## session, so it is not the call to use here.
    forked <- tryCatch(dialogue$session$clone(), error = function(e) NULL)
    if (is.null(forked)) {
        rlang::abort(c(
            "This session cannot be cloned.",
            i = "The session implementation does not expose clone().",
            i = "Without a copy the branch would share history with the original."
        ))
    }
    out$session <- forked
    out$branch <- name %||% "fork"
    out
}

#' Chat about an evidence bundle, turn by turn
#'
#' Opens an interactive conversation: the findings are printed, then each line
#' typed is sent to the model until the caller leaves. This is the same dialogue
#' [ask()] drives, with the loop written for you.
#'
#' Two prefixes are commands rather than questions:
#'
#' * `!` records what you type as a caller assertion -- something you know that
#'   the data does not. It is stored at level `U` and never changes a verdict.
#' * `:q` leaves, returning the dialogue so the transcript can be kept.
#'
#' Anything else is a question. The bounds are the same as [ask()]: a regulator
#' the evidence does not contain is reported, not discussed.
#'
#' Outside an interactive session this reads from stdin, so it also works when
#' piping a script of questions into `Rscript`.
#'
#' @param dialogue An `evidenceDialogue` from [interpret_dialogue()].
#' @param prompt Text shown before each input line.
#' @param echo Print the model's reply as it arrives. `FALSE` returns the
#'   dialogue silently, for testing.
#' @return The updated `evidenceDialogue`.
#' @seealso [interpret_dialogue()], [ask()], [say()]
#' @examples
#' \dontrun{
#' ev <- evidence(genes, upstream = c("tf", "perturbagen"),
#'                design = list(treatment = "knockdown"))
#' chat(interpret_dialogue(ev, context = "airway smooth muscle"))
#' }
#' @export
chat <- function(dialogue, prompt = "? ", echo = TRUE) {
    if (!inherits(dialogue, "evidenceDialogue")) {
        rlang::abort("`dialogue` must come from interpret_dialogue().")
    }

    if (echo) {
        cat("Discussing ", nrow(dialogue$findings), " finding(s).\n", sep = "")
        cat("  !your note   record something you know (level U, no effect on verdicts)\n")
        cat("  :q           leave\n\n")
    }

    ## One input source for the whole loop; a fresh connection per line would
    ## drop everything after the first line when reading a pipe.
    state <- .dialogue_input_state()
    if (!is.null(state$con)) on.exit(close(state$con), add = TRUE)

    repeat {
        line <- .dialogue_read(prompt, state)
        if (is.null(line)) break                      # EOF
        line <- trimws(line)
        if (!nzchar(line)) next

        ## Any case, and any of the usual spellings: a mistyped command used to
        ## be sent to the model as a question and cost an API call.
        if (tolower(line) %in% c(":q", ":quit", ":exit", ":q!")) break
        if (startsWith(line, ":") && echo) {
            cat("  (:q leaves, !text records an assertion)\n\n")
            next
        }

        if (startsWith(line, "!")) {
            ## Strip every leading "!", then check for content: a bare "!" used
            ## to call say("") which aborted and swallowed the rest of the input.
            note <- trimws(sub("^!+", "", line))
            if (!nzchar(note)) {
                if (echo) cat("  (! needs some text after it)\n\n")
                next
            }
            dialogue <- say(dialogue, note)
            if (echo) cat("  (recorded as an assertion, level U)\n\n")
            next
        }

        dialogue <- ask(dialogue, line)
        if (echo) {
            turn <- dialogue$transcript[nrow(dialogue$transcript), ]
            cat("\n", .dialogue_wrap(turn$text), "\n", sep = "")
            if (!is.na(turn$unsupported)) {
                cat("  ! named without evidence: ", turn$unsupported, "\n", sep = "")
            }
            cat("\n")
        }
    }

    dialogue
}

## One line of input, or NULL at end of file.
##
## aisdk::console_input() gives line editing and history in an interactive
## session and falls back cleanly when there is none.
##
## The piped case must NOT call readLines(con = "stdin") once per line: that
## opens a fresh connection each time, and the second call finds the rest of the
## stream already consumed, so a scripted conversation stopped after one answer.
## One connection is held for the whole loop and closed with it.
.dialogue_read <- function(prompt, state = NULL) {
    if (is.null(state)) state <- .dialogue_input_state()

    if (state$interactive) {
        out <- tryCatch(aisdk::console_input(prompt), error = function(e) NULL)
        if (!is.null(out)) return(out)
    }
    if (is.null(state$con)) return(NULL)
    out <- tryCatch(readLines(con = state$con, n = 1L, warn = FALSE),
                    error = function(e) character(0))
    if (!length(out)) return(NULL)
    out
}

## The input source for one conversation, opened once.
.dialogue_input_state <- function() {
    if (interactive() && exists("console_input", where = asNamespace("aisdk"))) {
        return(list(interactive = TRUE, con = NULL))
    }
    list(interactive = FALSE,
         con = tryCatch(file("stdin", open = "r"), error = function(e) NULL))
}

#' @rdname interpret_dialogue
#' @param x An `evidenceDialogue`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.evidenceDialogue <- function(x, ...) {
    cat("<evidenceDialogue>\n")
    n <- nrow(x$findings)
    cat("  discussing:   ", n, " finding(s) from ", nrow(x$admitted),
        " admitted entit(ies)\n", sep = "")
    if (!is.null(x$context)) {
        cat("  context:      ", substr(x$context, 1, 70),
            if (nchar(x$context) > 70) "..." else "", "\n", sep = "")
    }
    ## A fork is a separate line of enquiry; without the label two branches
    ## print identically and cannot be told apart in a console log.
    if (!is.null(x$branch)) {
        cat("  branch:       ", x$branch, "\n", sep = "")
    }
    if (.nrow0(x$assertions)) {
        cat("  assertions:   ", .nrow0(x$assertions),
            " (level U, do not affect verdicts)\n", sep = "")
    }
    if (.nrow0(x$transcript)) {
        asked <- sum(x$transcript$kind == "question")
        told <- sum(x$transcript$kind == "assertion")
        cat("  turns:        ", asked, " question(s)", sep = "")
        if (told) cat(", ", told, " assertion(s)", sep = "")
        cat("\n\n")

        for (i in seq_len(nrow(x$transcript))) {
            kind <- x$transcript$kind[i] %||% "question"
            ## Assertions and questions both have role "user"; the prefix is
            ## what tells a reader which one they are looking at.
            who <- if (identical(kind, "assertion")) "!" else
                   if (x$transcript$role[i] == "user") "Q" else "A"

            if (identical(kind, "failed")) {
                cat("  A: (no reply -- the model call failed)\n\n")
                next
            }
            cat("  ", who, ": ", .dialogue_wrap(x$transcript$text[i]), "\n\n", sep = "")
            if (!is.na(x$transcript$unsupported[i])) {
                cat("     ! named without evidence: ", x$transcript$unsupported[i], "\n\n", sep = "")
            }
        }
    } else {
        cat("  turns:        0 (use ask() to start)\n")
    }
    invisible(x)
}

# --- internals --------------------------------------------------------------

## nrow() of NULL is NULL, and `if (NULL)` is an error.
.nrow0 <- function(x) if (is.null(x)) 0L else nrow(x)

.dialogue_empty_transcript <- function() {
    data.frame(role = character(), text = character(), kind = character(),
               unsupported = character(), stringsAsFactors = FALSE)
}

.dialogue_assertions <- function(assertions) {
    if (!length(assertions)) {
        return(data.frame(statement = character(), level = character(),
                          stringsAsFactors = FALSE))
    }
    data.frame(statement = as.character(assertions), level = "U",
               stringsAsFactors = FALSE)
}

## The rules the conversation runs under, stated as rules rather than as advice.
.dialogue_system_prompt <- function(x, admitted, context) {
    lines <- c(
        "You are discussing the mechanism hypotheses supported by one enrichment",
        "analysis. The analysis is deterministic and has already decided which",
        "entities the evidence supports; your job is to explain and interrogate",
        "those decisions with the user, not to add to them.",
        "",
        "Rules, all of them binding:",
        "1. Only discuss entities in the list below. If asked about anything else,",
        "   say that this analysis did not examine it, and stop there. Do not",
        "   answer from background knowledge: an unexamined factor that sounds",
        "   plausible is indistinguishable from an invented one.",
        "2. When the user supplies information, treat it as their stated premise.",
        "   It is context, not evidence, and it never changes a verdict.",
        "3. Say which verdict you are relying on when you make a claim.",
        "4. If the evidence cannot settle a question, say so plainly.",
        "",
        "Entities you may discuss, with the verdict the analysis reached:"
    )

    ## Each entity is given with the numbers behind its verdict. Without them a
    ## reply can only restate the verdict, and the user's first question is
    ## always going to be "how strong is that".
    evidence <- .dialogue_evidence_lines(x, admitted)

    rows <- vapply(seq_len(nrow(admitted)), function(i) {
        ## direction is NA for every collection-sourced finding, and paste0()
        ## would render that as the literal "NA" in the prompt.
        dir <- admitted$direction[i]
        head <- paste0("  - ", admitted$regulator[i],
                       ": verdict=", admitted$verdict[i],
                       ", level=", admitted$max_level[i],
                       ", sources=", admitted$sources[i],
                       if (!is.na(dir)) paste0(", direction=", dir) else "")
        extra <- evidence[[admitted$regulator[i]]]
        if (length(extra)) paste0(head, "\n", paste0("      ", extra, collapse = "\n"))
        else head
    }, character(1))

    lines <- c(lines, rows)

    if (!is.null(context)) {
        lines <- c(lines, "", "The user described the experiment as:", context)
    }
    paste(lines, collapse = "\n")
}

## ChatSession$send() returns an aisdk GenerateResult, not a string. The reply
## is its `text` field; anything else reaching here is flattened defensively so
## the bounds check always has a string to work on.
.dialogue_reply_text <- function(reply) {
    if (is.null(reply)) return("")
    if (is.character(reply)) return(paste(reply, collapse = "\n"))
    if (inherits(reply, "GenerateResult") || is.environment(reply)) {
        txt <- tryCatch(reply$text, error = function(e) NULL)
        return(if (is.null(txt)) "" else paste(txt, collapse = "\n"))
    }
    if (is.list(reply)) {
        for (field in c("text", "content", "message")) {
            if (!is.null(reply[[field]])) return(.dialogue_reply_text(reply[[field]]))
        }
        return(paste(unlist(reply), collapse = "\n"))
    }
    ""
}

## Entity-shaped tokens in a reply that no admitted regulator accounts for.
##
## Three things this has to get right, each learned from a bypass:
##
## * Case. Gene names are written "tp53", "TP53" and "Nfkb1" in ordinary prose,
##   so matching capitals only missed most real mentions.
## * Aliases. The literature says "p65" (RELA), "c-Jun" (JUN) and "NF-kB"
##   (NFKB1) far more often than the official symbol. Checking only SYMBOL let
##   every one of those through.
## * Words that are also symbols. WAS, MET and SET are genes; "was", "met" and
##   "set" are words. A lower-case form is read as the word.
##
## The token set comes from the annotation rather than a hand-written pattern:
## a name is a gene name when the OrganismDbi knows it and the analysis did not
## admit it. Anything else is guessing at what a symbol looks like.
.dialogue_out_of_bounds <- function(text, admitted, species = "Homo sapiens") {
    if (!nzchar(text)) return(character(0))
    vocab <- toupper(admitted$regulator)

    ## Multi-part names ("NF-kB", "c-Jun", "NF-kappaB") are matched before
    ## splitting, because a plain word regex would return "NF" and "kB".
    ## Tokens may contain hyphen and a dot, as both occur in aliases.
    toks <- unique(unlist(regmatches(text,
        gregexpr("\\b[A-Za-z][A-Za-z0-9.-]{0,13}[A-Za-z0-9]\\b|\\b[A-Za-z][A-Za-z0-9]{0,13}\\b",
                text))))
    if (!length(toks)) return(character(0))

    names_known <- mechanism_known_names(species)
    if (!length(names_known)) {
        ## A checker that silently passes everything is worse than no checker,
        ## because the caller believes the boundary was enforced.
        rlang::abort(c(
            paste0("No annotation is available for species '", species, "'."),
            i = "The out-of-bounds check cannot run without it.",
            i = "Install the matching OrgDb, or do not open a dialogue for this species."
        ))
    }

    ## Symbols are matched case-insensitively: prose writes tp53, TP53 and
    ## Tp53. Aliases are matched exactly, because they are mixed-case by
    ## design ("p65", "c-Jun") and folding them would collide with words.
    ## A name is either.
    upper_all <- toupper(names_known)
    is_name <- (toupper(toks) %in% upper_all | toks %in% names_known |
                toupper(toks) %in% toupper(vocab))
    candidate <- toks[is_name]
    if (!length(candidate)) return(character(0))

    ## A lower-case form of a symbol that is also an English word is the word:
    ## "was", "met", "set" are not the genes WAS, MET, SET. The same rule keeps
    ## a sentence-initial "The" from being read as a symbol.
    word_like <- grepl("^[a-z]+$", candidate) & !grepl("-", candidate)
    candidate <- candidate[!word_like]

    ## Sentence-case words whose folded form is a symbol: "Was", "Met".
    sentence_case <- grepl("^[A-Z][a-z]{2,}$", candidate)
    candidate <- candidate[!sentence_case]

    ## Some symbols are also common words in capitals: WAS, SET, MET, MAX, IS.
    ## A reply that says "the sample WAS collected" is not naming a gene, and a
    ## warning on every such sentence trains the reader to ignore warnings.
    ## These few are dropped by name; the list is short because it only needs
    ## the words a technical paragraph actually uses in capitals.
    candidate <- candidate[!toupper(candidate) %in% .dialogue_word_symbols()]
    if (!length(candidate)) return(character(0))

    upper <- toupper(candidate)
    ## A single capital plus digits ('A4', 'Q4') is a list marker, not a symbol.
    ## Test the token as written, not its folded form: folding turns the alias
    ## "p65" into "P65", which the pattern would then discard -- exactly the
    ## kind of name this check exists to catch.
    keep <- !grepl("^[A-Z][0-9]+$", candidate)
    candidate <- candidate[keep]
    upper <- upper[keep]
    if (!length(candidate)) return(character(0))

    unexplained <- candidate[!upper %in% vocab]
    upper <- upper[match(unexplained, candidate)]
    ## An admitted family spelling is acceptable; a different member is not.
    unexplained <- unexplained[!vapply(upper,
        function(tok) any(startsWith(vocab, tok)), TRUE)]
    unique(unexplained)
}

## Wrap a paragraph for the console without splitting words.
.dialogue_wrap <- function(text, width = 76) {
    paste(strwrap(text, width = width), collapse = "\n    ")
}

## The numbers behind each admitted entity, grouped by source.
##
## A verdict with no figures behind it invites a conversation that can only
## repeat the verdict, so the per-set detail travels with the prompt.
.dialogue_evidence_lines <- function(x, admitted) {
    r <- x$regulators
    if (is.null(r) || !nrow(r)) return(list())

    keep <- r$regulator %in% admitted$regulator &
        !is.na(r$padj_within_source) & r$padj_within_source < 0.05
    r <- r[keep, , drop = FALSE]
    if (!nrow(r)) return(list())

    by_reg <- split(r, r$regulator, drop = TRUE)
    lapply(by_reg, function(one) {
        ## Take the best row from each source before filling by significance.
        ## Plain head(4) could show four target-set rows while the prompt still
        ## announced "sources=perturbagen,tf", so the model saw a claim of
        ## measured evidence and no measured row to go with it.
        one <- one[order(match(one$source, c("perturbagen", "tf", "mirna")),
                         one$padj_within_source), , drop = FALSE]
        one <- one[!duplicated(one$source), , drop = FALSE]
        rest <- r[r$regulator %in% one$regulator[1] &
                      !r$term %in% one$term, , drop = FALSE]
        if (nrow(rest)) {
            rest <- rest[order(rest$padj_within_source), , drop = FALSE]
            one <- rbind(one, utils::head(rest, 4 - nrow(one)))
        }
        one <- utils::head(one, 4)
        vapply(seq_len(nrow(one)), function(j) {
            ## The dataset's cell line and method decide whether a measured row
            ## is relevant to the caller's system; without them a reply can
            ## only repeat an accession code.
            where <- if (!is.na(one$dataset[j] %||% NA)) {
                paste0(" [", one$dataset[j],
                       if (!is.na(one$system[j] %||% NA)) paste0(", ", one$system[j]) else "",
                       if (!is.na(one$method[j] %||% NA)) paste0(", ", one$method[j]) else "",
                       "]")
            } else ""
            paste0("[", one$source[j], "] ", one$term[j], where,
                   ": ", one$Count[j], "/", one$setSize[j], " genes, padj=",
                   signif(one$padj_within_source[j], 3),
                   if (!is.na(one$direction[j]))
                       paste0(", direction=", one$direction[j]) else "")
        }, character(1))
    })
}

## Symbols that are also ordinary words. Lower-case forms are already handled by
## the word rule; this covers the capitalised and emphasised uses ("the sample
## WAS collected", "SET the value"), which carry no gene meaning.
.dialogue_word_symbols <- function() {
    c("WAS", "SET", "MET", "IS", "ARE", "BE", "IN", "ON", "AT", "TO", "OF",
      "AND", "OR", "NOT", "FOR", "ALL", "ANY", "CAN", "HAS", "HAD", "NOW",
      "ONE", "TWO", "TOP", "MAX", "MIN", "NEW", "OLD", "USE", "USED", "HAS",
      "THE", "THIS", "THAT", "WITH", "FROM", "SUCH", "ONLY", "BOTH", "EACH",
      "IF", "AS", "BY", "DO", "NO", "SO", "UP", "WE", "IT")
}

## Refuse to continue a dialogue whose transcript no longer matches its session.
##
## The session is a reference object and the transcript is a value, so copying
## the dialogue gives two objects that share one conversation. The copy's
## transcript then under-reports what the model has seen, and nothing about the
## output would reveal it.
.dialogue_check_in_sync <- function(dialogue) {
    n <- tryCatch(length(dialogue$session$get_history()), error = function(e) NA_integer_)
    if (is.na(n)) return(invisible(TRUE))

    ## Each ask() writes one user row and one assistant row; each say() writes
    ## one user row. Comparing counts rather than exact contents tolerates a
    ## session that summarises its history.
    expected <- .nrow0(dialogue$transcript)
    if (n > expected) {
        rlang::abort(c(
            "This dialogue's transcript is out of step with its session.",
            x = paste0("The model has seen ", n,
                       " message(s) but the transcript records ", expected, "."),
            i = "This happens when a dialogue is copied with `<-`: the session is shared, the transcript is not.",
            i = "Use fork_dialogue() to branch, or start a new dialogue."
        ))
    }
    invisible(TRUE)
}

## Append rows to a transcript.
##
## rbind() copies the whole table on every turn, which is quadratic in the
## number of turns. The table is small at the default limit, so this keeps the
## simple form and puts the one place that would change if that stops being
## true.
.dialogue_add_turn <- function(transcript, rows) {
    if (is.null(transcript)) transcript <- .dialogue_empty_transcript()
    for (col in setdiff(names(transcript), names(rows))) {
        rows[[col]] <- NA_character_
    }
    rbind(transcript, rows[, names(transcript), drop = FALSE])
}
