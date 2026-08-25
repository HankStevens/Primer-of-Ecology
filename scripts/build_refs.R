#!/usr/bin/env Rscript
# scripts/build_refs.R
#
# Quarto pre-render hook (see _quarto.yml: project > pre-render).
# Runs automatically before every `quarto render` / `quarto preview`.
#
# What it does:
#   1. Scans every .qmd source file in the project for pandoc-style
#      citation keys (things like @smith2020, [@jones1999; @lee2001]).
#   2. Reads the full personal reference library (Zotero, kept current on
#      disk via Better BibTeX's auto-export) as PLAIN TEXT.
#   3. Copies out, verbatim, only the raw BibTeX entries whose key is
#      actually cited somewhere in the book, as book-refs.bib, inside the
#      project folder.
#
# Deliberately does NOT use RefManageR (or any BibTeX parser) to read or
# rewrite entries. Re-serializing entries through a parser tends to
# enforce biblatex-strict field rules (e.g. requiring "journaltitle"
# instead of "journal") that a perfectly valid, Zotero-exported, citeproc-
# renderable entry may not satisfy -- which aborts the whole build over
# a field-naming technicality unrelated to whether Quarto can actually
# cite it. Plain-text block copying sidesteps that class of problem
# entirely: whatever citeproc already knows how to render, it keeps
# knowing how to render.
#
# book-refs.bib is what _quarto.yml's `bibliography:` field points at, and
# it IS meant to be committed to git -- it's a small, portable, always-
# current snapshot of just the references this book uses, safe to hand to
# a collaborator or a CI runner without exposing your whole library.
#
# No package dependencies -- base R only.

# --- Configuration ----------------------------------------------------

# Path to your full Zotero-backed library (Better BibTeX auto-export).
# Update this if tidy_zotero.bib ever moves.
master_bib_path <- path.expand("~/tidy_zotero.bib")

# Where the generated subset gets written (relative to the project root).
out_bib_path <- "book-refs.bib"

# Folders to skip when scanning for citations (build/output/cache dirs).
skip_dirs <- c("_book", "_bookdown_files", "_main_files", "docs",
               ".quarto", "second_ed_prop", "slides", "0X-chapters",
               "Chap09", "renv")

# --- 1. Find source files ----------------------------------------------

# .qmd only -- the leftover pre-conversion .Rmd files sitting alongside
# some chapters are no longer part of the book (see _quarto.yml's
# book > chapters list) and only add noise (stale citation keys) here.
all_files <- list.files(
  path = ".",
  pattern = "\\.qmd$",
  full.names = TRUE,
  recursive = TRUE
)

skip_pattern <- paste0("(^|/)(", paste(skip_dirs, collapse = "|"), ")/")
source_files <- all_files[!grepl(skip_pattern, all_files)]

if (length(source_files) == 0) {
  warning("build_refs.R: no .qmd source files found; ",
          "writing an empty book-refs.bib.")
}

# --- 2. Extract citation keys from the source files ----------------------

extract_keys <- function(file) {
  txt <- paste(readLines(file, warn = FALSE), collapse = "\n")

  # Drop fenced R code chunks first, so things like obj@slot or roxygen's
  # @param don't get mistaken for citation keys.
  txt <- gsub("```\\{[^}]*\\}.*?```", "", txt, perl = TRUE)

  # Pandoc citation keys: an @ followed by alphanumerics/._-:etc.
  # (matches @key, [@key], [@key1; @key2], [-@key, p. 12], etc.)
  m <- gregexpr("(?<=@)[a-zA-Z0-9_][a-zA-Z0-9_:.#$%&+?<>~/-]*", txt, perl = TRUE)
  keys <- regmatches(txt, m)[[1]]

  # A narrative citation at the end of a sentence -- "...shown by @tilman1982."
  # with no enclosing brackets -- has no character to stop the match at the
  # real key boundary, so the sentence-ending period gets swallowed into the
  # key (yielding "tilman1982." instead of "tilman1982"), which then can't
  # match anything in tidy_zotero.bib. No real Better BibTeX key ends in
  # trailing punctuation, so it's always safe to trim it off here.
  sub("[.:,;]+$", "", keys)
}

used_keys <- unique(unlist(lapply(source_files, extract_keys)))

# Drop Quarto/bookdown cross-reference labels, which use the same @-syntax
# as citations but aren't bibliography entries (e.g. @fig-mymap, @tbl-1,
# @sec-intro, or the bare "ref" left behind by bookdown's \@ref(...)).
# This only cleans up the diagnostic output below -- it has no effect on
# which real citations get pulled into book-refs.bib, since none of these
# would ever match an entry in tidy_zotero.bib anyway.
crossref_prefixes <- c("fig", "tbl", "eq", "sec", "thm", "lem", "cor",
                        "prp", "exm", "exr", "def", "rem", "sol", "lst",
                        "apx", "nte", "ref")
crossref_pattern <- paste0("^(", paste(crossref_prefixes, collapse = "|"), ")-")
used_keys <- used_keys[!grepl(crossref_pattern, used_keys) & used_keys != "ref"]

# --- 3. Split the master library into raw entry blocks --------------------

if (!file.exists(master_bib_path)) {
  stop(sprintf("build_refs.R: master library not found at %s.",
               master_bib_path), call. = FALSE)
}

lib_lines <- readLines(master_bib_path, warn = FALSE, encoding = "UTF-8")

# Entries always start at the beginning of a line: @type{key, ...
start_idx <- grep("^@[A-Za-z]+\\s*\\{", lib_lines)
n_lines <- length(lib_lines)

get_block <- function(i) {
  from <- start_idx[i]
  to <- if (i < length(start_idx)) start_idx[i + 1] - 1 else n_lines
  block <- lib_lines[from:to]
  # trim trailing blank lines so entries don't accumulate extra spacing
  while (length(block) > 0 && grepl("^\\s*$", block[length(block)])) {
    block <- block[-length(block)]
  }
  block
}

get_key <- function(first_line) {
  # perl = TRUE matters here, not just style: R's default (TRE) regex engine
  # does NOT treat \s as a shorthand class inside a bracket expression --
  # [^,\s] is parsed as "not a comma, not a literal backslash, not a
  # literal letter s", silently failing to parse the key of any entry
  # whose citekey contains the letter "s" anywhere (the vast majority).
  # PCRE mode (perl = TRUE) interprets \s inside brackets correctly.
  m <- regmatches(first_line,
                   regexec("^@[A-Za-z]+\\s*\\{\\s*([^,\\s]+)\\s*,", first_line,
                           perl = TRUE))[[1]]
  if (length(m) >= 2) m[2] else NA_character_
}

# Key of every library entry, in order (parallel to start_idx).
lib_keys <- vapply(seq_along(start_idx), function(i) get_key(get_block(i)[1]),
                    character(1))

# --- 4. Match each used citation key to a library entry --------------------
#
# Better BibTeX auto-generates keys as author+year+shorttitle (e.g.
# "tilman1994bsg"), but citations already written in the book use whatever
# key was current when they were typed -- often an older, shorter form
# ("tilman1994"), or an even older Mendeley-style "author:yearXXsuffix"
# form ("hastings:1980kx"). Rather than requiring every citation in the
# text to be manually updated whenever the library's key-naming changes,
# match more flexibly:
#
#   1. Exact match on the full key.
#   2. The used key is a strict prefix of exactly one library key, with
#      only lowercase letters after it (i.e. it's the same author+year
#      with the Better BibTeX shorttitle suffix chopped off).
#   3. Same idea, but first strip the used key down to its leading
#      "author" + 4-digit "year" (handling a Mendeley-style ":" or "_"
#      separator, or an old single-letter disambiguator like "brown1977a")
#      before trying the prefix match.
#
# If a used key matches more than one library entry this way, it's
# genuinely ambiguous (the author published more than one paper that
# year) and is reported rather than guessed at.
#
# Whichever library entry is chosen, it's written into book-refs.bib
# under the KEY AS CITED IN THE TEXT, not the library's own key -- citeproc
# matches literally against whatever's cited in the .qmd source, so the
# .bib entry has to carry that same key regardless of what Better BibTeX
# calls it internally.

find_candidates <- function(stem) {
  remainder <- substring(lib_keys, nchar(stem) + 1)
  which(!is.na(lib_keys) & startsWith(lib_keys, stem) & grepl("^[a-z]*$", remainder))
}

rewrite_key <- function(first_line, new_key) {
  sub("^(@[A-Za-z]+\\s*\\{\\s*)[^,\\s]+(\\s*,)", paste0("\\1", new_key, "\\2"),
      first_line, perl = TRUE)
}

kept <- character(0)
kept_keys <- character(0)
ambiguous <- character(0)
missing_keys <- character(0)

for (used_key in used_keys) {
  idx <- which(!is.na(lib_keys) & lib_keys == used_key)

  if (length(idx) == 0) {
    idx <- find_candidates(used_key)
  }

  if (length(idx) == 0) {
    m <- regmatches(used_key,
                     regexec("^([a-zA-Z]+)[:_]?([0-9]{4})", used_key))[[1]]
    if (length(m) == 3) {
      stem <- paste0(m[2], m[3])
      if (!identical(stem, used_key)) {
        idx <- find_candidates(stem)
      }
    }
  }

  if (length(idx) == 1) {
    block <- get_block(idx)
    block[1] <- rewrite_key(block[1], used_key)
    kept <- c(kept, block, "")   # blank line between entries
    kept_keys <- c(kept_keys, used_key)
  } else if (length(idx) > 1) {
    ambiguous <- c(ambiguous, sprintf("%s -> {%s}", used_key,
                                       paste(lib_keys[idx], collapse = ", ")))
  } else {
    missing_keys <- c(missing_keys, used_key)
  }
}

if (length(ambiguous) > 0) {
  warning("build_refs.R: ", length(ambiguous),
          " citation key(s) matched MORE THAN ONE library entry, so none ",
          "was included -- pick the intended one and update the .qmd ",
          "source to cite it by its full key: ",
          paste(ambiguous, collapse = "; "))
}
if (length(missing_keys) > 0) {
  warning("build_refs.R: ", length(missing_keys),
          " citation key(s) used in the text were not found in ",
          master_bib_path, ": ", paste(missing_keys, collapse = ", "))
}

# --- 5. Write the project-local subset --------------------------------------

writeLines(kept, out_bib_path)

message(sprintf(
  "build_refs.R: wrote %d reference(s) to %s (scanned %d source file(s)).",
  length(kept_keys), out_bib_path, length(source_files)
))
