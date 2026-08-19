#!/usr/bin/env Rscript
# scripts/build_refs.R
#
# Quarto pre-render hook (see _quarto.yml: project > pre-render).
# Runs automatically before every `quarto render` / `quarto preview`.
#
# What it does:
#   1. Scans every .qmd/.Rmd source file in the project for pandoc-style
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
# Update this if tidy.bib ever moves.
master_bib_path <- path.expand("~/tidy.bib")

# Where the generated subset gets written (relative to the project root).
out_bib_path <- "book-refs.bib"

# Folders to skip when scanning for citations (build/output/cache dirs).
skip_dirs <- c("_book", "_bookdown_files", "_main_files", "docs",
               ".quarto", "second_ed_prop", "slides", "0X-chapters",
               "Chap09", "renv")

# --- 1. Find source files ----------------------------------------------

all_files <- list.files(
  path = ".",
  pattern = "\\.(qmd|Rmd)$",
  full.names = TRUE,
  recursive = TRUE
)

skip_pattern <- paste0("(^|/)(", paste(skip_dirs, collapse = "|"), ")/")
source_files <- all_files[!grepl(skip_pattern, all_files)]

if (length(source_files) == 0) {
  warning("build_refs.R: no .qmd/.Rmd source files found; ",
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
  regmatches(txt, m)[[1]]
}

used_keys <- unique(unlist(lapply(source_files, extract_keys)))

# Drop Quarto/bookdown cross-reference labels, which use the same @-syntax
# as citations but aren't bibliography entries (e.g. @fig-mymap, @tbl-1,
# @sec-intro, or the bare "ref" left behind by bookdown's \@ref(...)).
# This only cleans up the diagnostic output below -- it has no effect on
# which real citations get pulled into book-refs.bib, since none of these
# would ever match an entry in tidy.bib anyway.
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
  m <- regmatches(first_line,
                   regexec("^@[A-Za-z]+\\s*\\{\\s*([^,\\s]+)\\s*,", first_line))[[1]]
  if (length(m) >= 2) m[2] else NA_character_
}

# --- 4. Keep only the blocks that are actually cited -----------------------

kept <- character(0)
kept_keys <- character(0)

for (i in seq_along(start_idx)) {
  block <- get_block(i)
  key <- get_key(block[1])
  if (!is.na(key) && key %in% used_keys) {
    kept <- c(kept, block, "")   # blank line between entries
    kept_keys <- c(kept_keys, key)
  }
}

missing_keys <- setdiff(used_keys, kept_keys)
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
