#!/usr/bin/env bash
#
# scripts/publish-chapter.sh <chapter-file>.qmd
#
# Publishes ONE chapter's source from `dev` into `master` without dragging
# in whatever else is in-progress on dev, then re-syncs dev so it doesn't
# drift from what's actually live. Mirrors the manual workflow in
# book_github_workflow.txt, with three changes from that doc:
#   - uses `master` (this repo's actual default branch -- there is no
#     `main` here, check with `git branch -a`)
#   - renders explicitly to `docs/` (the GitHub Pages folder) only here,
#     at publish time, on master -- never as a side effect of a plain
#     `quarto render` on dev (see the render step below for why)
#   - pauses twice for a manual look before anything is committed or
#     pushed, so a bad chapter or a missed dependency is easy to back out
#     of instead of already being on GitHub
#
# What this script is FOR: a normal single-chapter content edit -- you
# changed one chapter's .qmd and nothing else, and want just that chapter
# live on master.
#
# What this script is NOT for: a change that touches shared files
# (_quarto.yml, book-refs.bib, scripts/build_refs.R, CSS, more than one
# chapter at once). The `git checkout dev -- "$CHAPTER"` step below pulls
# in exactly one file, by design -- anything else you changed on dev stays
# on dev and quietly does NOT reach master through this script. For that
# kind of change, use scripts/publish-global.sh instead (see that file,
# or Primer_github_workflow.txt, for the full path).
#
# Usage:
#   scripts/publish-chapter.sh 07-direct_comp_mut.qmd
#
# Run this from the repo root, on the `dev` branch, with a clean working
# tree (commit or stash whatever you're mid-edit on first -- the script
# refuses to run otherwise, see the check below).

set -euo pipefail
# -e:  stop immediately if any command fails, instead of plowing on and
#      possibly committing/pushing a half-finished publish
# -u:  treat use of an unset variable as an error (catches typos in the
#      variable names below before they do something surprising)
# -o pipefail: a failure inside a pipeline (e.g. `cmd1 | cmd2`) counts as
#      a failure of the whole pipeline, not just silently swallowed

DEV_BRANCH="dev"
MAIN_BRANCH="master"   # this repo's default/published branch. If you ever
                        # rename it (or start using a different repo where
                        # it really is `main`), this is the only line to
                        # change -- everything below refers to it by name.

# --- Sanity checks -----------------------------------------------------
# Everything in this section only reads state and exits early on problems;
# nothing below this point runs until all three checks pass.

if [ $# -ne 1 ]; then
  echo "Usage: $0 <chapter-file>.qmd" >&2
  exit 1
fi
CHAPTER="$1"

if [ ! -f "$CHAPTER" ]; then
  # Catches typos in the filename before we've switched branches on you --
  # much easier to fix now than after a `git checkout master`.
  echo "Error: '$CHAPTER' not found in $(pwd)." >&2
  exit 1
fi

current_branch=$(git branch --show-current)
if [ "$current_branch" != "$DEV_BRANCH" ]; then
  # This script always starts from dev and ends back on dev. Starting it
  # from anywhere else (master, a stray feature branch, detached HEAD)
  # means the `git checkout dev -- "$CHAPTER"` step later would pull the
  # chapter from the wrong place.
  echo "Error: you're on '$current_branch', not '$DEV_BRANCH'. Switch to $DEV_BRANCH first." >&2
  exit 1
fi

if [ -n "$(git status --porcelain)" ]; then
  # A dirty working tree here is dangerous for two reasons: `git checkout
  # master` a few lines down would carry uncommitted changes onto master
  # with you, and if any of those changes touch files this script doesn't
  # expect to move, they'd get swept into the publish commit by the later
  # `git add -A`. Commit or stash on dev first, then re-run.
  echo "Error: working tree isn't clean. Commit or stash your changes on $DEV_BRANCH first." >&2
  git status --short
  exit 1
fi

echo "== Publishing $CHAPTER from $DEV_BRANCH into $MAIN_BRANCH =="

# --- Human checkpoint 1: did this chapter drag anything else along? ----
# `git checkout dev -- "$CHAPTER"` below is a narrow, single-file pull --
# it does NOT know or care whether your recent commits on this chapter
# also touched book-refs.bib, _quarto.yml, a shared script, etc. If they
# did, those changes are invisible to this script and won't reach master.
# This is just a nudge to look; it can't actually detect the problem for
# you.
echo "-- Does this chapter's recent work touch anything outside its own .qmd?"
echo "   (shared bib entries, _quarto.yml, scripts/, etc. won't come along with a single-file checkout)"
git log --oneline -5 -- "$CHAPTER"
read -r -p "Press Enter to continue, or Ctrl-C to stop and check manually... "

# --- The actual single-file publish ------------------------------------

# Switch the working tree over to master. Everything currently on disk
# now reflects master's last published state, not dev's.
git checkout "$MAIN_BRANCH"

# Pull ONLY this one file's content from dev into the master working
# tree, staged as a change against master. This is not a merge -- it
# doesn't bring dev's commit history, and it doesn't touch any other
# file. Every other chapter on disk right now is still master's
# currently-published version.
git checkout "$DEV_BRANCH" -- "$CHAPTER"

echo "-- Rendering full book on $MAIN_BRANCH (using $MAIN_BRANCH's current state for every other chapter)..."
# A full `quarto render` here rebuilds every chapter, not just $CHAPTER --
# that's what makes this the book as it will actually look once
# published: this chapter's new content next to every other chapter's
# already-published content, cross-references and all.
#
# --output-dir docs is written explicitly, here, and ONLY here. docs/ is
# the git-tracked folder GitHub Pages actually serves, so a render that
# targets it is a real publish action, not a preview. Contrast this with
# a plain `quarto render` while working on dev (no --output-dir flag),
# which goes to the default, gitignored _book/ and touches nothing that
# git or GitHub Pages cares about. If you ever run `quarto render
# --output-dir docs` while still on dev to "just take a quick look," you
# will dirty the live docs/ folder with dev's in-progress, possibly
# unfinished state -- and the moment you `git checkout master` afterward,
# that mess follows you and shows up in the diff below as if it were
# meant to be published. Preview on dev without the flag; only ever add
# the flag here, on master, as part of this deliberate publish step.
quarto render --output-dir docs

# --- Human checkpoint 2: does the diff actually match what you expect? -
# This is the last look before anything becomes permanent. `git status`
# shows which files changed (the chapter's own .qmd, its rendered
# docs/<chapter>.html, its docs/<chapter>_files/ figures, and -- because
# it's a full render -- possibly small incidental diffs in shared files
# like docs/search.json or cross-reference numbers in other chapters if
# this chapter's edits renumbered a shared figure/equation). `git diff
# --stat` gives the shape of the change (file counts and line counts) so
# a render that touched far more than expected stands out at a glance.
echo "-- Diff to be published:"
git status --short
git diff --stat

read -r -p "Looks right? Press Enter to commit + push, or Ctrl-C to bail out (checkout $DEV_BRANCH manually to unwind)... "
# If you Ctrl-C here: the working tree is still on master with these
# changes staged/unstaged but nothing committed. To back out cleanly:
#   git checkout .          # discard the working-tree changes
#   git checkout $DEV_BRANCH
# master itself is untouched on GitHub either way, since nothing has
# been pushed yet at this point.

# `-A` stages everything that changed in the render, not just $CHAPTER --
# deliberately, since the rendered HTML/figures for this chapter (and any
# incidental shared-file diffs noted above) need to be part of the same
# commit as the source change that produced them.
git add -A
git commit -m "Publish: $CHAPTER"
git push origin "$MAIN_BRANCH"

# --- Re-sync dev so it doesn't drift from what's actually live ---------
# dev's copy of $CHAPTER is currently a plain-text-identical copy of what
# was just published (that's what the single-file checkout above
# guaranteed), so this merge is always a fast-forward/no-conflict merge
# on the .qmd itself. It exists to bring master's newly-committed
# rendered output (docs/*.html, docs/*_files/) back onto dev too, so a
# `git status` on dev doesn't show a phantom docs/ diff next time you
# render there, and so dev and master never quietly diverge over time.
#
# --no-edit: if this ever isn't a clean fast-forward (e.g. dev picked up
# other commits since you started), git needs a merge-commit message and
# would otherwise open $EDITOR (Vim) here. --no-edit accepts the
# auto-generated message instead so this never blocks on an interactive
# prompt.
git checkout "$DEV_BRANCH"
git merge "$MAIN_BRANCH" --no-edit

echo "== Done. $CHAPTER is published on $MAIN_BRANCH and $DEV_BRANCH is back in sync. =="
