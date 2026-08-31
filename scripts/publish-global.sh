#!/usr/bin/env bash
#
# scripts/publish-global.sh ["commit message"]
#
# Publishes dev's current state into master as a single "everything as it
# stands on dev" merge -- the counterpart to publish-chapter.sh, for
# changes that touch shared files (_quarto.yml, book-refs.bib,
# scripts/build_refs.R, CSS) or more than one chapter at once, where
# publish-chapter.sh's single-file trick would silently leave part of the
# change behind.
#
# What this script assumes: you've already made your change on dev,
# rendered a full-book preview (plain `quarto render`, no --output-dir, so
# it went to the gitignored _book/, not the live docs/) and eyeballed it,
# and committed the source change on dev yourself. This script picks up
# from there and does the master-side publish: merge dev into master,
# render the live docs/, show you the diff, and wait for you to confirm
# before committing and pushing.
#
# Why this needs its own script and isn't just a one-line `git merge`:
# a merge between two branches that have each moved forward on their own
# is not a fast-forward, so git needs a commit message for the merge and
# opens $EDITOR (Vim, by default) to ask for one -- even when the merge
# itself is conflict-free. Every merge below passes --no-edit to accept
# git's auto-generated message instead, so that prompt never happens.
#
# Usage:
#   scripts/publish-global.sh
#   scripts/publish-global.sh "short description of the change"
#
# Run this from the repo root, on the `dev` branch, with a clean working
# tree (commit your global change on dev first -- the script refuses to
# run otherwise, same guard as publish-chapter.sh).

set -euo pipefail

DEV_BRANCH="dev"
MAIN_BRANCH="master"   # this repo's default/published branch -- see the
                        # matching note in publish-chapter.sh.

# --- Sanity checks -----------------------------------------------------

current_branch=$(git branch --show-current)
if [ "$current_branch" != "$DEV_BRANCH" ]; then
  echo "Error: you're on '$current_branch', not '$DEV_BRANCH'. Switch to $DEV_BRANCH first." >&2
  exit 1
fi

if [ -n "$(git status --porcelain)" ]; then
  # Same reasoning as publish-chapter.sh: a dirty tree here would ride
  # along onto master with `git checkout master` and get swept into the
  # publish commit by the later `git add -A`.
  echo "Error: working tree isn't clean. Commit your global change on $DEV_BRANCH first." >&2
  git status --short
  exit 1
fi

echo "== Publishing $DEV_BRANCH's current state into $MAIN_BRANCH =="

# --- Human checkpoint 1: does dev actually look right? ------------------
echo "-- Last few commits on $DEV_BRANCH going into this publish:"
git log --oneline -5
read -r -p "Press Enter to continue, or Ctrl-C to stop and check manually... "

# --- Merge + full-book render on master ---------------------------------

git checkout "$MAIN_BRANCH"

# --no-edit: see the note up top -- without it, this stops here waiting
# on an interactive editor for a merge-commit message.
git merge "$DEV_BRANCH" --no-edit

echo "-- Rendering full book on $MAIN_BRANCH into docs/ (the live GitHub Pages folder)..."
# Explicit --output-dir docs, here and only here, same as
# publish-chapter.sh -- this is the deliberate publish action, never a
# side effect of a plain `quarto render` on dev.
quarto render --output-dir docs

# --- Human checkpoint 2: does the diff match what you expect? -----------
echo "-- Diff to be published:"
git status --short
git diff --stat

read -r -p "Looks right? Press Enter to commit + push, or Ctrl-C to bail out (checkout $DEV_BRANCH manually to unwind)... "

if [ $# -ge 1 ]; then
  COMMIT_MSG="$1"
else
  read -r -p "One-line commit message describing this publish: " COMMIT_MSG
fi

git add -A
git commit -m "$COMMIT_MSG"
git push origin "$MAIN_BRANCH"

# --- Re-sync dev so it doesn't drift from what's actually live ---------
git checkout "$DEV_BRANCH"
git merge "$MAIN_BRANCH" --no-edit

echo "== Done. $MAIN_BRANCH is published and $DEV_BRANCH is back in sync. =="
