# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

prolfquasaint provides SAINT-based protein interaction analysis and
registers the `saint` modelling facade for the prolfqua ecosystem.

## Changelog

Record every change to this package’s code or user-visible behaviour in
`NEWS.md` at the repository root. Add a bullet under the heading for the
current `DESCRIPTION` version, creating a new version heading when a
change opens one. Log the user-visible effect, not the implementation
detail, and match the existing format. Do not log changes that only
touch agent-instruction or repository meta files (`AGENTS.md`,
`CLAUDE.md`, and similar) unless explicitly asked. `make new-version`
bumps the version, commits, and tags but does not write `NEWS.md`, so
add the entry as part of the change itself.
