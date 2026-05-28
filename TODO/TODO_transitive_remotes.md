# TODO: Verify Transitive Remotes Resolution for prolfquapp Docker

## Context

After the SAINTexpress refactor ([TODO_saint_refactor.md](TODO_saint_refactor.md)),
`prolfquasaint/DESCRIPTION` declares:

```
Imports: saintexpress
Suggests: saintexpressbin
Remotes:
    github::prolfqua/saintexpress,
    github::prolfqua/saintexpressbin
```

`prolfquapp` imports `prolfquasaint` and installs it via
`pak::pkg_install("github::prolfqua/prolfquasaint")` (or `remotes::install_github`)
inside its Docker build. For that to succeed, the installer must walk
`prolfquasaint`'s `Remotes:` field and pull
`prolfqua/saintexpress` (and optionally `prolfqua/saintexpressbin`, since it is
a `Suggests`) from GitHub.

## What to verify

1. **pak** (used in the prolfquapp Dockerfile) honors transitive `Remotes:` when
   installing a GitHub package. Confirm by inspecting the build log or running
   a local dry-run:

   ```bash
   Rscript -e 'pak::pkg_deps("github::prolfqua/prolfquasaint")' \
     | grep -E 'saintexpress|saintexpressbin'
   ```

2. The `prolfquapp` Docker build advances past the `prolfquasaint` install step
   and that `library(prolfquasaint)` works in the resulting image.

3. `runSaint(..., engine = "r")` works without `saintexpressbin` installed (this
   is the soft-fallback path) so the Docker image can omit the binaries if the
   image's licensing posture or size budget excludes them.

## If pak does not resolve transitive Remotes

Options, in order of preference:

- Add explicit `Remotes:` for `saintexpress` and `saintexpressbin` to
  `prolfquapp/DESCRIPTION` so they appear at the top level.
- Pin both packages by SHA or release tag in `prolfquapp`'s install command.
- Pre-install `saintexpress` (and optionally `saintexpressbin`) earlier in the
  Dockerfile before installing `prolfquasaint`.

Do **not** vendor the source of `saintexpress`/`saintexpressbin` into
`prolfquapp` — keep them as standalone GitHub-hosted packages.

## Related

- `prolfquasaint/TODO/TODO_dependencies.md` — the earlier cycle fix
  (`prolfquasaint` must not import `prolfquapp`).
- `prolfquasaint/TODO/TODO_saint_refactor.md` — the SAINTexpress split that
  introduced these new Remotes.
