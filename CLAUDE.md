# `CLAUDE.md` - nf-metabarcoding


## The Golden Rule

When unsure about implementation details, ALWAYS ask the developer.


## Project Context

nf-metabarcoding is an open-source nextflow workflow for microbiome
analysis.

Authoritative documents:
- [`SPECIFICATIONS.md`](SPECIFICATIONS.md) — what the workflow must
  do, with stable `[Sxx]` IDs
- [`DECISIONS.md`](DECISIONS.md) — open questions that block one or
  more `[Sxx]` IDs
- [`tests/README.md`](tests/README.md) — how to run the suite
- [`tests/COVERAGE.md`](tests/COVERAGE.md) — `[Sxx]` → test mapping
- [`README.md`](README.md) — user-facing entry point


### Guidelines

- **always ask a human for instruction** if you think a comment should
  be modified or removed
- in rare occasions the specified behaviour (see
  [`SPECIFICATIONS.md`](SPECIFICATIONS.md)) and the actual behaviour
  of the workflow may differ. If you identify such a situation, it is
  very important that you describe the case and **ask** for a human
  review
- if the README or SPECIFICATIONS lack precision regarding the actual
  behaviour, you can suggest improvements to these files (in chat —
  do not edit `CLAUDE.md` directly without explicit authorization)
- test-driven development (write tests first, each specification must
  be covered by tests)


### Test-driven development

The TDD cycle for this repo:

1. Find or add the `[Sxx]` bullet in
   [`SPECIFICATIONS.md`](SPECIFICATIONS.md). If the spec is unclear,
   open a `Dxx` entry in [`DECISIONS.md`](DECISIONS.md) and stop —
   tests cannot be written against an undefined behaviour.
2. Add or update the row in [`tests/COVERAGE.md`](tests/COVERAGE.md)
   so the `[Sxx]` ID maps to a test file. Status starts as `red`.
3. Write a **failing** test under `tests/processes/<process>.nf.test`
   or `tests/main.nf.test`. Tag it with a comment of the form
   `// COVERAGE: [Sxx]` (one or more IDs, comma-separated). Pending
   tests use `tag "pending"` — see [`tests/README.md`](tests/README.md).
4. Implement (or fix) the workflow code until the test passes.
5. Run `nf-test test` locally. The coverage gate
   (`bash tests/coverage-gate.sh`) verifies that every `[Sxx]` in
   SPECIFICATIONS is mentioned in COVERAGE.md and that every
   `// COVERAGE:` comment references a known ID.
6. Commit test and code together.

Fixtures must stay as small as possible — a handful of reads is
enough to exercise a behaviour. Add new fixtures to
`tests/data/generate.sh` rather than committing binary files.


### Refactoring instructions

In this environment, you have access to `git`, and you have
branch-creation rights:
- create a new git branch, stemming from the `dev` branch, using the
  naming pattern `tmp_$(date +%Y%m%d%H%M%S)` (if your work depends on
  another in-flight `tmp_*` branch that hasn't merged to `dev` yet,
  branch from that one and call it out in chat)
- checkout that new branch
- commit often, it is ok to have small changes in each commit, it
  makes bug-hunting easier
- add "Co-Authored-By" trailers to commit messages, for my co-worker
  Florian Filloux:
  `Co-Authored-By: Florian Filloux <ffillouxdev@users.noreply.github.com>`
- do not add yourself (Claude) as co-author in commit messages (sorry)
- when done, ask for a human review. Humans are in charge of merging
  your work into the dev branch
- check modified bash files with `shellcheck` and fix reported issues
- check modified python files with `flake8` and fix reported issues
- a refactor is "done" when:
  - `nf-test test --tag '!pending'` passes
  - `bash tests/coverage-gate.sh` passes
  - `shellcheck` / `flake8` report nothing on touched files
  - `tests/COVERAGE.md` reflects any spec changes in the same commit


### Bash conventions

- declare one variable per line — do not group declarations such as
  `local a="x" b="y"` or `local x; x="..."` on a single line
- mark variables read-only when possible:
  - top-level constants use `readonly NAME=value` (or `readonly -a`
    for arrays)
  - single-shot function-scope variables use `local -r`
  - inside loop bodies, use plain `local` — `local -r` fails on the
    second iteration because the variable is already read-only
- when assigning the result of a command substitution to a local
  variable, split it onto two lines so the command's exit status is
  not masked by `local`'s (see shellcheck SC2155):

  ```bash
  local q
  q="$(some_command)"
  ```


## What AI Must always Do

AI must always leave these untouched **unless explicitly authorized
by the developer in the current conversation**:
1. **comments** — they're there for a reason
2. **already existing tests** — they've been reviewed by a human
3. **`CLAUDE.md`** — suggest changes in chat; do not edit directly

Remember: We prefer maintainability, correctness, and readability over
cleverness. When in doubt, choose the boring solution or ask for help.
