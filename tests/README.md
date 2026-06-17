# Tests

Three test runners cover three layers of the workflow:

| Layer                              | Runner    | Lives in           | Speed       |
|------------------------------------|-----------|--------------------|-------------|
| Workflow + process integration     | nf-test   | `tests/*.nf.test`  | seconds     |
| `bin/*.sh` and `bin/*.awk` units   | bats      | `tests/bin/`       | milliseconds|
| `bin/*.py` units                   | pytest    | `tests/python/`    | milliseconds|

All three reference the same `[Sxx]` IDs from
[`../SPECIFICATIONS.md`](../SPECIFICATIONS.md), and the coverage gate
audits all three. The full TDD cycle is documented in
[`../CLAUDE.md`](../CLAUDE.md).

## Prerequisites

- `nextflow >= 23.04.0`
- `nf-test` (install with `curl -fsSL https://code.askimed.com/install/nf-test | bash`)
- `bats` (install with `apt-get install bats` or from
  [bats-core](https://github.com/bats-core/bats-core))
- `python >= 3.10` with `pytest` (and `flake8` for lint)
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm` available on `PATH` (or run with a
  container profile) — only required for nf-test

## Running the tests

```bash
# generate fixtures (one-time, or whenever generate.sh changes)
bash tests/data/generate.sh

# nf-test: only the CI-tagged tests (this is what CI runs)
nf-test test --tag ci

# nf-test: everything, including red TDD-phase tests
nf-test test

# nf-test: a single test
nf-test test tests/processes/part_a/merge_fastq_pairs.nf.test

# bats: unit tests for bin/ shell + awk helpers
bats tests/bin/

# pytest: unit tests for bin/ python helpers
pytest

# audit the spec <-> test mapping (covers all three runners)
bash tests/coverage-gate.sh
```

## Layout

```
tests/
  README.md         <- this file
  COVERAGE.md       <- [Sxx] -> test mapping
  coverage-gate.sh  <- audit script (scans .nf.test, .bats, test_*.py)
  main.nf.test      <- nf-test, workflow-level smoke test
  functions/        <- .nf.test for the Groovy helpers (functions.nf)
  processes/        <- one .nf.test per process, mirroring
    part_a/            modules/local/<part>/<process>.nf
    part_b/
    part_c/
  bin/              <- bats unit tests for bin/*.sh and bin/*.awk
  python/           <- pytest unit tests for bin/*.py (conftest.py
                       adds bin/ to sys.path)
  data/
    README.md       <- description of each fixture
    generate.sh     <- reproducible fixture generation
```

## When to pick which runner

- **nf-test** for anything that needs nextflow's channel topology,
  process wiring, or `vsearch`/`cutadapt`/`swarm`. One nf-test per
  spec `[Sxx]` at minimum; the workflow-level test in
  `tests/main.nf.test` covers integration.
- **bats** for shell and awk helpers in `bin/`. Use it for fast
  feedback on argument handling, edge cases, and stdout shape —
  things you don't want to spin up a nextflow process to check.
- **pytest** for python helpers in `bin/`. Importable modules are
  preferred (faster, more debuggable); CLI-only scripts can be
  exercised via `subprocess.run`. See `tests/python/test_harness.py`
  for both patterns.

A helper covered by bats or pytest should still be reached by at
least one nf-test, so a regression in how nextflow invokes it isn't
silently missed.

## Writing a new test

1. Find the relevant `[Sxx]` ID in
   [`../SPECIFICATIONS.md`](../SPECIFICATIONS.md). If the behaviour
   isn't specified yet, open a `Dxx` entry in
   [`../DECISIONS.md`](../DECISIONS.md) and stop.
2. Add (or update) the row in [`COVERAGE.md`](COVERAGE.md). Status
   starts as `red`.
3. Write a failing test under `tests/processes/<part>/` (Part A, B,
   or C) or `tests/`. Tag it with `// COVERAGE: [Sxx]` (one or more
   IDs).
4. Implement (or fix) the workflow code until the test passes.
5. Move the COVERAGE row from `red` to `done` in the same commit.

Fixtures should be **as small as possible** — a handful of reads is
usually enough to exercise a behaviour. Extend
`tests/data/generate.sh` instead of committing binary files.

## Tagging convention

`nf-test`'s `--tag` filter is **include-only** — there is no built-in
"exclude tag X" mode — so we explicitly opt tests into CI:

- **`tag "ci"`** — the test is expected to pass on every commit. CI
  runs `nf-test test --tag ci`, so only these tests gate the build.
- **`tag "pending"`** — the test is ahead of the implementation
  (red TDD phase) **or** is blocked on a `Dxx` decision. Pending
  tests never run in CI.

Example pending test:

```groovy
test("describes the eventual behaviour") {
    // COVERAGE: [S09]
    tag "pending"
    when { ... }
    then {
        assert workflow.success
        // TODO: assertions describing the intended behaviour
    }
}
```

Once the implementation lands, swap `tag "pending"` for `tag "ci"`
and fill in the real assertions in the same commit that flips the
COVERAGE row from `red` to `done`.

Developers see the full picture (CI tests + reds) by running plain
`nf-test test` locally.
