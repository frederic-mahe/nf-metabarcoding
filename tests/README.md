# Tests

This workflow uses [nf-test](https://www.nf-test.com/) for both
process-level and workflow-level tests. The full TDD cycle is
documented in [`../CLAUDE.md`](../CLAUDE.md).

## Prerequisites

- `nextflow >= 23.04.0`
- `nf-test` (install with `curl -fsSL https://code.askimed.com/install/nf-test | bash`)
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm` available on `PATH` (or run with a
  container profile)

## Running the tests

```bash
# generate fixtures (one-time, or whenever generate.sh changes)
bash tests/data/generate.sh

# run all non-pending tests (this is what CI runs)
nf-test test --tag '!pending'

# run absolutely everything, including red TDD-phase tests
nf-test test

# run a single test
nf-test test tests/processes/merge_fastq_pairs.nf.test

# audit the spec <-> test mapping
bash tests/coverage-gate.sh
```

## Layout

```
tests/
  README.md         <- this file
  COVERAGE.md       <- [Sxx] -> test mapping
  coverage-gate.sh  <- audit script
  main.nf.test      <- workflow-level smoke test
  processes/        <- one .nf.test per process in main.nf
  data/
    README.md       <- description of each fixture
    generate.sh     <- reproducible fixture generation
```

## Writing a new test

1. Find the relevant `[Sxx]` ID in
   [`../SPECIFICATIONS.md`](../SPECIFICATIONS.md). If the behaviour
   isn't specified yet, open a `Dxx` entry in
   [`../DECISIONS.md`](../DECISIONS.md) and stop.
2. Add (or update) the row in [`COVERAGE.md`](COVERAGE.md). Status
   starts as `red`.
3. Write a failing test under `tests/processes/` or `tests/`.
   Tag it with `// COVERAGE: [Sxx]` (one or more IDs).
4. Implement (or fix) the workflow code until the test passes.
5. Move the COVERAGE row from `red` to `done` in the same commit.

Fixtures should be **as small as possible** — a handful of reads is
usually enough to exercise a behaviour. Extend
`tests/data/generate.sh` instead of committing binary files.

## Pending-test idiom

When a test is written ahead of the implementation (the red phase of
TDD) **or** is blocked on a `Dxx` decision, tag it as `pending`:

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

CI runs `nf-test test --tag '!pending'`, so pending tests do not
break the build. Developers see the red by running plain
`nf-test test` locally.

Once the implementation lands, remove the `tag "pending"` line and
fill in the real assertions in the same commit that flips the
COVERAGE row from `red` to `done`.
