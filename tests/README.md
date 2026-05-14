# Tests

This workflow uses [nf-test](https://www.nf-test.com/) for both
process-level and workflow-level tests.

## Prerequisites

- `nextflow >= 23.04.0`
- `nf-test` (install with `curl -fsSL https://code.askimed.com/install/nf-test | bash`)
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm` available on `PATH` (or run with a container profile)

## Running the tests

```bash
# generate fixtures (one-time, or whenever generate.sh changes)
bash tests/data/generate.sh

# run all tests
nf-test test

# run a single test
nf-test test tests/processes/merge_fastq_pairs.nf.test
```

## Layout

```
tests/
  README.md         <- this file
  COVERAGE.md       <- spec bullet -> test mapping
  main.nf.test      <- workflow-level smoke test
  processes/        <- one .nf.test per process in main.nf
  data/
    README.md       <- description of each fixture
    generate.sh     <- reproducible fixture generation
```

## Writing a new test

1. Find the relevant specification bullet in `../SPECIFICATIONS.md`.
2. Add a row to `COVERAGE.md` (or update an existing one).
3. Write a failing test under `tests/processes/` or `tests/`.
4. Implement (or fix) the workflow code until the test passes.
5. Commit test + code together.

Fixtures should be **as small as possible** — a handful of reads is
usually enough to exercise a behaviour. Prefer extending
`tests/data/generate.sh` over committing binary files.
