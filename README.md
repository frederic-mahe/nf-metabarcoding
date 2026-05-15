# A swarm-based metabarcoding pipeline

`nf-metabarcoding` is a fast, user-friendly, and scalable eDNA
workflow built around [swarm](https://github.com/torognes/swarm). It
is a first attempt at converting the [reference
pipeline](https://github.com/frederic-mahe/fred-metabarcoding-pipeline)
to [Nextflow](https://www.nextflow.io/). Its modular structure allows
for meta-studies grouping several independent datasets.


## Project layout

```
main.nf                <- workflow definition
nextflow.config        <- manifest, profiles
nf-test.config         <- nf-test wiring
SPECIFICATIONS.md      <- behaviour spec (single source of truth)
DECISIONS.md           <- open spec questions blocking [Sxx] IDs
tests/
  README.md            <- how to run the tests
  COVERAGE.md          <- [Sxx] -> test mapping
  coverage-gate.sh     <- audit script (run in CI)
  main.nf.test         <- workflow-level test
  processes/           <- one .nf.test per process in main.nf
  data/
    README.md          <- description of each fixture
    generate.sh        <- reproducible fixture generation
```


## Prerequisites

Required to run the workflow:

- `nextflow >= 25.04.0`
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm` available on `PATH` (or run with a
  container profile)

Required to run the test suite and the linters:

- `nf-test >= 0.9.0` (install with
  `curl -fsSL https://code.askimed.com/install/nf-test | bash`)
- `bats >= 1.10` (install with `apt-get install bats` or from
  [bats-core](https://github.com/bats-core/bats-core)) — for
  `tests/bin/*.bats`
- `python >= 3.10` with `pytest` — for `tests/python/test_*.py`
- `flake8` — lints `bin/*.py`
- `shellcheck` — lints `bin/*.sh`


## How to run

```bash
# generate test fixtures (one-time)
bash tests/data/generate.sh

# run on the bundled fixtures
nextflow run main.nf \
    --fastq_folder  tests/data \
    --fastq_pattern '/paired_merge_ok_{1,2}.fastq.gz' \
    --threads       1

# run on your own data
nextflow run main.nf \
    --fastq_folder  /path/to/fastq_dir \
    --fastq_pattern '/*_R{1,2}_001.fastq.gz'
```

See [`SPECIFICATIONS.md`](SPECIFICATIONS.md) for the full list of
supported file-name patterns and parameters.


## How to test

```bash
# run the CI-tagged tests (pending tests excluded)
nf-test test --tag ci

# run a single test file
nf-test test tests/processes/part_a/merge_fastq_pairs.nf.test

# audit spec <-> test coverage
bash tests/coverage-gate.sh
```

See [`tests/README.md`](tests/README.md) for the TDD workflow and the
pending-test convention.


## Cleaning up after a run

`cleanup = true` in [`nextflow.config`](nextflow.config) removes
per-task `work/` directories on success, but each runner leaves some
state behind:

```bash
# nextflow runs (in the project root or wherever you launched from)
rm -rf work/ .nextflow.log* .nextflow/

# published outputs land in --fastq_folder; remove them from the
# bundled fixture dir after demo runs
rm -f tests/data/*.fas tests/data/*.qual tests/data/*.stats

# nf-test runs
rm -rf .nf-test/ .nf-test.log

# pytest cache
rm -rf .pytest_cache/ tests/python/__pycache__/
```


## Releasing

The version number lives in two files that must be kept in sync:
[`nextflow.config`](nextflow.config) (`manifest.version`) and
[`CITATION.cff`](CITATION.cff) (`version:`). When bumping, update
both and also refresh `CITATION.cff`'s `date-released`.


## Status

### Implemented
- Part A: per-sample processing — merge pairs, trim primers, derep,
  local clustering with swarm
- Per-process nf-test scaffolding chained via `setup{}` blocks
- Reproducible fixture generation (`tests/data/generate.sh`)

### Partial / in progress
- Workflow-level smoke test (currently `red` on `[S01]`, `[S03]`,
  `[S09]`)
- Auto-detection of fastq naming patterns (`[S11]`, `[S12]`)
- Empty-sample passthrough (`[S09]`)

### Planned
- Part B: occurrence table assembly (`[S15]`)
- Part C: taxonomic assignment
- HPC / slurm profile (`[S07]`)
- Container profiles (`[S08]`)
- Multiplexed-input subworkflow
- Reference database auto-deduction from primers
- Force-rerun controls (see "Re-run policy" below)


## Re-run policy

These inputs should trigger a complete or partial re-run when they
change (currently a developer note — not yet enforced by the
workflow):

- new versions of `vsearch`, `cutadapt`, `swarm`
- new versions of external python scripts
- new version of the reference database
- new set of fastq files
- user-requested global or partial re-run
