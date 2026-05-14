# A swarm-based metabarcoding pipeline

`nf-metabarcoding` is a fast and scalable eDNA workflow built around
[swarm](https://github.com/torognes/swarm). It is a first attempt at
converting the [reference
pipeline](https://github.com/frederic-mahe/fred-metabarcoding-pipeline)
to [Nextflow](https://www.nextflow.io/).


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

- `nextflow >= 23.04.0`
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm` available on `PATH` (or run with a
  container profile)
- `nf-test` to run the test suite (install with
  `curl -fsSL https://code.askimed.com/install/nf-test | bash`)


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
nf-test test tests/processes/merge_fastq_pairs.nf.test

# audit spec <-> test coverage
bash tests/coverage-gate.sh
```

See [`tests/README.md`](tests/README.md) for the TDD workflow and the
pending-test convention.


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
