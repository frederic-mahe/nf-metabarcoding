# A swarm-based metabarcoding pipeline

`nf-metabarcoding` is a fast, user-friendly, and scalable eDNA
workflow built around [swarm](https://github.com/torognes/swarm). It
is a first attempt at converting the [reference
pipeline](https://github.com/frederic-mahe/fred-metabarcoding-pipeline)
to [Nextflow](https://www.nextflow.io/). Its modular structure allows
for meta-studies grouping several independent datasets.

The pipeline aims to produce single-nucleotide resolution clustering
results, equivalent to Amplicon Sequence Variants (ASVs). Clustering
results, chimera detection results, and taxonomic assignment results
are grouped in a ready-to-use occurrence table.


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
- `vsearch`, `cutadapt`, `swarm`, `mumu` available on `PATH` — or
  resolved automatically with `-profile conda` (bioconda) or
  `-profile modules` (see "Running on an HPC cluster (slurm)" below)

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

# run on the bundled fixtures (paired_merge_ok matches canonical row 7
# — no --fastq_pattern needed)
nextflow run main.nf \
    --fastq_folder  tests/data \
    --threads       1

# run on your own data — auto-detect via the canonical pattern table
nextflow run main.nf \
    --fastq_folder  /path/to/fastq_dir

# multiple folders (comma-separated)
nextflow run main.nf \
    --fastq_folder  /data/run17,/data/run18
```

### Input discovery

The pipeline globs every fastq file in the listed folders (`*.fastq`,
`*.fq`, with optional `.gz` / `.bz2`) and identifies paired-end pairs
using the canonical pattern table documented in
[`SPECIFICATIONS.md`](SPECIFICATIONS.md). Files that match no
paired-end pattern (or whose R2 partner is missing) are processed as
single-end samples that skip the merging step.

### Adding a custom pattern (`--fastq_pattern`)

If your file names don't match any canonical row, supply your own
glob via `--fastq_pattern`. Format: a shell-style glob containing the
literal token `{1,2}` (or `{R1,R2}`) that marks the R1/R2
discriminator. A `*` before the discriminator captures the sample ID.

```bash
# study42-mateA.fastq.gz + study42-mateB.fastq.gz → sample "study42"
nextflow run main.nf \
    --fastq_folder  /data/run17 \
    --fastq_pattern '*-mate{A,B}.fastq.gz'

# fixed prefix, no '*': the literal prefix (minus trailing _/./-)
# becomes the sample ID
nextflow run main.nf \
    --fastq_folder  tests/data \
    --fastq_pattern 'paired_merge_ok_{1,2}.fastq.gz'
```

User patterns take precedence over the canonical table; supported
glob meta-characters are limited to `*` (any chars) and the `{R1,R2}`
brace token — anything else is matched literally.


## Running on an HPC cluster (slurm)

`[S07]`/`[S08]`. Combine the `slurm` executor profile with **one**
dependency profile:

```bash
# tools from bioconda (Nextflow builds the env once and caches it)
nextflow run main.nf -profile slurm,conda \
    --fastq_folder /scratch/me/run17 \
    --slurm_queue  normal \
    --threads      8

# tools from the cluster's module system (set the module names/versions)
nextflow run main.nf -profile slurm,modules \
    --fastq_folder    /scratch/me/run17 \
    --slurm_queue     normal \
    --module_vsearch  vsearch/2.31.0 \
    --module_swarm    swarm/3.1.5 \
    --module_cutadapt  cutadapt/4.9 \
    --module_mumu     mumu/1.1.1 \
    --threads         8
```

What the `slurm` profile does:

- submits every process as an sbatch job (`process.executor =
  'slurm'`); tune the queue/account/concurrency with `--slurm_queue`,
  `--slurm_account`, `--slurm_queue_size`.
- sets per-process resources by tier. `--threads` is the single knob
  for cores: it feeds `task.cpus`, which is what the tools actually
  request (`vsearch/swarm/cutadapt --threads`/`--cores`). The
  single-threaded `mumu` and the bash/awk/python glue stay at one
  core; `chimera_detection` reserves two (it is a
  `filtering | uchime_denovo` pipe).
- scales the two memory-bound global steps off the dataset size:
  `vsearch` global dereplication asks for ~1× and `swarm` fastidious
  global clustering ~3× of `--dataset_size_gb` (fixed 16 GB / 64 GB
  fallbacks when it is not set). On an OOM or timeout kill, a process
  retries up to twice with proportionally more memory and wall-time.
- keeps the `[S49]` stampa scatter at its slurm default
  (`stampa_chunk_size = 1000`); the `local` profile sets `0` to feed
  the whole fasta to a single `vsearch` instead.

Notes:

- Run `work/` on shared scratch visible to every compute node, and
  point `NXF_CONDA_CACHEDIR` at shared storage so the conda env is
  built once and reused.
- `publishDir` defaults to hard links (`--publish_mode link`). Switch
  to `--publish_mode copy` when `--results_folder`/`--fastq_folder`
  live on a different filesystem than `work/`.
- This path is not part of automated CI; it is exercised by a manual
  smoke test on the cluster (`[S07]`/`[S08]` in
  [`tests/COVERAGE.md`](tests/COVERAGE.md)).


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
- HPC / slurm profile + conda / module dependency profiles (`[S07]`,
  `[S08]`) — see "Running on an HPC cluster (slurm)" below

### Partial / in progress
- Workflow-level smoke test (currently `red` on `[S01]`, `[S03]`,
  `[S09]`)
- Empty-sample passthrough (`[S09]`)

### Planned
- Part B: occurrence table assembly (`[S15]`)
- Part C: taxonomic assignment
- Container profiles — docker / singularity (`[S08]`)
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
