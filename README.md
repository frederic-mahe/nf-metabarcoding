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
main.nf                <- entry workflow + Part B/C standalone routers
nextflow.config        <- manifest, parameter defaults, profiles
nf-test.config         <- nf-test wiring
modules/local/         <- one .nf per process; functions.nf holds the
                          Groovy helpers (normalize_path, usage, ...)
subworkflows/local/    <- part_a / part_b / part_c subworkflow wiring
SPECIFICATIONS.md      <- behaviour spec (single source of truth)
DECISIONS.md           <- open spec questions blocking [Sxx] IDs
tests/
  README.md            <- how to run the tests
  COVERAGE.md          <- [Sxx] -> test mapping
  coverage-gate.sh     <- audit script (run in CI)
  main.nf.test         <- workflow-level test
  functions/           <- .nf.test for the Groovy helpers
  processes/           <- one .nf.test per process (mirrors modules/local)
  data/
    README.md          <- description of each fixture
    generate.sh        <- reproducible fixture generation
```


## Prerequisites

Required to run the workflow:

- `nextflow >= 25.04.0`
- `bash >= 4`
- `vsearch`, `cutadapt`, `swarm`, `mumu` available on `PATH` — or
  resolved automatically with `-profile conda` (bioconda),
  `-profile modules`, or a container engine
  (`-profile docker`/`podman`/`singularity`/`apptainer`; see "Running
  with containers" below). A container engine only needs the engine
  itself installed — the tools are built into the image.

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
# — no --fastq_pattern needed). --outdir is where every result lands.
nextflow run main.nf \
    --fastq_folder  tests/data \
    --outdir        results \
    --threads       1

# run on your own data — auto-detect via the canonical pattern table
nextflow run main.nf \
    --fastq_folder  /path/to/fastq_dir \
    --outdir        /path/to/results

# multiple folders (comma-separated)
nextflow run main.nf \
    --fastq_folder  /data/run17,/data/run18 \
    --outdir        /path/to/results
```

### Output layout (`--outdir`)

All published artefacts land under `--outdir` (default `results`), in a
fixed layout:

```
<outdir>/
├── per_sample/         per-sample .fas / .qual / .stats + Part A logs
├── occurrence_table/   the occurrence table(s), step logs, and the
│                       taxonomy-annotated tables (Part B + Part C)
└── pipeline_info/      software_versions.yml (tool versions for the run)
```

Inputs are never written to — `--fastq_folder` / `--fasta_folder` are
read-only.

> **Migration note (0.1.0, breaking).** Earlier versions published Part
> A artefacts back into `--fastq_folder` and Part B/C artefacts into
> `--results_folder`. As of 0.1.0 everything goes under `--outdir` in
> the layout above. `--results_folder` still works as a **deprecated
> alias** for `--outdir` (with a warning), but the files now live in the
> `occurrence_table/` sub-directory rather than at the top level — update
> any scripts that read `<results_folder>/*.tsv` to
> `<outdir>/occurrence_table/*.tsv`.

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


## Running with containers

`[S08]`. Four engine profiles run every tool inside a container so the
only thing you install is the engine itself:

```bash
nextflow run main.nf -profile docker      ...   # Docker
nextflow run main.nf -profile podman      ...   # Podman
nextflow run main.nf -profile singularity ...   # Singularity
nextflow run main.nf -profile apptainer   ...   # Apptainer
```

You do **not** build or pull an image by hand. Each profile turns on
[Seqera Wave](https://docs.seqera.io/wave), which builds the image on
the fly from [`environment.yml`](environment.yml) — the same pinned
dependency list the `conda` profile uses — and caches it for reuse.
`environment.yml` stays the single source of truth: bump a version
there and the next run rebuilds the image automatically (see
[`DECISIONS.md`](DECISIONS.md) D10).

Requirements and notes:

- **Outbound network** must be reachable from wherever tasks run (the
  login node *and* the compute nodes), because Wave resolves and builds
  the image at task start.
- Combine with the executor for HPC, e.g. `-profile slurm,singularity`
  (see below). `singularity` / `apptainer` set `autoMounts` so the work
  directory and inputs are visible inside the container.
- Pick **one** dependency mechanism: a container profile, `conda`, or
  `modules` — not several at once.
- Container *execution* is validated by a manual cluster smoke test;
  the profile *wiring* is checked in CI by
  `tests/check-container-profiles.sh`.


## Running on an HPC cluster (slurm)

`[S07]`/`[S08]`. Combine the `slurm` executor profile with **one**
dependency profile (`conda`, `modules`, or a container engine —
`singularity` / `apptainer` are the usual HPC choices, see "Running
with containers" above):

```bash
# tools in a container, built on the fly from environment.yml via Wave
nextflow run main.nf -profile slurm,singularity \
    --fastq_folder /scratch/me/run17 \
    --slurm_queue  normal \
    --threads      8

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
- scales the memory-bound steps off `--dataset_size_gb`: `swarm`
  fastidious global clustering asks for ~3× the dataset size (the
  heaviest step), the other whole-dataset steps (`vsearch` global
  dereplication, chimera detection, the all-vs-all self-search,
  `mumu` and occurrence-table assembly) ~1×. The taxonomic-assignment
  steps (`vsearch --usearch_global`/`--sintax`) are reference-bound
  instead — they load the reference database and a k-mer index, so
  they scale off `--reference_size_gb` (~4× to cover the index).
  Fixed fallbacks apply when those params are unset. On an OOM or
  timeout kill, a process retries up to twice with proportionally
  more memory and wall-time.
- keeps the `[S49]` stampa scatter at its slurm default
  (`stampa_chunk_size = 1000`); the `local` profile sets `0` to feed
  the whole fasta to a single `vsearch` instead.

Launching the run:

- The `nextflow` driver process is long-lived (it stays up submitting
  and reaping sbatch jobs for the whole run), so start it from a login
  node inside `tmux` or `screen` — not as an sbatch job itself, and
  not on a session that will disconnect:

  ```bash
  tmux new -s metabarcoding
  export NXF_OPTS='-Xms512m -Xmx4g'   # cap the driver JVM heap
  nextflow run main.nf -profile slurm,conda \
      --fastq_folder /scratch/me/run17 \
      --slurm_queue  normal \
      --threads      8 \
      -resume                          # reuse cached tasks after a stop
  ```

- `-resume` lets an interrupted run pick up where it left off, but it
  needs the per-task `work/` directories to still exist. The default
  `cleanup = true` deletes them on success, which defeats `-resume`
  across separate invocations — set `cleanup = false` (e.g. a small
  `-c` override) for long or flaky runs, and clean `work/` by hand
  afterwards.

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

The authoritative, per-specification status lives in
[`tests/COVERAGE.md`](tests/COVERAGE.md) (each `[Sxx]` ID → its
test(s) → `done` / `TODO` / `blocked`). The CI suite
(`nf-test test --tag ci`, plus `bats` and `pytest`) is green and the
coverage gate maps every `[Sxx]` in SPECIFICATIONS to at least one
test. Summary at the time of writing:

### Implemented
- Part A — per-sample processing (merge pairs, trim primers,
  dereplicate, extract quality, local clustering with swarm), plus the
  experimental shadow pipeline for unmergeable pairs (`[S04]`)
- Part B — occurrence-table assembly (global dereplication + swarm
  clustering, cluster cleaving, chimera detection, substring-OTU
  merging, mumu post-clustering curation)
- Part C — taxonomic assignment (the stampa scatter-gather and the
  sintax shadow path), with optional majority-rule assignment (`[S66]`)
- `--input` samplesheet input (`[S70]`) and the unified `--outdir`
  output layout (`[S71]`)
- HPC / slurm profile + `conda` / `modules` dependency profiles
  (`[S07]`, `[S08]`) — see "Running on an HPC cluster (slurm)" below.
  The `conda` profile now resolves every tool — including `mumu` —
  from bioconda
- Container profiles — `docker` / `podman` / `singularity` /
  `apptainer` (`[S08]`), built on the fly from `environment.yml` via
  Seqera Wave; see "Running with containers" below
- Multi-layer test suite: nf-test (process + workflow), bats, pytest,
  with shellcheck / flake8 linting and a spec ↔ test coverage gate

### Planned / not yet done
- Two-table occurrence output (`--split-occurrence-table`, `[S15]`)
- Per-sample marker for unmerged-pair clusters (`[S05]`, blocked on
  [`DECISIONS.md`](DECISIONS.md) D02)
- Multiplexed-input subworkflow (demultiplexing is out of scope,
  `[S16]`)
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
