# Specifications

Each testable bullet has a stable identifier (`[Sxx]`). Tests reference
these IDs via `// COVERAGE: [Sxx]` comments and via the rows of
[`tests/COVERAGE.md`](tests/COVERAGE.md). **Never re-number an
existing `[Sxx]`** — add new IDs at the end. Open questions that
block one or more `[Sxx]` IDs live in [`DECISIONS.md`](DECISIONS.md).


## Tests

- `[S00]` test-driven development: every other `[Sxx]` bullet must be
  covered by at least one test
  - **Pass when:** the coverage gate (see
    [`tests/README.md`](tests/README.md)) reports every `[Sxx]` from
    this file mapped to at least one test in `tests/COVERAGE.md`
- tests are minimalistic input files (as small as possible, just
  enough to test a specification)
- the pipeline relies on `vsearch`, `cutadapt`, and `swarm`. These are
  well-tested upstream and **we do not re-test their behaviour** — we
  only test the glue around them (parameter wiring, channel topology,
  output naming, error propagation)


## Workflow structure

- `[S01]` three-part workflow:
    1. **Part A** — fastq files → dereplicated fasta files (merge
       reads, trim primers, extract quality, local-clustering with
       swarm). See `../fred-metabarcoding-pipeline/` for the reference
       implementation
    2. **Part B** — dereplicated fasta files → occurrence table
       (vsearch, swarm, and python scripts). See
       `../fred-metabarcoding-pipeline/` for the reference
       implementation
    3. **Part C** — taxonomic assignment (stampa or sintax): update
       the occurrence table
  - **Pass when:** running the full pipeline on a paired-end fixture
    produces, in order, per-sample `.fas` (Part A), an occurrence
    table (Part B), and a taxonomy-annotated occurrence table (Part C)
- `[S02]` each part can be run separately, or all at once
  - **Pass when:** `--only=A|B|C` (or equivalent flag) runs that part
    in isolation; default runs all three
- `[S03]` fastq files can be paired-end or single-end, compressed
  (`.gz`, `.bz2`) or not
  - **Pass when:** the same input data in any of these forms produces
    byte-identical Part A outputs
- `[S04]` when processing paired-end fastq files, reads that can be
  merged are processed normally; reads that cannot be merged follow a
  parallel pipeline (joined with Ns, Ns converted to As when passed
  to swarm, then converted back), yielding a second occurrence table
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) — N↔A round-trip
    rules need to be defined to avoid rewriting legitimate As
- `[S05]` unmerged-pair clusters appear in the occurrence table with a
  per-sample marker (working name: `sampleID_partial`)
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) — final marker name


## Interfaces between parts

Each part is a directory contract so parts can be tested in
isolation.

- **Part A → Part B**
  - `<sampleId>.fas` — dereplicated FASTA, vsearch
    `--sizeout`-style abundance annotation, `--fasta_width 0`
  - `<sampleId>.qual` — TSV with columns `sha1\tee\tlength`, sorted
    by length then SHA1
  - `<sampleId>.stats` — swarm `--statistics-file` output, filtered
    to clusters with > 2 reads
- **Part B → Part C**
  - one occurrence table (see schema below)
- **Part C output**
  - the occurrence table with taxonomy column updated (no rows
    added or removed)


## Occurrence table schema

- **single-table mode** (default)
  - columns: `cluster_id`, `sequence`, `abundance_total`, one column
    per sample (`<sampleId>` — integer read count), `taxonomy` (Part
    C, empty until then)
  - rows: one per cluster; empty samples (`[S08]`) contribute a
    zero-filled column
- **two-table mode** (`--split-occurrence-table`)
  - `occurrences.tsv` (long format): `clusterId`, `sampleId`,
    `abundance`
  - `clusters.tsv`: `clusterId`, `sequence`, `abundance_total`,
    `taxonomy`


## Workflow requirements

- `[S06]` configurable via a config file or command-line parameters
  - **Pass when:** the same run can be reproduced from either a
    `nextflow.config`/`params.yml` or from `--key value` flags
- `[S07]` runs locally or on HPC with slurm
  - Not part of automated CI; covered by manual smoke tests on the
    cluster. Profiles live in `nextflow.config`.
- `[S08]` runs local or containerized applications (vsearch, swarm,
  cutadapt)
  - Not part of automated CI; covered by a `-profile docker` /
    `-profile singularity` manual smoke test.
- `[S09]` empty input samples must travel through and appear in the
  occurrence table (but not in the two-table mode long-format)
  - **Pass when:** an empty fastq pair produces a row in the
    occurrence table with `abundance = 0` for that sample (and no
    error)
- `[S10]` accepts a directory, or a list of directories (absolute or
  relative paths)
  - **Pass when:** `--fastq_folder a,b,c` or repeated `--fastq_folder`
    arguments are all walked
- `[S11]` searches listed directories automatically using the common
  fastq name patterns listed below
  - **Pass when:** for each pattern in the table, the matching fixture
    is discovered without an explicit `--fastq_pattern`
- `[S12]` automatically deduces sample names from fastq file names
  - **Pass when:** the sample-name derivation strips the matched
    pattern suffix and is idempotent across compression variants
- `[S13]` warns if two or more samples share a derived name
  - **Pass when:** running with two same-named inputs prints a warning
    to stderr and continues (does not abort)
- `[S14]` collision policy for same-named samples
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) — merge / refuse /
    suffix-disambiguate?
- `[S15]` users can choose between a single occurrence table or a
  two-part table (long-format + per-cluster metadata)
  - **Pass when:** `--split-occurrence-table` toggles between the two
    schemas defined above
- `[S16]` expects demultiplexed fastq files
  - **Pass when:** documented in README; demultiplexing is out of
    scope (could be added later as a subworkflow)
- `[S17]` per-cluster minimum-read threshold defaults to > 2 reads
  - **Pass when:** `list_local_clusters` emits no row with
    `reads <= 2` on a fixture that contains a singleton cluster


## Common fastq file-name patterns

The pattern detector recognises (paired-end shown; single-end drops
the `_{1,2}` segment):

| Pattern                                  | Source                          |
|------------------------------------------|---------------------------------|
| `_L001_R{1,2}_001.fastq.gz`              | MiSeq default                   |
| `_L001_R{1,2}.fastq.gz`                  | MiSeq variant                   |
| `_[1-9]_{1,2}.fastq.gz`                  | numeric-lane variant            |
| `_[1-9]_{1,2}_.*.fastq.gz`               | numeric-lane with tail          |
| `_L001_.*_R{1,2}.fastq.bz2`              | bz2 variant of `_L001_R{1,2}…`  |
| `_L005_R{1,2}.fastq.gz`                  | non-default lane number         |
| `_L001_R{1,2}_002.fastq.bz2`             | non-default trailing segment    |
| `_R{1,2}.fastq.gz`                       | minimal R1/R2                   |
| `_{1,2}.fastq.gz`                        | minimal 1/2                     |

- `fastq` may also be `fq`
- compression may be absent: full extension is
  `.(fastq|fq)(.gz|.bz2)?`
- this table is the **single source of truth**; README links here
