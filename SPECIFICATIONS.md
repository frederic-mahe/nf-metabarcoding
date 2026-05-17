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
  relative paths). `--fastq_folder` accepts a single path or a
  comma-separated list (`--fastq_folder a,b,c`); a `nextflow.config`
  may instead supply a Groovy list (`fastq_folder = ['a', 'b']`).
  - **Pass when:** every fastq file in every listed folder is
    discovered; identical inputs spread across two folders produce
    identical artefacts to the single-folder run.
- `[S11]` builds the input list by globbing every fastq file
  (`*.fastq`, `*.fq`, with optional `.gz` / `.bz2`) in the listed
  directories, then identifies paired-end pairs by matching the R1
  basename against the canonical pattern table below. Users may
  prepend a custom pattern via `--fastq_pattern` (a glob containing
  the `{1,2}` brace token that marks the R1/R2 discriminator —
  e.g. `*_run17_{1,2}.fastq.gz`); custom patterns take precedence
  over the canonical list.
  - **Pass when:** for each row of the canonical table, an R1 fixture
    is paired with its R2; an R1 fixture matching `--fastq_pattern`
    but not any canonical pattern is also paired correctly.
- `[S12]` derives the sample ID for a paired-end pair from the R1
  basename by stripping the **first matching pattern** suffix (the
  same pattern table that drives `[S11]`). Single-end files derive
  their sample ID by stripping the `.(fastq|fq)(.gz|.bz2)?`
  extension (`[S21]`).
  - **Pass when:** for every row of the canonical table, the test
    `R1 filename → expected sample ID` produces the documented value;
    sample IDs are stable across compression variants of the same
    base name.
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
- `[S18]` required parameters (`forward_primer`, `reverse_primer`,
  `fastq_folder`) have no default and must be supplied via CLI or
  project config; the workflow aborts at startup with a message
  naming the missing parameter. `forward_primer` and `reverse_primer`
  are not required when `--no_trimming` is set (see `[S20]`).
  - **Pass when:** running the workflow without `--forward_primer`
    (or `--reverse_primer`, or `--fastq_folder`) exits non-zero and
    the stderr/log identifies the missing parameter; supplying the
    parameter via either `-params-file` or `--key value` lets the
    run proceed
- `[S19]` each Part A step emits a per-sample log file alongside its
  data output, published to `params.fastq_folder`:
    - merging       → `<sampleId>_merging.log`
    - trimming      → `<sampleId>_trimming.log` (only when the
      trimming step runs — see `[S20]`)
    - dereplicating → `<sampleId>_dereplicating.log`
    - clustering    → `<sampleId>_clustering.log`
  - **Pass when:** running Part A on any sample produces all four
    log files in `params.fastq_folder`, each non-empty (three when
    `--no_trimming` is set: no `_trimming.log`)
- `[S20]` `--no_trimming` toggle (default: `false`) skips the primer
  trimming step. The toggle and the primer parameters are mutually
  exclusive:
  - when `--no_trimming` is `true`, `forward_primer` and
    `reverse_primer` must be empty (not set on the CLI or in the
    config). The minimum-length and max-N filters that the trimming
    step used to enforce are taken over by
    `filter_and_convert_to_fasta` (vsearch `--fastq_minlen 32`,
    `--fastq_maxns <caller-supplied>`).
  - when `--no_trimming` is `false` (the default),
    `forward_primer` and `reverse_primer` are required (see
    `[S18]`).
  - **Pass when:** (a) running with `--no_trimming true` and empty
    primers succeeds and does **not** publish a `_trimming.log`;
    (b) running with `--no_trimming true` together with a non-empty
    `forward_primer` or `reverse_primer` exits non-zero and the
    error names the conflicting parameter.
- `[S21]` Part A collects every fastq file in the listed directories
  (`[S10]`/`[S11]`); a file whose R1 basename does not match any
  canonical pattern row nor the user-provided `--fastq_pattern`, or
  whose paired R2 partner is missing, is processed as a single-end
  sample by skipping `merge_fastq_pairs`. Single-end sample IDs come
  from stripping the `.(fastq|fq)(.gz|.bz2)?` extension; paired-end
  sample IDs are derived per `[S12]`.
  - **Pass when:** running Part A on a directory containing only an
    unpaired fastq file produces the expected per-sample artefacts
    (`<sampleId>.fas`, `_dereplicating.log`, `_clustering.log`) and
    the workflow trace does **not** mention `merge_fastq_pairs`.
- `[S22]` Part B's first step re-cleaves global swarm clusters by
  detecting alternative ("sub-") seeds that appear in a configurable
  fraction of samples (default 5 %, exposed as `--percentage`). For
  each global cluster that contains at least one sub-seed, the
  cluster is partitioned along the internal-structure father-son
  tree so each sub-seed inherits its reachable amplicons; clusters
  that contain no sub-seed are **not** re-emitted (the legacy
  pipeline concatenates the original and the cleaved triples
  downstream). The step consumes the four swarm outputs
  (`--seeds`, `--statistics-file`, `--internal-structure`,
  `--output-file`) plus a concatenated per-sample stats file, and
  emits an augmented `(stats, swarms, representatives)` triple. The
  fastidious flag used by the upstream swarm run is propagated to
  this step (default: `--fastidious`) so the representatives output
  filename matches the swarm parameters.
  - **Pass when:** golden-file characterization tests for
    `bin/cluster_cleaver.py` reproduce the byte-exact output of the
    legacy `tmp/OTU_cleaver.py` on a fixture covering: (a) a cluster
    that splits on a sub-seed, (b) a cluster that does not split,
    (c) a candidate sub-seed below the per-sample-presence threshold,
    (d) a sort tie within a sub-cluster.


## Common fastq file-name patterns

The pattern detector walks this table top-to-bottom and uses the
first row that matches the R1 basename. `<ext>` stands for
`(fastq|fq)(\.(gz|bz2))?`.

| #  | R1 pattern (anchored to end of basename) | Source / variants                                  |
|----|------------------------------------------|----------------------------------------------------|
| 1  | `_L00[1-9]_R1_00[1-9]\.<ext>`            | MiSeq default (e.g. `_L001_R1_001.fastq.gz`)       |
| 2  | `_L00[1-9]_.*_R1\.<ext>`                 | MiSeq with extra middle segment                    |
| 3  | `_L00[1-9]_R1\.<ext>`                    | MiSeq variant without trailing `_001`              |
| 4  | `[._][1-9]_1_.*\.<ext>`                  | numeric-lane with tail (e.g. `_1_1_junk.fastq.gz`) |
| 5  | `[._][1-9]_1\.<ext>`                     | numeric-lane (e.g. `_1_1.fastq.gz`)                |
| 6  | `[._]R1\.<ext>`                          | minimal `R1` (e.g. `_R1.fastq.gz`, `.R1.fastq.gz`) |
| 7  | `[._]1\.<ext>`                           | minimal `1` (e.g. `_1.fastq.gz`, `.1.fastq.gz`)    |

The R2 file name is derived by replacing the R1 discriminator with
its R2 counterpart (`R1`→`R2`, `_1_`→`_2_`, etc.) within the matched
span. The sample ID is the R1 basename **with the matched span and
everything after it stripped** — e.g. `A_L001_R1_001.fastq.gz` → `A`.

- `fastq` may also be `fq`; compression may be absent
- this table is the **single source of truth**; both the pattern
  matcher in `bin/discover_fastq.py` and the user-facing docs in
  [`README.md`](README.md) consult it

### Custom user pattern (`--fastq_pattern`)

`--fastq_pattern` lets users add a pattern when their file names
don't fit any canonical row. Format: a shell-style glob containing
the literal token `{1,2}` (the R1/R2 discriminator). The portion
before `{1,2}` is the sample-ID-bearing prefix; the portion after is
the suffix. Examples:

| `--fastq_pattern` glob               | Matches R1                       | Sample ID for `myrun_demoX_1.fq.gz` |
|--------------------------------------|----------------------------------|--------------------------------------|
| `*_{1,2}.fq.gz`                      | `*_1.fq.gz`                      | `myrun_demoX`                        |
| `*_run17_{1,2}.fastq.gz`             | `*_run17_1.fastq.gz`             | prefix before `_run17_`              |

A custom pattern is checked **before** the canonical table; if it
matches an R1 file, the canonical table is bypassed for that file.
