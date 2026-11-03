# Changelog

All notable changes to `nf-metabarcoding` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/), and the
project adheres to [Semantic Versioning](https://semver.org/). The
version here must match `manifest.version` in
[`nextflow.config`](nextflow.config) and `version` in
[`CITATION.cff`](CITATION.cff) (enforced by `[S77]` / `[S80]`).

## [0.8.0] - 2026-11-04

### Added

- `bin/prepare_silva.py`, a standalone tool that builds the stampa and
  sintax references from a SILVA SSU release (taxonomy, taxmap and
  fasta files). The pipeline does not call it. It was tested with
  SILVA 138.1, 138.2 and 144, and both of its outputs pass the
  reference check (`D21`, `[S119]`–`[S122]`). In the sintax output,
  placeholder names (`Incertae Sedis`, `uncultured`, ...) are prefixed
  with their nearest named ancestor, the organism becomes a species
  (`s:Genus_epithet`) only when it agrees with the lineage, and the
  remaining homonyms can be listed with `--homonyms`
  (`[S123]`–`[S125]`).
- Experimental support for Nanopore runs (`D23`). With
  `--pool_fastq_by_folder`, each subfolder of `--fastq_folder` (e.g.
  `fastq_pass/barcode01/`) becomes one single-end sample, and its
  fastq chunks are pooled into one file (`[S128]`, `[S129]`). Without
  the flag, discovery warns about a numbered file series (`_0`, `_1`,
  `_2`, ...) instead of pairing `_1` with `_2` silently (`[S130]`).
  An `--input` samplesheet may repeat a sample on several single-end
  rows to pool them (`[S132]`). `--fastq_recursive` scans subfolders
  at any depth (`[S131]`). The two flags are mutually exclusive and
  cannot be combined with `--input` (`[S133]`). Default settings stay
  tuned for Illumina reads.

### Changed

- Clearer task names in the log, trace and reports. The run-mode
  workflows in `main.nf` are now named after their input:
  `from_fastq`, `from_fasta`, `from_table` and `from_representatives`
  (were `part_a`, `part_b`, `part_c` and `part_c_fasta`), and Part C's
  internal `part_C_assign` is now `assign`. A task that used to appear
  as `part_a:part_C:part_C_assign:assign_taxonomy_stampa` now appears
  as `from_fastq:part_C:assign:assign_taxonomy_stampa`. Outputs are
  unchanged; site configs with `withName:` selectors on full task
  names need the new names (bare process names are unaffected).

### Removed

- Part B's `fake_taxonomic_assignment` and `fake_taxonomic_assignment2`
  processes (`[S33]`, `[S36]`, now retired). The occurrence table
  builder writes the placeholder identity / taxonomy / references
  (`0.0` / `NA` / `NA`) itself, and its `--assignments` option is
  now optional (`[S127]`). Part B output is unchanged, and each run
  has up to four fewer tasks.

### Fixed

- Part C standalone runs (and every other run without `--fastq_folder`)
  no longer log `WARN: Access to undefined parameter fastq_folder`.
  `fastq_folder`, `forward_primer` and `reverse_primer` are now
  declared `null` in `nextflow.config`, like the other input selectors;
  Part A still requires them (`[S18]`, `[S118]`).

## [0.7.0] - 2026-10-31

### Added

- Part C checks every header of its reference before any assignment
  (`[S114]`–`[S117]`, `D22`). A new `check_reference` task reads the
  whole file (plain, `.gz` or `.bz2`) with the rules of the method that
  uses it (stampa or sintax, shadow path included). Errors stop the run
  before the assignment: `format`, `invalid_utf8`, `control_char`,
  `duplicate_accession`, `comma_in_accession` (with
  `--majority_assignment`) and `whitespace` (extra or trailing spaces
  around a stampa lineage, which used to change results silently).
  Warnings (`empty_rank`, `star_in_lineage`) are printed and the run
  continues. Every finding goes to
  `logs/part_c/<basename>_reference_check.tsv`, published even when the
  run stops. Previously, a bad header after the first one crashed the
  run after the vsearch search, or silently changed the results. The
  check cannot be turned off; `-resume` reuses its result.
- `--majority_assignment` aborts at startup when the `--reference_dataset`
  accessions contain commas, as in a hybrid `>id;tax=...; lineage`
  header (`[S111]`). The majority step used to crash at the very end of
  the run.
- README: a "Reference datasets and the reference check" section
  describes the two reference formats, hybrid files and the check.
- Every `bin/*.py` helper now checks the Python version before anything
  else and, on Python < 3.10, exits with status 1 and a one-line message
  naming the script, the required and the found version, and suggesting
  an engine profile (`[S109]`). Previously such a run crashed with a
  `SyntaxError` or `TypeError` deep in the file. This happens when no
  container or conda profile is enabled and tasks use the host's
  `python3`, e.g. `-profile genotoul,demo`, where the compute nodes ship
  Python 3.6. The helpers no longer use
  `from __future__ import annotations`, which is itself a `SyntaxError`
  on Python < 3.7 and would hide the check.

### Changed

- An empty rank in a stampa lineage (`A||C`) is now specified: it stays
  an empty name, `*` when the top hits disagree, and the majority step
  counts it like any other name (`[S113]`). Behaviour is unchanged.
- `[S64]` now explains why two reference flags exist and documents
  hybrid references.

### Fixed

- `stampa_merge.py`, `update_occurrence_table.py` and
  `majority_assignment.py` no longer crash on reference headers that
  are not valid UTF-8 (e.g. a Latin-1 `é`); they pass those bytes
  through unchanged, whatever the locale (`[S112]`).
- `main.nf` failed to compile on Nextflow 26.10+ (`Return statement
  cannot be used in a workflow body`). The entry workflow now selects
  its run mode with an `if` / `else if` / `else` chain, and the Part A
  end-to-end path moves into a named `workflow part_a` ([S110]). A new
  CI check (`tests/check-compile.sh`) compiles `main.nf` on Nextflow
  25.10.2 and on the latest stable release.

## [0.6.0] - 2026-10-10

### Changed

- Require `vsearch` 2.32.0: the exact pins in `environment.yml` and the
  CI workflow (`.github/workflows/test.yml`) move from `2.31.0` to
  `2.32.0` ([S69]). The modules-profile examples in `README.md` and
  `conf/site.config.example` are bumped to match.
- Require `swarm` 3.2.0: the exact pins in `environment.yml` and the CI
  workflow move from `3.1.7` to `3.2.0` ([S69]), with the
  modules-profile examples in `README.md` and `conf/site.config.example`.
  swarm 3.2.0 fixes a data race in the fastidious phase that could
  silently drop light amplicons (making `--fastidious` results depend on
  thread timing), reports failed writes with a non-zero exit status,
  respects container memory limits, and is two to three times faster
  with `d` > 1.
- `merge_fastq_pairs`, `strip_reads` and `filter_and_convert_to_fasta`
  no longer pass `--fastq_qmax`: vsearch 2.32.0 defaults it to the
  highest quality the encoding can represent, which is the value the
  modules derived (`[S106]`). `merge_fastq_pairs` keeps its derived
  `--fastq_qmaxout`, whose default stays at 41 for `--fastq_mergepairs`.
- sintax classifications with a non-zero `--sintax_randseed` are now
  identical whatever the thread count (`[S91]`). vsearch 2.32.0 seeds a
  random substream per query, so bootstrap values, and occasionally the
  assignment retained at `--sintax_cutoff`, differ from those produced
  with vsearch 2.31.0 for the same seed.

### Fixed

- stampa (`assign_taxonomy_stampa`) now passes `--minwordmatches 0` to
  `vsearch --usearch_global` ([S108]). `--maxaccepts 0 --maxrejects 0`
  alone left vsearch's k-mer pre-filter on: a reference sharing fewer
  than 12 distinct 8-mers with the amplicon was never aligned, so short
  or divergent amplicons could miss their best reference or come out as
  `No_hit`. On 167 real 18S V4 amplicons against the V4 region of PR2
  5.0.0, 10 got a better top hit (5 of them previously `No_hit`), none
  got worse, and run time grew by about 4 %. Assignments can therefore
  change for such amplicons.
- `join_notmerged` wrote the join-padding quality as a fixed `I`, which
  is Q40 at `--fastq_encoding 33` but Q9 at `64`: the eight padding
  bases then added about 1.0 expected error per read, enough to drop
  every shadow cluster at the `max_ee` filter. The padding quality is
  now Q40 under either encoding (`h` at offset 64) ([S63]).

## [0.5.0] - 2026-07-23

### Changed

- Require `swarm` 3.1.7: the exact pins in `environment.yml` and the CI
  workflow (`.github/workflows/test.yml`) move from `3.1.6` to `3.1.7`
  ([S69]). The modules-profile examples in `README.md` and
  `conf/site.config.example` are bumped to match.

## [0.4.0] - 2026-07-10

### Added

- Optional, terminal post-mumu re-clustering pass for divergent markers
  (`[S102]`–`[S105]`, D20). `--recluster_id` (a real in `(0, 1]`, default
  `null` = OFF) is the master switch and the `vsearch --cluster_size`
  identity threshold; `--recluster_iddef` (int `[0, 4]`, default `2`)
  tunes the identity definition and is rejected without `--recluster_id`.
  When enabled, `recluster_search` runs abundance-based greedy clustering
  on the post-mumu FASTA and `recluster_merge`
  (`bin/recluster_otu_table.py`) folds each member OTU onto its centroid
  (samples + total summed, `spread` recomputed, metadata from the
  centroid, `cloud` left `NA`, OTUs renumbered `1..N`, read count
  asserted conserved). The reclustered table **replaces** Part B's
  emitted table and is fed to Part C, with a matching coarse
  `<basename>_table.fas` and a `<basename>_reclustering.log`. Default OFF
  keeps Part B's output byte-identical. Gated symmetrically on the shadow
  Part B path.

## [0.3.0] - 2026-07-02

### Fixed

- `--fastq_folder` / `--fasta_folder` now accept a Groovy list
  (`fastq_folder = ['a', 'b']`) as `[S10]` / `[S27]` document, not only
  a comma-separated string. The nf-schema migration had typed these
  params as `string`, so `validateParameters()` rejected an array
  (`Value is [array] but should be [string]`) before the discovery
  processes — which already handle a list — ever ran. The schema type is
  widened to `["string", "array"]`; the comma-separated string form is
  unaffected.

## [0.2.0] - 2026-06-30

### Added

- `--sintax_randseed` (`[S91]`, default `0`): seeds vsearch's PRNG for
  the sintax classifier so the bootstrap classification is reproducible
  across runs. The default `0` is vsearch's "use a random data source"
  sentinel, preserving the previous non-deterministic behaviour. Shared
  by the regular Part C sintax path and the shadow path.
- Stub blocks on every tool-invoking process (`[S85]`): the whole
  Part A→B→C pipeline now runs under `nextflow run -profile demo
  -stub-run` with none of vsearch / swarm / cutadapt / mumu installed,
  producing placeholder outputs in seconds. Validates the channel
  topology without tools or real data — a fast topology-CI job and an
  onboarding smoke check. The input-discovery / samplesheet processes are
  exempt (pure-Python glue that bootstraps the sample channel).
- Nextflow execution reports (`[S84]`): every run now writes
  `execution_timeline.html`, `execution_report.html`,
  `execution_trace.txt`, and `pipeline_dag.html` to
  `<outdir>/pipeline_info/`. The trace/report give per-task `peak_rss` /
  `realtime` against the requested CPU/memory — the primary input for
  sizing `--dataset_size_gb` / `--reference_size_gb` and per-step
  overrides on a new cluster.
- Documented an air-gapped / offline container path (`[S83]`): build the
  Wave image once on a connected node and reuse it from the shared
  container cache, or supply a pre-pulled `process.container` image via a
  `-c site.config` composed with `-profile slurm` (no Wave). See the
  README "Air-gapped clusters" section and `conf/site.config.example`.

### Fixed

- A Part A-only run (`--fastq_folder` without `--project_name`) now
  publishes `pipeline_info/software_versions.yml`, as `[S68]` requires
  on every entry point. The version dump was previously skipped outside
  the Part B / Part C paths.
- Every piped process script now runs under `set -euo pipefail`
  (`[S81]`), so a failure in any stage of a pipe fails the task instead
  of being masked by the last stage's exit status (e.g. a crashed
  `vsearch --fastx_filter` feeding `vsearch --uchime_denovo` in
  `chimera_detection` no longer silently truncates the `.uchime` table).
  Fixed a latent SIGPIPE abort in `chimera_detection_post_cleave`'s
  minimum-size computation (`sort -n | head -n 1` →
  `sort -n | sed -n '1p'`) that would bite large datasets, and guarded
  `fake_taxonomic_assignment`'s `grep` so a header-less input still
  yields an empty `.results`.

### Changed

- `run_mumu` now passes `--minimum_relative_cooccurrence` explicitly,
  coupled to the cleaving threshold `--percentage` (`[S43]` / `[S22]`)
  as `1 - percentage`. The default cleaving `0.05` yields `0.95` —
  mumu's own default — so behaviour at the default is unchanged, but the
  threshold now tracks `--percentage` and no longer silently inherits a
  future change to mumu's built-in default.
- **`cleanup` now defaults to `false`** (`[S82]`, was `true`). A
  successful run keeps its per-task `work/` directories, so `-resume`
  works across separate invocations and a run stays inspectable. Clean
  `work/` by hand when done, or set `cleanup = true` in a `-c` override
  for throwaway runs that should auto-reclaim it.
- Repository metadata made internally consistent (`[S80]`):
  `manifest.homePage` and the `nextflow_schema.json` `$id` now point at
  the canonical `frederic-mahe/nf-metabarcoding` repository on its
  `main` default branch (added `manifest.defaultBranch`), matching
  `CITATION.cff`.

## [0.1.0] - 2026-06-18

Initial release.

- Part A — per-sample processing (read merging, primer trimming,
  dereplication, quality extraction, local clustering with swarm), plus
  the opt-in experimental shadow pipeline for unmergeable pairs
  (`--recover_unmerged`).
- Part B — occurrence-table assembly (global dereplication + swarm
  clustering, cluster cleaving, chimera detection, substring-OTU
  merging, mumu post-clustering curation).
- Part C — taxonomic assignment (stampa scatter-gather and the sintax
  shadow path), with optional majority-rule assignment.
- `--input` samplesheet input and the unified `--outdir` output layout.
- Execution profiles: `slurm`, `conda`, `modules`, container engines
  (`docker` / `podman` / `singularity` / `apptainer` via Seqera Wave),
  and a self-contained `demo` profile.
- `nf-schema` parameter validation driven by `nextflow_schema.json`.
- Multi-layer test suite (nf-test, bats, pytest) with shellcheck /
  flake8 linting and a spec ↔ test coverage gate.
