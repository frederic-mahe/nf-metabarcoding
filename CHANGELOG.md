# Changelog

All notable changes to `nf-metabarcoding` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/), and the
project adheres to [Semantic Versioning](https://semver.org/). The
version here must match `manifest.version` in
[`nextflow.config`](nextflow.config) and `version` in
[`CITATION.cff`](CITATION.cff) (enforced by `[S77]` / `[S80]`).

## [Unreleased]

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

### Added

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
