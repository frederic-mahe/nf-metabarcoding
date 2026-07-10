# Changelog

All notable changes to `nf-metabarcoding` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/), and the
project adheres to [Semantic Versioning](https://semver.org/). The
version here must match `manifest.version` in
[`nextflow.config`](nextflow.config) and `version` in
[`CITATION.cff`](CITATION.cff) (enforced by `[S77]` / `[S80]`).

## [Unreleased]

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
