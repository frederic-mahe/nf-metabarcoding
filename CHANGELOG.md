# Changelog

All notable changes to `nf-metabarcoding` are documented here. The
format follows [Keep a Changelog](https://keepachangelog.com/), and the
project adheres to [Semantic Versioning](https://semver.org/). The
version here must match `manifest.version` in
[`nextflow.config`](nextflow.config) and `version` in
[`CITATION.cff`](CITATION.cff) (enforced by `[S77]` / `[S80]`).

## [Unreleased]

### Fixed

- A Part A-only run (`--fastq_folder` without `--project_name`) now
  publishes `pipeline_info/software_versions.yml`, as `[S68]` requires
  on every entry point. The version dump was previously skipped outside
  the Part B / Part C paths.

### Changed

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
