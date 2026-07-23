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
    3. **Part C** — taxonomic assignment (stampa or sintax,
       `[S47]`–`[S51]`): update the occurrence table
  - **Pass when:** running the full pipeline on a paired-end fixture
    produces, in order, per-sample `.fas` (Part A), an occurrence
    table (Part B), and a taxonomy-annotated occurrence table (Part C)
- `[S02]` each part can be run separately, or all at once. The
  entry point is selected implicitly by which input param is set
  (no explicit `--only=` flag): `--occurrence_table` runs table-input
  Part C standalone (and the shadow path when its preconditions hold,
  see [S62]); `--representatives_fasta` runs fasta-input Part C
  standalone ([S48]); `--fasta_folder` (or an `--input` samplesheet in
  the fasta profile, [S70]) runs Part B standalone (with Part C chained
  on when `--reference_dataset` is set); `--fastq_folder` (or an
  `--input` samplesheet in the fastq profile) runs Part A end-to-end
  (with Part B / Part C chained on when `--project_name` /
  `--reference_dataset` are set); `--accession` runs the standalone
  fetch mode ([S97]).
  The six input-mode selectors `--occurrence_table` /
  `--representatives_fasta` / `--input` / `--fasta_folder` /
  `--fastq_folder` / `--accession` are **mutually exclusive**: setting
  more than one aborts at startup with a message listing the ones that
  were set, rather than silently picking one and ignoring the rest.
  (This supersedes the earlier "first match wins" dispatch and subsumes
  the `--input`-vs-folder exclusivity of [S70].)
  - **Pass when:** running with only `--fastq_folder` set invokes
    every Part A process and no Part B / Part C process; the existing
    Part B standalone (`--fasta_folder`) and Part C standalone
    (`--occurrence_table`) tests cover the other two cases; setting two
    mode selectors (e.g. `--occurrence_table` with `--fastq_folder`)
    aborts at startup naming both.
- `[S03]` fastq files can be paired-end or single-end, compressed
  (`.gz`, `.bz2`) or not
  - **Pass when:** the same input data in any of these forms produces
    byte-identical Part A outputs
- `[S04]` when processing paired-end fastq files, reads that can be
  merged are processed normally. The handling of reads that **cannot**
  be merged is gated by `--recover_unmerged` ([S78], default `false`):
  when it is unset the not-merged reads are dropped (no shadow
  artefacts, the default), and only when it is set are they collected
  via `vsearch --fastq_mergepairs --fastqout_notmerged_fwd` and
  `--fastqout_notmerged_rev` (reads stay in sync) and routed through a
  parallel **shadow** Part A pipeline:
    1. join R1/R2 with `vsearch --fastq_join`, passing
       `--join_padgap` set to a run of `A`s of length
       `params.join_padding_length` (default 8; see [S63]) and
       `--join_padgapq` set to the same number of `I` quality
       characters (Phred 40). `A` is used instead of vsearch's
       default `N` so the joined sequence carries only A/C/G/T —
       swarm accepts the padding without complaint and no mask /
       restore round-trip is needed. The shadow path has no merging
       step (by definition the reads in this branch could not be
       merged) and therefore **no `_merging.log`** is published for
       shadow samples — only the per-step logs from later stages
       (trimming, dereplicating, clustering) appear under
       `<sampleId>_notmerged_*.log`;
    2. trim primers (when `--no_trimming` is false) and convert to
       fasta with `vsearch --fastx_filter --fastq_maxns 0`; the
       A-padding contributes no `N`s, so the shadow path uses the
       same max-N budget as the regular path;
    3. dereplicate and extract ee values per the normal pipeline; the
       published `.fas` retains the run of `A`s at the join site.
       The run of `A`s is **artificial padding**, not biological
       sequence — it is the user's responsibility to keep that in
       mind when interpreting shadow-pipeline output (the shadow
       pipeline is experimental).
  - shadow-pipeline sample IDs are `<sampleId>_notmerged`; published
    artefacts are `<sampleId>_notmerged.{fas,qual,stats}` (data, under
    `<outdir>/per_sample/`) plus the per-step logs from the stages that
    actually run in the shadow path (under `<outdir>/logs/part_a/per_sample/`,
    [S71]/D16): `<sampleId>_notmerged_trimming_forward.log` and
    `<sampleId>_notmerged_trimming_reverse.log` (the two cutadapt
    passes, see `[S19]`), `<sampleId>_notmerged_dereplicating.log`,
    and `<sampleId>_notmerged_clustering.log`. No
    `<sampleId>_notmerged_merging.log` is produced — see step 1
    above.
  - **Pass when:** running Part A on a paired-end fixture whose reads
    cannot overlap **with `--recover_unmerged true`** produces non-empty
    `<sampleId>_notmerged.{fas,qual,stats}`
    in `<outdir>/per_sample/`, the published `.fas` contains a run of
    `params.join_padding_length` `A`s (default 8), no shadow output
    sequence line contains `N`, the trimming / dereplicating /
    clustering shadow logs are present and non-empty in
    `<outdir>/logs/part_a/per_sample/`, and no
    `<sampleId>_notmerged_merging.log` exists. Running the same fixture
    **without** `--recover_unmerged` (the default) produces no
    `_notmerged` artefacts and runs no `join_notmerged` / `strip_reads`
    process ([S78]).
- `[S78]` `params.recover_unmerged` (boolean, default `false`) is the
  master switch for the experimental shadow pipeline ([S04]). When
  `false` (the default) the whole shadow path is off: Part A drops the
  not-merged reads (no `join_notmerged` / `strip_reads`, no
  `<sampleId>_notmerged.*` artefacts), and the three entry points
  (end-to-end Part A→B→C, Part B standalone, Part C standalone) never
  invoke `part_B_shadow` / `part_C_shadow` — so no `_notmerged`
  occurrence table or taxonomy is produced, and no empty shadow
  artefacts appear when every pair merges. When `true` the shadow path
  runs as described in [S04] / [S56] / [S50], subject to its other
  preconditions (e.g. `--reference_dataset_sintax` for shadow Part C,
  [S64]). The `[S23]` reserved-suffix guard is **independent** of this
  flag (a sample whose ID ends in `notmerged` is always rejected, so it
  cannot collide with shadow naming if the flag is later enabled). The
  `stripright` ([S24]) and `join_padding_length` ([S63]) knobs only take
  effect when this flag is set.
  - **Pass when:** an end-to-end run on a paired-end fixture whose reads
    cannot merge produces `_notmerged` per-sample artefacts only when
    `--recover_unmerged true` is set; the same run at the default
    (unset) produces none and schedules no `join_notmerged` /
    `strip_reads` / `part_B_shadow` process.
- `[S05]` unmerged-pair clusters appear in the occurrence table with a
  per-sample marker (working name: `sampleID_partial`)
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) — final marker name


## Interfaces between parts

Each part is a directory contract so parts can be tested in
isolation.

- **Part A → Part B**
  - `<sampleId>.fas` — dereplicated FASTA, vsearch
    `--sizeout`-style abundance annotation, `--fasta_width 0`
  - `<sampleId>.qual` — three space-separated columns
    `sha1 ee length`, sorted by length (numeric asc), then SHA1
    (ascii asc), then ee (numeric asc). One row per unique SHA1 —
    `uniq --check-chars=40` collapses duplicates onto the lowest-ee
    row.
  - `<sampleId>.stats` — swarm `--statistics-file` output, filtered
    to clusters with > 2 reads
- **Part B → Part C**
  - one occurrence table (see schema below)
- **Part C output**
  - the occurrence table with taxonomy column updated (no rows
    added or removed)


## Occurrence table schema

The Part B occurrence table (`<basename>_table.tsv`, [S46]) and Part C's
annotated sibling (`<basename>_table_assigned.tsv`, [S51]) carry 13
metadata columns followed by one column per sample, in this order:

| #  | Column       | Notes                                                                                          |
|----|--------------|------------------------------------------------------------------------------------------------|
| 1  | `OTU`        | renumbered 1..N after mumu ([S44])                                                             |
| 2  | `total`      | sum of the sample columns; size=0 → 1 hotfix ([S44]) so vsearch `--sizein` consumers can read it |
| 3  | `cloud`      | swarm `--statistics-file` cloud size; forced to `"NA"` after mumu ([S44])                      |
| 4  | `amplicon`   | SHA1 seed ID; the join key against `.qual` / `.distr`                                          |
| 5  | `length`     | seed sequence length                                                                           |
| 6  | `abundance`  | seed abundance                                                                                 |
| 7  | `chimera`    | uchime `Y` / `N` / `NA`; forced to `"N"` after mumu ([S44])                                    |
| 8  | `spread`     | count of samples with `abundance > 0`; recomputed at [S39] / [S44]                             |
| 9  | `quality`    | ee/length ratio from [S28]'s per-sample `.qual` merge                                          |
| 10 | `sequence`   | seed nucleotides, single-line                                                                  |
| 11 | `identity`   | top-hit identity %, or `0.0`; placeholder until Part C runs                                    |
| 12 | `taxonomy`   | LCA lineage, or `NA`; placeholder until Part C runs                                            |
| 13 | `references` | top-hit accessions, or `NA`; placeholder until Part C runs                                     |
| 14–N | `<sampleId>` | integer read count per sample, sorted; empty samples ([S09]) contribute a zero-filled column  |

Part C never adds or removes columns or rows; it only overwrites
columns 11–13 (`identity`, `taxonomy`, `references`) on rows whose
amplicon appears in the taxonomy assignments.

### Two-table mode

When `--split-occurrence-table` is set (default `false`, see [S15]),
the single-table `<basename>_table.tsv` is **replaced** by a pair of
sibling files under `--results_folder`:

- `<basename>_clusters.tsv` — columns 1–13 above (every metadata
  column; no per-sample columns). One row per cluster.
- `<basename>_occurrences.tsv` — long format with three columns:
  `OTU`, `sampleId`, `abundance`. One row per non-zero
  `<OTU, sampleId, abundance>` triple — zero-abundance cells and
  empty samples ([S09]) contribute no rows. The `OTU` column is the
  join key against `<basename>_clusters.tsv`.

Only one of `_table.tsv` / (`_clusters.tsv` + `_occurrences.tsv`) is
published per run; the [S59] whitelist picks the active set. The
split applies to both the regular Part B path ([S46]) and the shadow
Part B path ([S56]); the basename carries the `_notmerged` token in
the latter case.


## Workflow requirements

- `[S06]` configurable via a config file or command-line parameters
  - Not part of automated CI: the workflow accesses every user-facing
    value through `params.X`, a hashmap Nextflow populates from
    `nextflow.config`, `-params-file params.yml`, and `--key value`
    flags indistinguishably. There is nothing in the workflow code to
    test — the equivalence is structurally guaranteed upstream. `[S18]`
    asserts the contract explicitly for required parameters
    ("supplying the parameter via either `-params-file` or
    `--key value` lets the run proceed"); that is the meaningful gate.
- `[S07]` runs locally or on HPC with slurm
  - Not part of automated CI; covered by manual smoke tests on the
    cluster. Profiles live in `nextflow.config`.
  - the `slurm` profile sets a `resourceLimits` ceiling
    (`params.max_cpus` / `params.max_memory` / `params.max_time`,
    defaults `16` / `128.GB` / `240.h`) so the per-retry escalation
    (`memory = N.GB * task.attempt` across the tiers) is clamped to the
    largest node a partition offers instead of requesting more than
    exists and looping to failure. `resourceLimits` is native Nextflow
    behaviour (>= 24.04), so — per `[S00]`'s "we do not re-test
    upstream tools" rule — the clamp itself is not unit-tested; the
    manual cluster smoke test confirms it.
  - the `slurm` profile submits at most `submitRateLimit` jobs per minute
    (default `50/1min`) and supports an **opt-in** `--slurm_array_size`:
    when set to `N`, `process.array = N` batches up to `N` ready tasks of
    a process into one `sbatch --array` submission instead of one sbatch
    per task — far less scheduler load for the per-sample fan-out, and the
    pattern HPC admins recommend at scale. Default `null` keeps one job
    per task (unchanged behaviour); a cluster profile may set a default
    (`genotoul` uses 50, `[S87]`).
    - **Pass when:** `nextflow config -profile slurm` resolves
      `process.array = null` (arrays off by default) and a non-empty
      `executor.submitRateLimit` (`tests/check-slurm-config.sh`).
- `[S87]` ships per-cluster *institutional* profiles so a known HPC site
  runs with no hand-written config. Each lives in
  `conf/clusters/<name>.config` and is exposed as a same-named profile in
  `nextflow.config`; that profile `includeConfig`s `conf/slurm.config`
  first (so it carries the full `[S07]` executor + resource tiers) then
  layers the site's overrides — partition/account routing, a
  `resourceLimits` ceiling clamped to that cluster's largest node, and
  the singularity/apptainer bind mounts + image cache. A cluster profile
  therefore *implies* slurm: a user selects a cluster **and** a
  dependency engine, e.g. `-profile meso,singularity` — `slurm` is never
  listed.
  - Configs are **vendored** (copied into the repo and pinned), not
    fetched from nf-core/configs at runtime, because the pipeline
    supports air-gapped compute nodes (`[S83]`) and an auditable, pinned
    config is the hard-to-misuse default. The shipped set: `abims`,
    `genotoul`, `ifb_core`, `meso`, `saga` (`saga` values are sourced
    from the Sigma2/NRIS docs, to be confirmed by a cluster smoke run).
    `conf/clusters/_template.config` documents the knobs for adding one.
  - The profile **wiring** is checked automatically by
    `tests/check-cluster-profiles.sh` (each cluster profile resolves the
    slurm executor + its `resourceLimits` ceiling via `nextflow config`,
    and composes with an engine profile). Per `[S00]`, actual submission
    on each cluster stays a manual smoke test — the partition/account
    *values* can only be confirmed on the hardware.
  - **Pass when:** `nextflow config -profile <cluster>` resolves
    `process.executor = 'slurm'` and the cluster's `resourceLimits`
    ceiling for each shipped cluster; `-profile <cluster>,singularity`
    additionally resolves `singularity.enabled = true`; and a plain
    `nextflow run` (no profile) does not set the slurm executor.
- `[S92]` a cluster profile may mark the slurm account as **mandatory**
  by setting `params.require_slurm_account = true`. On a site where every
  job must be charged to a project (e.g. ABiMS for the projects that
  enforce it), submitting without `-A` is rejected by the scheduler, so
  the workflow fails fast: when `require_slurm_account` is set and
  `--slurm_account` is empty/null, it aborts at startup naming the flag
  rather than letting every `sbatch` bounce mid-run. The flag defaults to
  `false`, so the generic `-profile slurm` and any cluster that leaves it
  unset keep `--slurm_account` optional ([S75]); the `abims` profile sets
  it `true`. The pure helper
  `slurm_account_requirement_error(require_account, slurm_account)`
  returns the error message (or `null` when the requirement is satisfied
  or not in force); the entry workflow throws when it is non-null.
  - **Pass when:** `slurm_account_requirement_error` unit tests show a
    non-null message naming `--slurm_account` when the requirement is on
    and the account is unset/blank, and `null` when the account is given
    or the requirement is off; and `nextflow config -profile abims`
    resolves `params.require_slurm_account = true`.
- `[S79]` under `-profile slurm`, the memory-bound steps scale their RAM
  request off `--dataset_size_gb` (the dataset-bound steps) and
  `--reference_size_gb` (the taxonomic-assignment steps), falling back to
  fixed defaults when those hints are unset ([S07]). Because a forgotten
  hint silently under-provisions a large run — surfacing only as an OOM
  kill deep in the pipeline — the workflow emits a **startup warning**
  (it does not abort) for each unset hint: one for `--dataset_size_gb`,
  and, when a reference is in use, one for `--reference_size_gb`. The
  warning fires only under `-profile slurm` (the resource tiers do not
  apply otherwise) and points the user at the flag (or at overriding the
  step's memory in a `-c site.config`, [S75]). The pure helper
  `resource_size_warnings(profile, dataset_size_gb, reference_size_gb,
  reference_in_use)` returns the messages; the entry workflow reads
  `workflow.profile` + the params and prints them.
  - **Pass when:** `resource_size_warnings` unit tests show: under
    `slurm` with `--dataset_size_gb` unset it returns a warning naming
    `dataset_size_gb`; with it set, none; the reference warning appears
    only when a reference is in use and `--reference_size_gb` is unset;
    a non-`slurm` profile returns no warnings regardless.
- `[S08]` runs tools either from the local environment (PATH /
  `-profile conda` / `-profile modules`) or inside containers. Four
  engine profiles — `docker`, `podman`, `singularity`, `apptainer` —
  each enable their engine plus Seqera Wave, which builds the image on
  the fly from `environment.yml` (the single pinned source of truth,
  `[S69]`) and caches it; no Dockerfile or registry is maintained
  (`DECISIONS.md` D10). The profiles compose with the executor and
  dependency profiles (e.g. `-profile slurm,singularity`). Wave needs
  outbound network from wherever tasks run.
  - The profile **wiring** is checked automatically by
    `tests/check-container-profiles.sh` (each profile resolves its
    engine + `wave.enabled` + the pinned conda env via
    `nextflow config`, and a plain `nextflow run` stays bare-PATH).
    Container **execution** is not run in automated CI — per `[S00]`'s
    "we do not re-test upstream tools" rule, running vsearch/swarm in a
    container exercises Wave/the engine, not our glue — so it stays a
    manual `-profile docker` / `-profile singularity` cluster smoke
    test.
  - **Pass when:** `nextflow config -profile <engine>` resolves
    `<engine>.enabled` and `wave.enabled` for each of docker / podman /
    singularity / apptainer; `-profile slurm,singularity` resolves the
    slurm executor alongside the engine; and a plain `nextflow run`
    (no profile) enables neither Wave nor conda.
- `[S09]` empty input samples must travel through and appear in the
  occurrence table (but not in the two-table mode long-format)
  - **Pass when:** an empty fastq pair produces a row in the
    occurrence table with `abundance = 0` for that sample (and no
    error)
- `[S10]` accepts a directory, or a list of directories (absolute or
  relative paths). `--fastq_folder` accepts a single path or a
  comma-separated list (`--fastq_folder a,b,c`); a `nextflow.config`
  may instead supply a Groovy list (`fastq_folder = ['a', 'b']`).
  Relative paths are resolved against the Nextflow launch directory
  (`launchDir`, i.e. the directory `nextflow run` is invoked from) —
  not against the config file's location — following standard Nextflow
  `file()` semantics.
  - **Pass when:** every fastq file in every listed folder is
    discovered; identical inputs spread across two folders produce
    identical artefacts to the single-folder run; a relative
    `fastq_folder` is discovered when its folder exists relative to
    `launchDir`.
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
- `[S67]` rejects a `--fastq_pattern` whose `{<r1>,<r2>}` brace token
  has identical sides (e.g. `{1,1}`, `{R1,R1}`). The two sides are the
  R1/R2 discriminator and must differ; were they equal, the derived
  R2 basename would equal the R1 basename and the R1 file would be
  paired with itself. The check fires at pattern-parse time, before
  any file is globbed, with a message naming the offending token.
  - **Pass when:** a pattern with equal sides (`*_{1,1}.fastq.gz`)
    raises an error at parse time; a pattern with differing sides
    (`{1,2}`, `{R1,R2}`) is accepted.
- `[S96]` `--fastq_pattern` may contain at most a small fixed number of
  `*` wildcards (3). Each `*` translates to a greedy `.*` in the
  discovery regex; an unbounded count of them lets a crafted pattern
  drive catastrophic regex backtracking — a denial-of-service hang — on
  the scanned file names. A legitimate pattern needs one `*` to capture
  the sample prefix (occasionally a second), so the cap never blocks a
  real pattern. A pattern exceeding the cap is rejected at pattern-parse
  time, before any file is globbed, with a message stating the limit.
  - **Pass when:** a `--fastq_pattern` carrying more than three `*`
    (e.g. `*a*a*a*a{1,2}.fastq`) raises at parse time naming the limit;
    patterns with up to three `*` are accepted.
- `[S12]` derives the sample ID for a paired-end pair from the R1
  basename by stripping the **first matching pattern** suffix (the
  same pattern table that drives `[S11]`). Single-end files derive
  their sample ID by stripping the `.(fastq|fq)(.gz|.bz2)?`
  extension (`[S21]`).
  - **Pass when:** for every row of the canonical table, the test
    `R1 filename → expected sample ID` produces the documented value;
    sample IDs are stable across compression variants of the same
    base name.
- `[S13]` aborts when two or more discovered samples share a
  derived sample ID; the error message names the colliding sample
  ID and lists every input file path involved in the collision
  (one group per ID). The check runs **before** any Part A / Part B
  process starts.
  - **Pass when:** running Part A with two fastq inputs that both
    derive to sample ID `A` exits non-zero, stderr names `A` and
    both fastq paths, and no Part A process appears in the
    nextflow trace.
- `[S14]` collision policy for same-named samples: **refuse**
  (sample IDs must be unique). Resolved per `D03` in
  [`DECISIONS.md`](DECISIONS.md); implemented by `[S13]`.
  - **Pass when:** `bin/discover_fastq.py` (Part A) and
    `bin/discover_fasta.py` (Part B) each expose a
    `check_unique_sample_ids()` helper that raises on duplicates,
    and each CLI exits non-zero with a stderr message listing every
    duplicate file path.
- `[S15]` users can choose between a single occurrence table and a
  two-part table by setting `--split-occurrence-table` (default
  `false`). When the flag is set, Part B publishes
  `<basename>_clusters.tsv` and `<basename>_occurrences.tsv` under
  `--results_folder` and **does not** publish
  `<basename>_table.tsv` — the split pair replaces the single table;
  the [S59] whitelist swaps accordingly. The two forms carry the
  same information modulo layout: every non-zero cell of the
  single-form table corresponds to exactly one row of
  `_occurrences.tsv` with matching `OTU` / `sampleId` / `abundance`
  values. The split applies to both the regular and the shadow
  ([S56]) Part B paths. Part C's output shape is unaffected by this
  flag for now (it continues to emit a single
  `<basename>_table_assigned.tsv` regardless) — see
  [`DECISIONS.md`](DECISIONS.md) D05 sub-question 6.
  - **Pass when:** an end-to-end Part B run with
    `--split-occurrence-table true` publishes
    `<basename>_clusters.tsv` and `<basename>_occurrences.tsv` (and
    no `<basename>_table.tsv`); the same run without the flag
    publishes only `<basename>_table.tsv` (no split pair); every
    non-zero cell of the single-form table appears as exactly one
    row in `_occurrences.tsv`.
- `[S16]` expects demultiplexed fastq files
  - **Pass when:** documented in README; demultiplexing is out of
    scope (could be added later as a subworkflow)
- `[S17]` per-cluster minimum-read threshold: `list_local_clusters`
  keeps only swarm clusters with strictly more than
  `--min_cluster_size` reads (default `2`, the legacy `> 2 reads`
  rule). The threshold reaches `filter_swarm_stats.awk` via
  `-v min_cluster_size=`; the awk script defaults to `2` when the
  variable is unset so it stays usable standalone. Validated as an
  integer `>= 0` by the schema.
  - **Pass when:** (a) with the default, `list_local_clusters` emits
    no row with `reads <= 2` on a fixture that contains a singleton
    cluster; (b) raising `--min_cluster_size` additionally drops the
    rows at or below the new threshold (a stale hard-coded `> 2`
    would keep them).
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
- `[S19]` every Part A per-sample artefact is published under
  `<outdir>` ([S71]; `--fastq_folder` is input-only). Data files go to
  `<outdir>/per_sample/`; the per-step logs go to the parallel
  `<outdir>/logs/part_a/per_sample/` (D16):
    - data files (`<outdir>/per_sample/`): `<sampleId>.fas`,
      `<sampleId>.qual`, `<sampleId>.stats`
    - per-step log files (`<outdir>/logs/part_a/per_sample/`):
        - merging       → `<sampleId>_merging.log` (regular path
          only — the shadow path has no merging step, see `[S04]`)
        - trimming      → `<sampleId>_trimming_forward.log` and
          `<sampleId>_trimming_reverse.log` (only when the trimming
          step runs — see `[S20]`). `trim_primers` runs cutadapt in
          two passes — the forward primer first, then the reverse
          primer — and each pass writes its own report so the
          per-primer trimming statistics stay separable.
        - dereplicating → `<sampleId>_dereplicating.log`
        - clustering    → `<sampleId>_clustering.log`
  - **Pass when:** running Part A on any sample produces all three
    data files under `<outdir>/per_sample/` and all five log files
    under `<outdir>/logs/part_a/per_sample/`, each non-empty (three logs when
    `--no_trimming` is set: no `_trimming_forward.log` /
    `_trimming_reverse.log`); no `*.log` is left in
    `<outdir>/per_sample/`.
- `[S20]` `--no_trimming` toggle (default: `false`) skips the primer
  trimming step. The toggle and the primer parameters are mutually
  exclusive:
  - when `--no_trimming` is `true`, `forward_primer` and
    `reverse_primer` must be empty (not set on the CLI or in the
    config). The minimum-length and max-N filters that the trimming
    step used to enforce are taken over by
    `filter_and_convert_to_fasta` (vsearch `--fastq_minlen
    <--fastq_minlen, default 32>` per `[S90]`, `--fastq_maxns 0`).
  - when `--no_trimming` is `false` (the default),
    `forward_primer` and `reverse_primer` are required (see
    `[S18]`).
  - **Pass when:** (a) running with `--no_trimming true` and empty
    primers succeeds and does **not** publish a `_trimming.log`;
    (b) running with `--no_trimming true` together with a non-empty
    `forward_primer` or `reverse_primer` exits non-zero and the
    error names the conflicting parameter.
- `[S88]` `trim_primers` passes `--error-rate <--primer_error_rate>`
  to both cutadapt passes (default `0.1`, cutadapt's own default).
  Validated as a real in `[0, 1]` by the schema. Only consulted when
  trimming runs.
  - **Pass when:** running `trim_primers` with a non-default
    `--primer_error_rate` echoes that value on cutadapt's
    `Command line parameters:` line in the per-sample forward
    trimming log.
- `[S89]` `trim_primers` sets cutadapt `--overlap` per primer to
  `int(primer_length * --primer_overlap_fraction)` (default `0.667`,
  reproducing the legacy `2/3` rule for the primer lengths seen in
  practice). Validated as a real in `(0, 1]` by the schema. Only
  consulted when trimming runs.
  - **Pass when:** running `trim_primers` with a non-default
    `--primer_overlap_fraction` echoes the computed `--overlap N` on
    cutadapt's `Command line parameters:` line in the per-sample
    forward trimming log.
- `[S90]` `filter_and_convert_to_fasta` passes `--fastq_minlen
  <--fastq_minlen>` to vsearch (default `32`, the legacy value):
  reads shorter than this (after trimming) are dropped before the
  fasta conversion. Validated as an integer `>= 1` by the schema.
  - **Pass when:** raising `--fastq_minlen` above every read length
    in a fixture yields an empty filtered fasta, while the default
    keeps the reads (a stale hard-coded `32` would keep them).
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
- `[S23]` `notmerged` is a reserved sample-ID suffix used by the
  shadow pipeline (`[S04]`). A sample whose ID ends with the literal
  string `notmerged` (e.g. `X_notmerged`) is rejected at discovery
  time, before any merging step runs, so shadow-pipeline artefacts
  cannot collide with user-supplied sample IDs.
  - **Pass when:** running Part A on a folder containing
    `<X_notmerged>_{1,2}.<ext>` (or any single-end variant) exits
    non-zero and the error message names the reserved keyword.
- `[S93]` a sample ID must contain only the safe character set
  `[A-Za-z0-9._-]` and must start with a letter, digit, or
  underscore. The ID becomes a shell token and an output-file
  basename in the Part A / Part B process scripts (`!{sampleId}`), so
  characters outside this set — whitespace, path separators, shell
  metacharacters (`$` `` ` `` `;` `|` `&` `(` `)` `<` `>` `*` `?`
  quotes …), a leading `-` (read as an option by downstream tools) or
  a leading `.` (hidden file / `..` traversal) — are rejected before
  any process starts. The rule applies to IDs from an `--input`
  samplesheet (`[S70]`) and to IDs derived from folder discovery
  (`[S12]`), and is enforced by a single shared validator so the two
  entry paths cannot diverge. It is independent of and additional to
  the uniqueness (`[S13]` / `[S14]`) and reserved-suffix (`[S23]`)
  checks.
  - **Pass when:** `bin/parse_samplesheet.py` on a samplesheet whose
    `sample` cell is `bad;id` (or contains a space, `/`, `$(...)`, or
    a leading `-`) exits non-zero naming the offending ID; and
    `bin/discover_fastq.py` / `bin/discover_fasta.py` on a folder
    holding a file that derives such an ID exit non-zero naming it.
    All three import the validator from the shared `bin/sample_id.py`
    module.
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
  emits an augmented `(stats, swarms, representatives)` triple.

  `params.fastidious` (default `true`) is the single propagation
  lever for the fastidious flag: when `true`, both the upstream
  `global_clustering` ([S32]) swarm run and the `cleaving` cluster
  cleaver run receive `--fastidious`, and every artefact carries
  the `_1f` swarm-parameters suffix (`_1f.swarms`, `_1f.stats`,
  `_1f.struct`, `_1f_representatives.fas`, `_1f.swarms2`,
  `_1f.stats2`, `_1f_representatives.fas2`). When `false`, neither
  process receives the flag, and the cleaver's representatives
  filename swaps from `_1f_representatives.fas2` to
  `_1_representatives.fas2` (the `.stats2`/`.swarms2` filenames
  follow the input naming, so they inherit the upstream `_1f`/`_1`
  suffix chosen by `global_clustering`). The toggle exists for
  very large datasets where `--fastidious` is too memory- or
  time-expensive; the default preserves the higher recall of the
  fastidious pass.
  - **Pass when:** golden-file characterization tests for
    `bin/cluster_cleaver.py` reproduce the byte-exact output of the
    legacy `tmp/OTU_cleaver.py` on a fixture covering: (a) a cluster
    that splits on a sub-seed, (b) a cluster that does not split,
    (c) a candidate sub-seed below the per-sample-presence threshold,
    (d) a sort tie within a sub-cluster; an end-to-end run with
    `--fastidious false` completes successfully and `cleaving`
    publishes `<basename>_1_representatives.fas2` (not the `_1f`
    variant).
- `[S24]` the shadow pipeline ([S04]) strips `params.stripright`
  nucleotides from the 3' end of each not-merged R1 and R2 read
  before joining, using `vsearch --fastx_filter --fastq_stripright`.
  The parameter is user-configurable and defaults to 30. Setting it
  to 0 disables the trim (vsearch runs as a no-op pass-through). The
  strip step sits between `merge_fastq_pairs` and `join_notmerged`,
  so the join padding falls inside higher-quality bases on both
  sides. The regular (merged) Part A path is unaffected.
  - **Pass when:** running `strip_reads` on a fastq with 70 nt reads
    produces 40 nt reads with the default `stripright = 30`, and 70 nt
    reads with `stripright = 0`.
- `[S25]` Part B requires `--project_name` (no default). The value
  is used as the filename prefix for every project-wide artefact
  (global dereplicated fasta, quality file, distribution file,
  per-sample-stats file, swarm outputs). The workflow aborts at
  startup if Part B is asked to run without `--project_name`. The
  parameter is **only** required for Part B; Part A and Part C do
  not consume it.
  - **Pass when:** running Part B without `--project_name` exits
    non-zero and stderr names the missing parameter; supplying it
    via `-params-file` or `--key value` lets the run proceed and
    every Part B artefact filename begins with the supplied value.
- `[S26]` **Superseded by `[S71]` (D09).** Part B's output location is
  now `--outdir` (default `results`), so no output-folder parameter is
  required — Part B no longer aborts when `--results_folder` is unset.
  `--results_folder` survives as a deprecated alias for `--outdir`
  ([S71]). Nextflow's `publishDir` materialises `<outdir>/<subdir>`
  (and any missing parents) on first publish; the workflow performs no
  filesystem I/O at parse time (D08).
  - **Pass when:** running Part B with neither `--outdir` nor
    `--results_folder` succeeds and publishes under `results/`
    (the default); `--outdir <d>` (non-existent, absolute or relative)
    succeeds, creates `<d>`, and every Part B artefact lands under
    `<d>/occurrence_table/` ([S71]).
- `[S27]` Part B builds its fasta channel by collecting every
  `.fas` file that satisfies one of:
    - produced by Part A in the same run (when Parts A and B run
      end-to-end), or
    - globbed from the list of `--fasta_folder` directories
      (Part B standalone; same single-path / comma-list / Groovy-list
      semantics as `[S10]` for `--fastq_folder`).

  Files whose basename ends in `_notmerged.fas` are excluded
  ([S04]'s shadow-pipeline artefacts are processed downstream by a
  dedicated path, not by the regular Part B pipeline). **Empty
  `.fas` files are kept** so that empty samples travel through to
  the occurrence table (`[S09]`) — they contribute a sample column
  filled with zeros and never a row. The Part B sample ID is the
  basename stripped of the `.fas` extension. Duplicate sample IDs
  are an error: the workflow aborts before any Part B process runs
  and the stderr lists every duplicated `.fas` path (see `[S13]`).
  - **Pass when:** running Part B on a folder containing `A.fas`,
    `A_notmerged.fas`, `B.fas`, and an empty `C.fas` builds a fasta
    channel of `{A.fas, B.fas, C.fas}` (notmerged excluded, empty
    sample retained); running on two folders that both contain
    `A.fas` exits non-zero and stderr lists both paths.
- `[S28]` Part B's `build_expected_error_file` merges every
  per-sample `<sampleId>.qual` from the `[S27]` fasta channel into
  a single project-wide quality file
  `<project_name>_<N>_samples.qual` (`N` = number of `.fas` in the
  channel). The merge keeps one row per SHA1 (the input rows are
  pre-sorted by length then SHA1, and `uniq --check-chars=40`
  collapses duplicates).
  - **Pass when:** running on two `.qual` fixtures whose SHA1 sets
    overlap emits one row per distinct SHA1 and no SHA1 appears
    twice; the output filename is
    `<project_name>_<N>_samples.qual`.
- `[S29]` Part B's `build_distribution_file` scans the FASTA
  headers of every `.fas` in the `[S27]` channel and writes the
  sequence ↔ sample map to `<project_name>_<N>_samples.distr`. Each
  row is `<sha1>\t<sampleId>\t<size>` (tab-separated); the
  `;size=N` annotation is stripped from the header.
  - **Pass when:** running on a fixture with samples `A` (two
    records) and `B` (one record) emits three rows, sample IDs
    match the fasta basenames, and the `<size>` column is the
    integer value of `;size=N`.
- `[S30]` Part B's `list_all_cluster_seeds_of_size_greater_than_2`
  concatenates every per-sample `.stats` (the swarm-stats files
  published by Part A, filtered to clusters > 2 reads per `[S17]`)
  into a single project-wide
  `<project_name>_<N>_samples_per_sample_OTUs.stats`. Each row is
  the original `.stats` row prefixed with `<sampleId>\t` (sample ID
  derived from the `.stats` basename).
  - **Pass when:** running on two `.stats` fixtures emits
    `len(rows_A) + len(rows_B)` rows; the first column carries the
    sample ID; the remaining columns preserve the original swarm
    `--statistics-file` layout.
- `[S31]` Part B's `global_dereplication` concatenates every
  `.fas` in the `[S27]` channel and runs
  `vsearch --fastx_uniques` with `--sizein --sizeout --fasta_width
  0`. The output is `<project_name>_<N>_samples.fas`; the vsearch
  log is published as `<project_name>_<N>_samples.log`.
  - **Pass when:** running on two `.fas` fixtures that share a
    sequence yields one record per distinct sequence; the
    `;size=N` value is the sum of the per-sample sizes.
- `[S32]` Part B's `global_clustering` runs swarm on the
  globally-dereplicated fasta from `[S31]` with `--differences 1
  --usearch-abundance` and the four output flags
  `--internal-structure`, `--output-file`, `--statistics-file`,
  `--seeds`. `--fastidious` is propagated from `params.fastidious`
  (see [S22]). The swarm-parameters suffix `<sfx>` is `1f` when
  `params.fastidious` is `true` (the default) and `1` when it is
  `false`. Output filenames follow the
  `<project_name>_<N>_samples_<sfx>.{swarms,stats,struct}` and
  `<project_name>_<N>_samples_<sfx>_representatives.fas` scheme;
  the run log is published as `<project_name>_<N>_samples_clustering.log`
  ([S45], the project-wide construct — independent of `<sfx>`).
  - **Pass when:** the four output files plus the log exist and are
    non-empty for the documented fixture; a run with
    `--fastidious false` produces `_1.{swarms,stats,struct}` /
    `_1_representatives.fas` (no `_1f` artefacts).
- `[S33]` Part B's `fake_taxonomic_assignment` writes a placeholder
  Part C input file by scanning the swarm representatives FASTA
  (`<basename>_1f_representatives.fas` from `[S32]`) and emitting,
  per record, the tab-separated row
  `<sha1>\t<size>\t0.0\tNA\tNA`. The output path is
  `<basename>_1f_representatives.results`. This step exists so the
  occurrence-table builder can run before the real taxonomic
  assignment lands (Part C); it is later overwritten by the
  stampa/sintax output.
  - **Pass when:** running on a representatives FASTA with three
    records emits three TSV rows, one per record, in the
    `<sha1>\t<size>\t0.0\tNA\tNA` shape.
- `[S34]` Part B's `chimera_detection` filters the global swarm
  representatives to abundance `>= params.chimera_minsize` (default
  2) with `vsearch --fastx_filter --minsize`, then pipes the kept
  records into `vsearch --uchime_denovo`. Outputs are
  `<basename>_1f_representatives.uchime` (the uchime hit table)
  and `<basename>_1f_representatives.log` (stderr of the uchime
  run).
  - **Pass when:** both output files exist and the `.log` is
    non-empty for the documented fixture (the `.uchime` table can
    legitimately be empty when no chimeras are found).
- `[S36]` Part B's `fake_taxonomic_assignment2` mirrors `[S33]` on
  the cleaver's representatives FASTA (`<basename>_1f_representatives.fas2`
  from `[S22]`). Output: `<basename>_1f_representatives.results2`.
  Empty input (no clusters got cleaved) yields an empty output
  file — downstream concatenation tolerates that.
  - **Pass when:** running on a fas2 with `K` records emits `K`
    placeholder rows; running on an empty fas2 emits an empty
    file.
- `[S37]` Part B's `chimera_detection_post_cleave` re-runs uchime_denovo on
  the concatenation of pre-cleave and cleaved representatives
  (`<basename>_1f_representatives.fas` from `[S32]` plus
  `<basename>_1f_representatives.fas2` from `[S22]`). The
  `--minsize` filter is dynamic: it drops to the smallest size
  observed in the cleaved file (so newly cleaved low-abundance
  clusters are still searched for chimeras), but never goes below
  `params.chimera_minsize` (default 2). If the cleaved file is
  empty, `--minsize` falls back to `params.chimera_minsize`. The
  concatenated stream is sorted by size before uchime sees it.
  Outputs are `<basename>_1f_representatives.uchime2` and
  `<basename>_1f_representatives.log2`.
  - **Pass when:** both output files exist and the `.log2` is
    non-empty for the documented fixture (the `.uchime2` table can
    legitimately be empty when no chimeras are found).
- `[S35]` Part B's `build_occurrence_table` runs
  `bin/build_filtered_contingency_table.py` to merge the
  swarm representatives, swarm stats, swarm output, uchime hits,
  per-amplicon quality (`build_expected_error_file`'s
  `[S28]` output), taxonomic assignments, and the sequence ↔
  sample distribution (`build_distribution_file`'s `[S29]` output)
  into a single occurrence table. Rows are filtered: a cluster is
  kept iff `chimera == "N" && ee/length <= --max_ee && (abundance >=
  --min_abundance || spread >= --min_spread)`. The three thresholds
  reach the script via `--max-ee` / `--min-abundance` / `--min-spread`
  and default to the legacy `0.0002` / `3` / `2` (`--max_ee` a real
  `>= 0`, `--min_abundance` and `--min_spread` integers `>= 1`).
  Sample columns appear in sorted order; empty samples (`[S09]`)
  contribute a zero-filled column.
  - **Pass when:** golden-file characterization tests for
    `bin/build_filtered_contingency_table.py` reproduce its byte-exact
    stdout on a fixture covering: (a) a cluster passing every
    filter, (b) each filter individually rejecting a cluster,
    (c) the `abundance >= 3` and `spread >= 2` alternatives,
    (d) a uchime row missing column 18 → status "NA",
    (e) a partial uchime row (`IndexError` on column 1),
    (f) duplicate sample rows in `.distr` summing onto the seed,
    (g) a sample present only in a filtered-out cluster surviving
    as a zero column, (h) `#` characters stripped from taxonomy
    strings.


- `[S38]` Part B's `search_for_terminal_gaps` extracts the
  representative sequences from the post-`build_occurrence_table`
  OTU table and runs `vsearch --cluster_smallmem --id 1.0
  --qmask none --usersort` to find OTU pairs that are 100% identical
  modulo terminal gaps (sub- and super-strings). Output: the `^H`
  lines from vsearch's `--uc` stream (the hits). The legacy bash
  name for this step is
  `extract_fasta_and_search_for_identical_sequences`; the workflow
  uses the shorter name.
  - **Pass when:** running on a small OTU table produces a non-empty
    `.uc`-style stream containing one row per identical-modulo-gaps
    OTU pair; rows start with the literal `H`.
- `[S39]` Part B's `merge_substring_otus` runs
  `bin/merge_substring_otus.py` to merge
  pupil OTUs into their masters (sample columns summed,
  ``spread`` recomputed from non-zero merged columns, ``total``
  summed, ``cloud`` incremented by ``pupil_cloud + 1`` per merged
  pupil), then sorts the resulting table by the OTU column and
  asserts that the total read count is conserved.
  - **Pass when:** golden-file characterization tests for
    `bin/merge_substring_otus.py` reproduce
    byte-exact output on a fixture covering: (a) a pass-through OTU,
    (b) a master with a single pupil, (c) a master with two pupils,
    (d) an overlap where one OTU is both master and pupil → script
    exits non-zero and prints WARNING to stderr.
- `[S40]` Part B exposes two FASTA-extraction processes that both
  read column 4 as the amplicon ID, column 2 as the abundance,
  and column 10 as the sequence (header
  `<amplicon>;size=<abundance>;`):
    - `extract_otu_fasta` — keeps every data row (`awk 'NR > 1'`)
      and is used before the mumu pass to feed
      `find_similar_sequences`;
    - `extract_mumu_fasta` — adds an `$2 != 0` filter so rows
      whose `total` is zero are dropped; used after `rebuild_post_mumu_table`
      to produce the project's final post-mumu fasta. The legacy
      bash uses `extract_fasta_sequences_from_occurrence_table`
      and `extract_fasta_sequences_from_occurrence_table2`
      respectively; these are the workflow's shorter names.
  - **Pass when:** the pre-mumu FASTA contains one record per data
    row (zero-abundance rows kept verbatim); the post-mumu FASTA
    contains one record per non-zero row.
- `[S41]` Part B's `trim_metadata_for_mumu` reduces the OTU table
  to two slices that mumu accepts: the amplicon column
  (column 4) followed by every sample column (columns 14 onward).
  - **Pass when:** the output has the amplicon column as column 1
    and the sample columns preserved verbatim from column 2 onward.
- `[S42]` Part B's `find_similar_sequences` runs
  `vsearch --usearch_global` self-search on the OTU FASTA with the
  legacy lulu-recommended parameters (`--id <--similar_id> --iddef 1
  --maxaccepts 0 --query_cov <--similar_query_cov> --maxhits
  <--similar_maxhits> --userfields query+target+id`). The three
  exposed thresholds default to `0.84` / `0.9` / `10` (`--similar_id`
  and `--similar_query_cov` reals in `(0, 1]`, `--similar_maxhits` an
  integer `>= 1`). The output is the userout stream with the
  `;size=N;` annotation stripped from every column.
  - **Pass when:** (a) the output is a 3-column TSV
    (`query\ttarget\tid`) with no `;size=` annotations; (b) tightening
    `--similar_id` to `1.0` drops the partial-identity matches the
    default `0.84` keeps.
- `[S43]` Part B's `run_mumu` invokes the `mumu` binary
  (`>= 1.1.1`) with `--otu_table`, `--match_list`,
  `--new_otu_table`, `--minimum_relative_cooccurrence`, and `--log`.
  The cleaned-up intermediate inputs (`_reduced.table`, `.match_list`)
  are not kept. `--minimum_relative_cooccurrence` is **not** an
  independent parameter: it is coupled to the cleaving threshold
  `--percentage` (`[S22]`) as `1 - percentage`, so the default cleaving
  `0.05` yields `0.95` (which also matches mumu's built-in default).
  The two thresholds are complementary — cleaving keeps a sub-seed that
  appears in at least `percentage` of samples, and mumu merges a child
  OTU only when it co-occurs with its parent in at least `1 -
  percentage` of the child's samples. Passing it explicitly also pins
  the value rather than inheriting a future change to mumu's default.
  - **Pass when:** the `_raw_mumu.table` and the `.mumu.log` are
    produced and the log is non-empty.
- `[S44]` Part B's `rebuild_post_mumu_table` runs
  `bin/rebuild_table_after_mumu.py` to splice the per-amplicon
  metadata (`length`, `abundance`, `quality`, `sequence`,
  `identity`, `taxonomy`, `references`) from the pre-mumu OTU
  table back onto every row mumu emits, then renumbers OTUs
  starting at 1 (`cloud` becomes `"NA"`; `chimera` is forced to
  `"N"`). A downstream awk hotfix replaces `total == 0` with `1`
  to satisfy `vsearch --sizein` consumers.
  - **Pass when:** golden-file characterization tests for
    `bin/rebuild_table_after_mumu.py` reproduce byte-exact output
    on a fixture covering: (a) header passthrough, (b) per-row
    metadata join keyed on amplicon ID, (c) new-`total` /
    new-`spread` recomputation, (d) a zero-abundance row surviving
    the join (size=0 → 1 hotfix is the bash wrapper's
    responsibility).
- `[S45]` Part B publishes one step-level log file per major step into
  `<outdir>/logs/part_b/` ([S71]/D16),
  named with the project-wide basename `<project>_<N>_samples`
  (the same construct as the final occurrence table, see `[S46]`)
  and a step suffix. The six step logs are:
    - `<basename>_dereplication.log` — `global_dereplication`
      vsearch log (`[S31]`)
    - `<basename>_clustering.log` — `global_clustering` swarm log
      (`[S32]`)
    - `<basename>_chimera_detection.log` — the concatenation of
      both chimera-detection runs' stderr: the pre-cleave
      `chimera_detection` (`[S34]`) followed by the post-cleave
      `chimera_detection_post_cleave` (`[S37]`). Each fragment is
      preceded by a one-line section header
      (`=== chimera_detection ===` /
      `=== chimera_detection_post_cleave ===`) so the two runs
      remain distinguishable
    - `<basename>_cleaving.log` — `cleaving`'s stderr from
      `bin/cluster_cleaver.py` (`[S22]`)
    - `<basename>_superstring_clustering.log` — combined stderr of
      `search_for_terminal_gaps` + `merge_substring_otus`
      (`[S38]`, `[S39]`)
    - `<basename>_post_clustering_curation.log` — `run_mumu`'s
      `--log` output (`[S43]`)
  - **Pass when:** running Part B on the documented fixture
    publishes the six step logs above into
    `<outdir>/logs/part_b/` ([S71]/D16); each file is
    non-empty; none of them appears in `<outdir>/occurrence_table/`.
- `[S46]` Part B publishes its final occurrence table as
  `<project>_<N>_samples_table.tsv` (after `rebuild_post_mumu_table`
  and the size=0→1 awk hotfix). Among the OTU tables produced
  along the way, only this final `_table.tsv` is published; the
  intermediate tables (`*.OTU.filtered.cleaved.table`,
  `*.nosubstringOTUs.table`, `*_raw_mumu.table`, `*.mumu.table`)
  stay in the work directory. Non-table artefacts (`.qual`,
  `.distr`, `_per_sample_OTUs.stats`, the global `.fas`, swarm
  outputs, `.uchime` hit tables, extracted `.fas` files) continue
  to be published.
  - **Pass when:** running Part B on the documented fixture
    publishes `<project>_<N>_samples_table.tsv` to
    `<outdir>/occurrence_table/` ([S71]); the file is non-empty and
    starts with the OTU-table header row; none of the intermediate
    `.OTU.filtered.cleaved.table`, `.nosubstringOTUs.table`,
    `_raw_mumu.table`, or `.mumu.table` files are present in the
    output.


### Optional post-mumu re-clustering (divergent markers)

For very variable / fast-mutating marker genes (e.g. viral
protein-coding markers) swarm `d=1` ([S32]) over-splits and mumu
([S43]) leaves a large tail of near-identical variants. An **opt-in,
terminal** coarse clustering pass (vsearch `--cluster_size`, abundance-
based greedy clustering) lumps that tail into broader OTUs, reducing
the dataset before Part C. It is **not** a replacement for swarm and is
never layered inside the curation chain — it runs last, on the fully-
curated post-mumu table. Default OFF; see D20 for the rationale and the
locked decisions D-a…D-d.

- `[S102]` The re-clustering pass is gated on a single master switch,
  `--recluster_id` (vsearch `--id`, a real in `(0, 1]`, default
  `null` = OFF). A second knob `--recluster_iddef` (vsearch `--iddef`,
  int `[0, 4]`, default `2`) tunes the identity definition for markers
  with indels / length variation. The `--recluster_iddef` knob is inert
  without `--recluster_id`: setting it to a non-default value while
  `--recluster_id` is unset aborts at startup naming `--recluster_id`,
  rather than silently ignoring the value.
  - **Pass when:** `--recluster_iddef` set to a non-default value with
    no `--recluster_id` exits non-zero and stderr names
    `--recluster_id`; the schema pins `recluster_id` to a real in
    `(0, 1]` and `recluster_iddef` to an int in `[0, 4]` with
    default `2`.
- `[S103]` When `--recluster_id` is set, Part B's `recluster_search`
  runs `vsearch --cluster_size` on the post-mumu FASTA
  (`extract_mumu_fasta`'s output, [S40], header `>amplicon;size=total;`)
  with `--sizein --sizeout --id <--recluster_id>
  --iddef <--recluster_iddef> --qmask none --fasta_width 0 --strand plus
  --maxaccepts 0 --maxrejects 0 --threads <task.cpus>` and captures the
  `^H` lines of the `--uc` stream (member → centroid connexions). The
  fixed flags (`--cluster_size`, `--sizein`/`--sizeout`, `--qmask none`,
  `--maxaccepts 0`/`--maxrejects 0`) are not user-facing; only `--id`
  and `--iddef` are exposed ([S102]). vsearch's `--log` captures the
  step's stderr.
  - **Pass when:** running `recluster_search` on a small FASTA whose
    records are identical modulo a chosen threshold produces a
    `.uc`-style stream whose rows all start with `H`; an empty stream
    (no member below the centroid) is a legitimate green outcome.
- `[S104]` Part B's `recluster_merge` runs `bin/recluster_otu_table.py`
  to fold each member OTU into its vsearch centroid: sample columns
  summed, `total` summed, `spread` recomputed from the non-zero merged
  columns, per-OTU metadata (`amplicon`, `length`, `abundance`,
  `quality`, `sequence`, `identity`, `taxonomy`, `references`) taken
  from the centroid (the most abundant member; D-c), `cloud` left at the
  post-mumu `"NA"` and `chimera` at `"N"` (no arithmetic — post-mumu
  they are already strings). Surviving centroid rows are emitted in
  input order and renumbered `1..N` (D-d). The `;size=N;` annotation is
  stripped from the `.uc` labels before matching on the amplicon column.
  As in [S39] the merge aborts non-zero if any OTU is both a centroid
  and a member, and asserts the total read count is conserved.
  - **Pass when:** golden-file characterization tests for
    `bin/recluster_otu_table.py` reproduce byte-exact output on a
    fixture covering: (a) a pass-through singleton, (b) a centroid with
    a single member, (c) a centroid with two members, (d) a member row
    that precedes its centroid in the table, (e) `cloud` staying `"NA"`,
    (f) contiguous `1..N` renumbering and read-count conservation; and a
    match list where one OTU is both centroid and member aborts
    non-zero with a WARNING on stderr.
- `[S105]` When `--recluster_id` is set the reclustered table
  **replaces** the post-mumu table as Part B's emitted output (D-a): it
  is published as `<basename>_table.tsv`, a coarse post-mumu FASTA
  sibling `<basename>_table.fas` is re-extracted from it (so table and
  FASTA stay consistent), and both the fine pre-recluster
  `rebuild_post_mumu_table` table and its FASTA stay unpublished in the
  work directory. Part C then runs on the reduced set. The combined
  vsearch-search + merge stderr is published as
  `<basename>_reclustering.log` under `<outdir>/logs/part_b/` (a
  conditional 7th step log alongside the six of [S45]). The
  [S59] `occurrence_table/` whitelist is unchanged — the two published
  filenames (`_table.tsv`, `_table.fas`) are the same, only their
  contents are the coarse set. The gating applies symmetrically to the
  shadow Part B path ([S56], `_notmerged` basename). When
  `--recluster_id` is unset (default) no re-clustering process is
  scheduled and Part B's output is byte-identical to today.
  - **Pass when:** a Part B run with `--recluster_id` unset schedules no
    `recluster_*` process and its `occurrence_table/` output is
    unchanged; a run with `--recluster_id` set publishes a
    `<basename>_table.tsv` with fewer data rows than the pre-recluster
    table, a matching `<basename>_table.fas`, a
    `<basename>_reclustering.log` under `logs/part_b/`, and Part C
    consumes the reduced table.


### Per-sample read/cluster tracking

- `[S107]` Part B publishes a per-sample **read/cluster tracking**
  summary `<basename>_read_counts.tsv` to `<outdir>/logs/part_b/`
  ([S71]/D16) — the Part B counterpart of Part A's `[S86]` summary.
  Part B is a **pooled** pipeline (every sample is merged at
  `global_dereplication`, [S31]), so its step logs ([S45]) carry only
  project-wide totals and cannot yield a per-sample view; the per-sample
  dimension survives instead in the abundance data. The summary is
  therefore reconstructed by `bin/build_part_b_read_counts.py` from the
  sequence-to-sample distribution ([S29] `.distr`) and the per-sample
  columns (14…N of the occurrence-table schema) of the intermediate OTU
  tables — **not** from the logs. One row per sample (sorted), then a
  `Total` row; every cell is a non-negative integer. The columns are:
    - `reads_in` — reads entering Part B, the sum of the sample's
      `.distr` sizes ([S29]);
    - `reads_kept` — reads surviving the occurrence-table filter
      ([S35]), the sum of the sample's column in `build_occurrence_table`'s
      output. Reads are conserved by every later curation step
      ([S39]/[S44]/[S104]), so this is also the sample's final read
      count and the per-sample loss `reads_in - reads_kept` is
      concentrated at the [S35] filter;
    - `clusters_kept` / `clusters_merged` / `clusters_mumu` — the number
      of clusters the sample appears in (cells `> 0`) after the filter
      ([S35]), after substring-merging ([S39]), and after mumu ([S44]);
    - `clusters_recluster` — appended **only** when the optional
      re-clustering pass runs (`--recluster_id`, [S105]): the same count
      on the reclustered table. Absent otherwise.
  The `Total` row sums each column: for the read columns this is the
  project read count, for the cluster columns the sum of per-sample
  cluster memberships (i.e. total `<cluster, sample>` incidences at that
  stage), not the number of distinct clusters. Empty samples ([S09])
  contribute no `.distr` rows, so they surface as an all-zero row — the
  authoritative sorted sample list is threaded in (as for [S35]'s
  `--samples`) so the row still appears. The summary is produced on both
  the regular Part B path and the shadow path ([S56], `_notmerged`
  basename); it is always-on (no gating flag).
  - **Pass when:** golden-style tests for
    `bin/build_part_b_read_counts.py` on a fixture covering (a) a sample
    losing reads at the [S35] filter (`reads_kept < reads_in`), (b) a
    sample whose cluster count drops at substring-merge and again at
    mumu, (c) an empty sample as an all-zero row, (d) a correct `Total`
    row, and (e) the `clusters_recluster` column present with
    `--recluster` and absent without it; and a process test publishes
    `<basename>_read_counts.tsv` under `logs/part_b/` with the documented
    header.


## Part C — taxonomic assignment

Part C re-implements the legacy
[`stampa`](./tmp/stampa/) method on top of nextflow. Inputs: either
the Part B occurrence table (`[S46]`) or a standalone fasta file
of representative sequences; plus a reference dataset (fasta,
optionally compressed). Output: an occurrence table whose
`taxonomy` / `identity` / `references` columns have been updated
from placeholder values to real taxonomic assignments.

- `[S47]` Part C requires a reference fasta whose **format
  matches the selected assignment method** ([S61]). The value is
  an absolute or relative path to a fasta file (plain
  `.fasta`/`.fas`, gzip-compressed `.gz`, or bzip2-compressed
  `.bz2`). Two flags expose the two formats:
    - `--reference_dataset` — stampa-formatted reference (header
      `>id space-separated lineage`); consumed by the regular
      path when `--taxonomy_method=stampa` (default).
    - `--reference_dataset_sintax` — sintax-formatted reference
      (header `>id;tax=d:Dom,p:Phy,...;`); consumed by the
      regular path when `--taxonomy_method=sintax` and by the
      shadow Part C path ([S50] / [S64]).
  The workflow aborts at startup if Part C is asked to run with
  the method-appropriate flag missing. Both flags are **only**
  required for Part C; Parts A and B do not consume them.
  - **Pass when:** running Part C with `--taxonomy_method=stampa`
    (default) but no `--reference_dataset` exits non-zero and
    stderr names `reference_dataset`; running Part C with
    `--taxonomy_method=sintax` but no `--reference_dataset_sintax`
    exits non-zero and stderr names `reference_dataset_sintax`;
    supplying the matching flag lets the workflow proceed.
- `[S48]` Part C accepts either an occurrence table (the
  `[S46]` `_table.tsv`, via `--occurrence_table`) or a representatives
  fasta (via `--representatives_fasta`) as its primary input; the two
  flags are mutually exclusive ([S02]).
  - **Table input** (`--occurrence_table`): the process
    `extract_fasta_sequences_from_occurrence_table` reads column 4
    (amplicon ID), column 2 (abundance), and column 10 (sequence)
    and emits a fasta of representatives with header
    `<amplicon>;size=<abundance>;`; the assignment is then spliced
    back onto the table as `<basename>_table_assigned.tsv` ([S51]).
  - **Fasta input** (`--representatives_fasta`): the extraction step
    and the occurrence-table join are skipped. The supplied fasta is
    fed straight to the assignment selected by `--taxonomy_method`
    ([S61]), and the sole deliverable is the standalone
    `<basename>_taxonomy_<method>.tsv` that the assignment step
    publishes ([S49]/[S61]) — `_taxonomy_stampa.tsv` (5-column,
    header) for stampa, the 4-column `_taxonomy_sintax.tsv` for
    sintax. No `_table_assigned.tsv` is produced. `<basename>` is
    derived from the fasta filename (extension stripped). The shadow
    path ([S50]) does not run in this mode, and `--majority_assignment`
    is rejected at startup ([S66]: no occurrence table to compute a
    per-OTU majority on). Both methods require the reference matching
    `--taxonomy_method` ([S47]/[S64]).
  - **Pass when:** a `--representatives_fasta` run with
    `--reference_dataset` (stampa) publishes
    `<basename>_taxonomy_stampa.tsv` and **no**
    `<basename>_table_assigned.tsv`; the same run with
    `--taxonomy_method=sintax` + `--reference_dataset_sintax`
    publishes the 4-column `<basename>_taxonomy_sintax.tsv`; setting
    both `--representatives_fasta` and `--occurrence_table` aborts at
    startup naming the input-mode selectors; combining
    `--representatives_fasta` with `--majority_assignment` aborts at
    startup naming `majority_assignment`.
- `[S49]` Part C's primary taxonomic-assignment path is a port of
  the legacy `stampa.sh` / `stampa_merge.py` pipeline, expressed
  as a Nextflow scatter-gather:
    1. the representatives fasta is split into chunks of
       `params.stampa_chunk_size` sequences via Nextflow's
       `splitFasta(by: N, file: true)` operator. The sentinel
       value `0` (set by the `local` profile) bypasses the
       split — the full fasta becomes a single chunk. The
       slurm-tuned default is `1000`.
    2. each chunk feeds `vsearch --usearch_global` against the
       reference dataset with `--top_hits_only --output_no_hits
       --maxaccepts 0 --maxrejects ${params.stampa_maxrejects}
       --id ${params.stampa_id} --notrunclabels --rowlen 0`
       and `--userfields query+id<iddef>+target`. The two
       tunables expose the vsearch flags directly:
       `params.stampa_maxrejects` (default `0` = vsearch's
       "no limit" sentinel: exhaustive search, every reference
       seq is considered as a potential top hit) and
       `params.stampa_id` (default `0.5`, identity threshold
       below which a hit is dropped — `stampa_merge.py` then
       emits `No_hit` for the amplicon). Each chunk is then
       parsed by `bin/stampa_merge.py`, which computes the
       last-common-ancestor taxonomy across top hits per
       amplicon. Each chunk emits a `stampa_chunk.tsv` slice
       with columns
       `amplicon\tabundance\tidentity\ttaxonomy\treferences`
       (same shape as `fake_taxonomic_assignment`'s `[S33]`
       output).
    3. the per-chunk slices are concatenated (via Nextflow's
       `collectFile`) and then sorted by the `sort_taxonomy`
       process running `LC_ALL=C sort -k2,2nr -k1,1d` (abundance
       desc, amplicon asc — the legacy stampa order), which
       publishes the merged `<basename>_taxonomy_stampa.tsv` to
       the results folder. That published table carries a header
       row with the `[S46]` column names
       (`amplicon\tabundance\tidentity\ttaxonomy\treferences`); the
       per-chunk `stampa_chunk.tsv` slices stay headerless. The sort
       is a real process rather than
       `collectFile`'s `sort:` closure because that closure orders
       whole entries (chunks), not the lines inside a chunk, and so
       cannot stabilise vsearch's thread-dependent within-chunk
       order.
     4. each chunk's `vsearch.log` is emitted alongside its
        `stampa_chunk.tsv` slice and the slices are gathered
        (via `collectFile`) into a single `<basename>_taxonomy.log`
        published to `<outdir>/logs/part_c/` ([S59]/D16), the stampa
        counterpart of the sintax path's `--log` artefact. No
        separate merge process.
  - **Pass when:** the scatter-gather smoke tests in
    `tests/main.nf.test` (one with `params.stampa_chunk_size=1`
    that forces a split, one with `=0` for the no-split branch)
    publish the expected `<basename>_taxonomy_stampa.tsv` sorted
    by abundance desc then amplicon asc; the per-chunk process
    test in `tests/processes/part_c/assign_taxonomy_stampa.nf.test`
    confirms one chunk's vsearch + `bin/stampa_merge.py` wiring and
    that it emits a `vsearch.log`; the subworkflow test in
    `tests/subworkflows/part_c.nf.test` confirms the stampa path
    gathers the per-chunk logs into a published
    `logs/part_c/<basename>_taxonomy.log`; and LCA correctness is
    pinned at the helper level by `tests/python/test_stampa_merge.py`.
- `[S50]` Part C's **shadow path** runs `vsearch --sintax`
  against `--reference_dataset_sintax` ([S64]) on the shadow
  occurrence table `<basename>_notmerged_table.tsv` produced by
  shadow Part B ([S56]) and emits
  `<basename>_notmerged_table_assigned.tsv` (the shadow sibling
  of `[S51]`'s `<basename>_table_assigned.tsv`). The shadow path
  has no method toggle — sintax is the only assignment method
  for the shadow path. The vsearch call uses
  `--sintax_cutoff ${params.sintax_cutoff}` (default `0.9`),
  `--dbmask none`, and `--tabbedout`; the assignment lifted into
  the occurrence table's `taxonomy` column is **column 2** of
  the tabbed output (the bootstrap-annotated lineage, e.g.
  `d:Bacteria(0.99),p:Proteobacteria(0.85)`). `params.sintax_cutoff`
  is passed to vsearch for completeness but only affects vsearch's
  internal column 4 (the cutoff-filtered lineage); since the workflow
  lifts column 2 (unfiltered, with per-rank bootstrap), the value is
  forward-compatible config — visible in the vsearch log but not in
  the published `taxonomy` column. The `identity` and
  `references` columns are not populated by sintax — they stay
  at their `[S33]`/`[S46]` placeholder values (`0.0` / `NA`).
  The shadow Part C workflow (`part_C_shadow`) is invoked
  alongside the regular `part_C` from all three entry points
  (end-to-end Part A→B→C, Part B standalone with a reference,
  and Part C standalone — see [S61]) **only when
  `--recover_unmerged` is set ([S78]) and
  `--reference_dataset_sintax` is set** ([S64]). When
  `--recover_unmerged` is unset (the default), the
  sintax-formatted reference is missing, or no shadow
  occurrence table is available (no `_notmerged` sample upstream,
  or no sibling file in standalone mode), the shadow Part C
  workflow is not invoked — silently, no error.
  - **Pass when:** an end-to-end run on a paired-end fixture
    that produces at least one `_notmerged` sample publishes
    both `<basename>_table_assigned.tsv` (stampa, regular path)
    and `<basename>_notmerged_table_assigned.tsv` (sintax,
    shadow path) under `--results_folder`; the shadow assigned
    table's `taxonomy` column carries the bootstrap-annotated
    lineage from vsearch's tabbed output.
- `[S91]` `assign_taxonomy_sintax` passes `--randseed
  <--sintax_randseed>` to `vsearch --sintax` (default `0`, vsearch's
  "use a random data source" sentinel — i.e. non-reproducible, the
  legacy behaviour). A non-zero value seeds vsearch's PRNG so the
  bootstrap classification is reproducible across runs. The same param
  feeds both the regular Part C sintax path
  (`--taxonomy_method=sintax`) and the shadow path (`[S50]`), since
  both share the one `assign_taxonomy_sintax` process. Validated as an
  integer `>= 0` by the schema (`[S72]`). Full thread-count-independent
  reproducibility additionally depends on the vsearch version; the flag
  is accepted by the pinned vsearch and is forward-compatible config.
  - **Pass when:** a non-default `--sintax_randseed` value is echoed by
    vsearch as `--randseed <value>` in the process `--log` (a dropped
    param or a hard-coded seed would leave the default in the log).
- `[S51]` Part C's `update_occurrence_table` splices the
  taxonomic assignment back onto the `[S46]` occurrence table by
  amplicon ID, overwriting the `identity`, `taxonomy`, and
  `references` columns. Rows are neither added nor removed. The
  assignments TSV it reads may carry a header row (the stampa
  `[S49]` `_taxonomy_stampa.tsv` does; the sintax
  `_assignments_sintax.tsv` intermediate and the empty-input
  fallback do not): a leading row whose first field is the literal
  `amplicon` column name is skipped, so the join is correct whether
  or not a header is present. The
  output filename is `<project>_<N>_samples_table_assigned.tsv`
  — a sibling of Part B's `<project>_<N>_samples_table.tsv`,
  **not** an overwrite (see [`DECISIONS.md`](DECISIONS.md) D04
  sub-question 2). When Part B and Part C run end-to-end, both
  files are present in the same results folder.
  - **Pass when:** running on a Part B-shaped occurrence table
    plus a stampa-shaped assignments TSV publishes
    `<basename>_table_assigned.tsv`; rows whose amplicon ID
    appears in the assignments get their identity / taxonomy /
    references columns overwritten; rows missing from the
    assignments pass through unchanged; the input `_table.tsv`
    is not modified.
- `[S52]` **Retired.** Was: Part A's `filter_and_convert_to_fasta`
  normalised every `U`/`u` to `T`/`t` to defend the shadow Part B
  U-masking scheme. With the A-padding redesign ([S04] / [S63]),
  no character is reserved any more, so the awk pre-pass was
  removed and the contract dropped. ID is preserved so prior git
  history stays readable.
- `[S53]` **Retired.** Was: Part B's `discover_fasta.py` rejected
  any input `.fas` whose sequence lines contained `U` or `u`. The
  check existed only to defend the U-masking round-trip; the
  A-padding redesign ([S04] / [S63]) leaves no reserved character
  to protect, so the check was removed. ID is preserved so prior
  git history stays readable.
- `[S54]` every fastq emitted by a `vsearch`-based process keeps
  the canonical 4-line layout (one record = `@header`, sequence,
  `+`, quality on four consecutive lines). vsearch produces this
  by default; the spec exists as a regression guard so that
  downstream consumers (further vsearch invocations, awk-based
  helpers, etc.) which assume the 4-line layout cannot drift
  unnoticed from the upstream contract. Every per-process nf-test
  of a vsearch fastq-emitting module (`merge_fastq_pairs`,
  `strip_reads`, `join_notmerged`) asserts the layout of each
  emitted fastq.
  - **Pass when:** every non-empty fastq output of those processes
    has a length divisible by 4 lines; lines (4k+1) start with
    `@`, lines (4k+3) start with `+`, and lines (4k+2) and (4k+4)
    have equal length for each record.
- `[S55]` every fasta emitted by a `vsearch`-based process keeps
  sequences on a single line (one record = `>header` + sequence on
  two consecutive lines). This is enforced upstream by passing
  `--fasta_width 0` to every fasta-emitting vsearch invocation and
  is part of the Part A→B interface contract. The downstream
  helpers `bin/cluster_cleaver.py` and
  `bin/build_filtered_contingency_table.py` parse fasta
  line-by-line and would silently truncate folded records, so the
  contract is asserted by the per-process nf-tests of every
  vsearch fasta-emitting module
  (`filter_and_convert_to_fasta`, `dereplicate_fasta`,
  `global_dereplication`).
  - **Pass when:** every non-empty fasta output of those processes
    has a length divisible by 2 lines, with even-indexed lines
    starting with `>` and odd-indexed lines containing no `>`
    character.
- `[S106]` every `vsearch` process that **reads** fastq accepts the
  full quality range representable under the active encoding by
  default, so high-accuracy reads — e.g. PacBio HiFi, whose per-base
  quality reaches the maximum — are not rejected. vsearch's default
  `--fastq_qmax` is 41 and aborts with a fatal error on any higher
  score; it also requires `offset + qmax <= 126` (126 is the last
  printable ASCII, `~`). The processes therefore derive
  `--fastq_qmax` as `126 - --fastq_ascii`, i.e. the highest score the
  encoding can represent: **93** at offset 33 (`fastq_encoding = 33`)
  and **62** at offset 64. The derived value is passed to
  `merge_fastq_pairs` (`--fastq_mergepairs`), `strip_reads`
  (`--fastx_filter`, both the R1 and R2 pass), and
  `filter_and_convert_to_fasta` (`--fastx_filter`). `merge_fastq_pairs`
  passes the same value as `--fastq_qmaxout` so the posterior quality
  of the merged region is written at up to that ceiling rather than
  clamped to the default 41. `join_notmerged` (`--fastq_join`) needs
  no flag — vsearch documents `--fastq_qmax` as ignored for that
  command — and `chimera_detection` (`--fastx_filter` on fasta
  representatives) reads no quality scores, so the option would be a
  no-op there.
  - **Pass when:** a fixture whose every base carries the maximum
    representable quality (`~`) is merged by `merge_fastq_pairs`
    without error and its merged quality is written at that ceiling
    (not clamped to Q41); the same reads are converted by
    `filter_and_convert_to_fasta` and stripped by `strip_reads`
    without vsearch aborting — under **both** `fastq_encoding = 33`
    (a hardcoded qmax of 93 would abort at offset 64 because
    `64 + 93 > 126`) and `fastq_encoding = 64`.
- `[S56]` shadow Part B is a separate workflow (`part_B_shadow`)
  called alongside `part_B` from both the end-to-end
  (`--project_name`) and the standalone (`--fasta_folder`) entry
  points, **only when `--recover_unmerged` is set** ([S78]); at the
  default (unset) it is never invoked, so no `_notmerged` Part B
  artefacts are produced even if `_notmerged.fas` files are present in
  a `--fasta_folder`. When enabled it runs the **same processes as
  regular Part B with no
  mask / restore wrapping around `global_clustering`** — the
  A-padding emitted by Part A ([S04]) is composed entirely of
  A/C/G/T, so swarm accepts the sequences as-is and the
  representatives carry the run of `A`s through every downstream
  step unchanged.
  The shadow basename is `<project_name>_<N>_samples_notmerged`,
  so every published artefact carries the `_notmerged` token (e.g.
  `<basename>_notmerged_table.tsv`,
  `<basename>_notmerged_1f.swarms`). The shadow path consumes the
  Part A shadow outputs (`<sampleId>_notmerged.{fas,qual,stats}`,
  see [S04]); when no shadow samples are present (no fastq pair
  failed to merge), the shadow workflow is not invoked. The run
  of `A`s in shadow representatives is artificial join padding,
  not biological sequence — see [S04] for the user-facing caveat.
  - **Pass when:** an end-to-end run on a paired-end fixture that
    produces at least one `_notmerged` sample publishes both the
    regular `<project>_<N>_samples_table.tsv` and the shadow
    `<project>_<M>_samples_notmerged_table.tsv` under
    `--results_folder`; the shadow representatives FASTA contains
    a run of `params.join_padding_length` `A`s on its sequence
    lines and no `N` or `U`.
- `[S57]` `nextflow run main.nf --help` (and the short form
  `--help` set via config) prints a usage block to stdout and
  exits cleanly without running any process. The usage block:
  documents the three entry points (end-to-end Part A→B[→C],
  Part B standalone via `--fasta_folder`, Part C standalone via
  `--occurrence_table`); lists every `params.*` flag grouped by
  the part that consumes it, with type and default; points to
  `README.md` and `SPECIFICATIONS.md` for the behaviour contract.
  The `--help` short-circuit fires before any parameter assertion
  so users can discover the interface without first supplying a
  mode-specific required flag.
  - **Pass when:** running the workflow with only `--help` set
    succeeds, stdout contains a `Usage:` line, references each of
    `--fastq_folder`, `--fasta_folder`, and `--occurrence_table`,
    and no process is executed.
- `[S58]` `params.publish_mode` is the publishDir mode used by every
  `publishDir` directive in the workflow. Default `'link'` (hard
  link) preserves the prior fast-publish behaviour on
  same-filesystem layouts; users whose results folder is on a
  different filesystem than the work directory must override it to
  `'copy'` (or another non-link mode) because hard links cannot
  cross devices. Allowed values: `'copy'`, `'copyNoFollow'`,
  `'link'`, `'move'`, `'rellink'`, `'symlink'` (the set Nextflow
  accepts on `publishDir mode:`). The value is validated at workflow
  startup; an unknown mode aborts with a clear message before any
  process runs. Note: `'symlink'` and `'rellink'` produce links that
  point into the work directory and break when `cleanup = true`
  removes work files; `'link'`, `'copy'`, `'copyNoFollow'`, and
  `'move'` are safe under cleanup.
  - **Pass when:** an end-to-end Part A run with
    `--publish_mode symlink` publishes a `*_merging.log` that is a
    symbolic link (per `Files.isSymbolicLink`); an invalid value
    such as `--publish_mode bogus` aborts the workflow with an
    error mentioning `publish_mode`.
- `[S59]` `<outdir>/occurrence_table/` ([S71]) holds a closed
  whitelist of **data** artefacts per regular / shadow run (the step
  logs moved to `<outdir>/logs/part_b/` and `<outdir>/logs/part_c/` —
  D15/D16):
    1. exactly one occurrence table: `<basename>_table.tsv` (from
       [S46]);
    2. exactly one FASTA: the post-mumu
       `<basename>_table.fas` emitted by `extract_mumu_fasta`
       (column 10 of the final table, header
       `>amplicon;size=total;`).
  The directory contains **no** `*.log`: the six canonical `[S45]` step
  logs are published to `<outdir>/logs/part_b/` and Part C's
  `<basename>_taxonomy.log` to `<outdir>/logs/part_c/` instead (D16).
  Both Part C paths publish that log: the sintax path from vsearch's
  `--log`, and the stampa path by gathering its per-chunk `vsearch.log`
  slices into the single `<basename>_taxonomy.log` (D16).
  When Part C runs, these additional data artefacts are published into
  the same directory (they are not Part B intermediates and are not
  covered by the Part B whitelist above):
    - `<basename>_table_assigned.tsv` — the joined occurrence table
      ([S51]); on the shadow path its sibling
      `<basename>_notmerged_table_assigned.tsv` ([S50]);
    - the standalone taxonomy table, named per assignment method
      `<basename>_taxonomy_<method>.tsv`: `_taxonomy_stampa.tsv` on the
      stampa path ([S49]) and `_taxonomy_sintax.tsv` on the regular
      sintax path ([S61], vsearch's verbatim 4-column `--tabbedout`
      with a header). The standalone sintax table is published **only**
      when the user explicitly selects `--taxonomy_method=sintax`;
      the shadow path ([S50]) never publishes it (its `_notmerged`
      basename is dropped by the process's `saveAs` closure);
    - `<basename>_taxonomy_stampa_majority.tsv` ([S66]) when
      `--majority_assignment` is set.
  Part C's per-amplicon **intermediates** stay in the work directory
  and are not published: the stampa path's chunk TSVs
  (`stampa_chunk.tsv`) and the sintax path's canonical 5-column join
  table `<basename>_assignments_sintax.tsv` (the file
  `update_occurrence_table` consumes — distinct from the published
  4-column `_taxonomy_sintax.tsv`). The run's
  `software_versions.yml`
  ([S68]) lives in the sibling `<outdir>/pipeline_info/`, **not** in
  `occurrence_table/`. All other Part B intermediates (`.qual`,
  `.distr`, `_per_sample_OTUs.stats`, the global dereplicated `.fas`,
  the swarm artefacts `.swarms` / `.stats` / `.struct` /
  `_representatives.fas` and their cleaver siblings,
  `.uchime` / `.uchime2`, `_representatives.results` /
  `.results2`, the pre-mumu `.nosubstringOTUs.fas`) stay in the
  Nextflow work directory and **must not** appear under
  `<outdir>/occurrence_table/`. The shadow Part B path follows the
  same contract with the `_notmerged` token in the basename.
  - **Pass when:** an end-to-end Part B-only run (no Part C) produces
    an `<outdir>/occurrence_table/` whose file set matches
    `{<basename>_table.tsv, <basename>_table.fas}` (plus the same set
    with `_notmerged` if the shadow path fired) and contains no `*.log`;
    the six `[S45]` step logs appear under
    `<outdir>/logs/part_b/` instead; none of the blacklisted
    intermediate filenames appear; `software_versions.yml` is in
    `<outdir>/pipeline_info/`, not `occurrence_table/`.
- `[S60]` Path-typed parameters are normalised at workflow startup so
  shell-style `~` prefixes that the shell did not expand (quoted on
  the CLI, or read from a `-params-file`) still resolve to the
  intended location instead of being joined to the launch directory.
  The parameters covered are `--reference_dataset`,
  `--reference_dataset_sintax` ([S64]), `--occurrence_table`,
  `--fastq_folder`, `--fasta_folder`, and `--results_folder`. A
  bare `~` or a leading `~/` is replaced with
  the current user's home directory (`$HOME`); a leading `~user/` is
  replaced with that user's home directory when it can be resolved
  (best-effort lookup via `getent passwd` then `/etc/passwd`).
  Unrecognised `~user` prefixes and paths without a leading `~`
  (absolute or relative) pass through unchanged. List- and
  comma-separated values (`--fastq_folder dir1,dir2`) are normalised
  element-wise.
  - **Pass when:** `normalize_path("~/foo.fas")` returns
    `${HOME}/foo.fas`; an end-to-end run with
    `--reference_dataset "~/<file>"` (quoted, so the shell does
    *not* expand) succeeds and produces the expected output, with
    no work-dir symlink containing a literal `~`.
- `[S61]` `params.taxonomy_method` controls the **regular** Part
  C assignment method only ([S49] / [S50]'s sibling sintax run on
  the regular table). Accepted values: `'stampa'` (default, the
  [S49] scatter-gather, consumes `--reference_dataset`) and
  `'sintax'` (use `vsearch --sintax` on the regular table with
  the same reshape rules as [S50]; consumes
  `--reference_dataset_sintax` — see [S64]). On the regular sintax
  path (and **only** there, never the shadow path) the standalone
  taxonomy table `<basename>_taxonomy_sintax.tsv` is published — the
  sintax counterpart of stampa's `<basename>_taxonomy_stampa.tsv`
  ([S49]). It carries vsearch's verbatim 4-column `--tabbedout`
  layout under a header row
  (`query\ttaxonomy\tstrand\tcutoff_taxonomy`, the manpage column
  order), not the canonical 5-column join format — that format lives
  in the unpublished `<basename>_assignments_sintax.tsv` intermediate
  ([S59]). The shadow path always uses sintax regardless of this flag
  and reads `--reference_dataset_sintax`. The value is validated at
  workflow startup; an unknown method aborts with a clear message
  before any process runs.
  - **Pass when:** `--taxonomy_method bogus` aborts the workflow
    with an error naming `taxonomy_method`; `--taxonomy_method
    sintax` with `--reference_dataset_sintax` set runs the
    regular path through `assign_taxonomy_sintax`, publishes
    `<basename>_table_assigned.tsv`, and also publishes a
    `<basename>_taxonomy_sintax.tsv` whose first line is the
    `query\ttaxonomy\tstrand\tcutoff_taxonomy` header;
    `--taxonomy_method sintax` without
    `--reference_dataset_sintax` aborts with a clear message.
- `[S62]` In Part C **standalone mode** (`--occurrence_table
  /path/to/<basename>_table.tsv`), the workflow looks for a
  shadow occurrence-table sibling at
  `<dirname>/<basename>_notmerged_table.tsv`. If `--recover_unmerged`
  is set ([S78]) **and** that file exists **and**
  `--reference_dataset_sintax` is set ([S64]),
  `part_C_shadow` ([S50]) runs on the sibling alongside the
  regular `part_C`, publishing both
  `<basename>_table_assigned.tsv` and
  `<basename>_notmerged_table_assigned.tsv` under
  `--results_folder`. If `--recover_unmerged` is unset (the default),
  the sibling does **not** exist, or `--reference_dataset_sintax` is
  unset, the shadow workflow is not invoked and only the regular
  `_table_assigned.tsv` is published — no error. The opt-in is
  `--recover_unmerged`; beyond it, the presence of the sibling file
  **and** the sintax reference are the further preconditions.

  The sibling toggle is expressed as **channel logic, not a
  parse-time disk probe** (D07): `part_C_shadow` is wired
  unconditionally whenever `--reference_dataset_sintax` is set (a
  config check), and the sibling is fed through a staged
  `Channel.fromPath(..., checkIfExists: false).filter { it.exists() }`
  that empties at runtime when the file is absent, so the branch
  self-suppresses. The DAG shape therefore depends only on the
  parameter, never on head-node disk state at parse time — keeping
  `-resume` correct and the workflow remote-`results_folder` safe.
  - **Pass when:** running Part C standalone with
    `<basename>_table.tsv`, `<basename>_notmerged_table.tsv`, and
    `--reference_dataset_sintax` set publishes both
    `_table_assigned.tsv` files; running with the sibling file
    on disk but no `--reference_dataset_sintax` publishes only
    the regular assigned table and exits cleanly; running with
    only `<basename>_table.tsv` publishes only the regular
    assigned table and exits cleanly.
- `[S63]` `params.join_padding_length` (positive integer, default
  `8`) sets the length of the run of `A`s that
  `vsearch --fastq_join` inserts between R1 and R2 in the shadow
  Part A path ([S04]). The same length is used for the matching
  `--join_padgapq` quality string (Phred 40, `I`). The parameter
  is the workflow's single knob over the join padding: increase
  it to make the artificial join site more visible to a human
  reader at the cost of slightly larger fastas; decrease it to
  shave bytes when downstream tooling already handles the
  metadata. Non-positive values, non-integers, and other invalid
  inputs abort the workflow at startup before any process runs.
  - **Pass when:** running Part A on a paired-end fixture whose
    reads cannot overlap, with `--join_padding_length 16`,
    publishes a shadow `.fas` whose sequence lines contain a run
    of 16 `A`s (and no run of 8 `A`s introduced by a stale
    default); running without the flag uses 8; running with
    `--join_padding_length 0` aborts the workflow before
    `merge_fastq_pairs` is scheduled, with stderr naming the
    parameter.
- `[S64]` `params.reference_dataset_sintax` (optional, default
  `null`) is a path to a **sintax-formatted** reference fasta
  (header `>id;tax=d:Dom,p:Phy,...;` — see vsearch's
  `--sintax` documentation). The same `.fasta`/`.fas`/`.gz`/`.bz2`
  shapes accepted by `[S47]`'s `--reference_dataset` apply. The
  flag is consumed by:
    1. the shadow Part C path ([S50]), which always uses
       `vsearch --sintax`; and
    2. the regular Part C path ([S49] / [S61]) when
       `--taxonomy_method=sintax`, where it replaces
       `--reference_dataset` (stampa- and sintax-formatted
       references differ in incompatible ways; a single file
       cannot satisfy both methods in production use).
  When `--recover_unmerged` is **unset** (the default, [S78]) or
  `--reference_dataset_sintax` is **unset**, the shadow Part
  C workflow is silently skipped from every entry point (no
  error, no published shadow artefact) and the regular path
  must run through stampa. Path-typed normalisation ([S60])
  applies — a leading `~` is expanded at workflow startup. This
  parameter is a stop-gap exposed to give users a way out of the
  one-format-fits-all bind until a better solution is identified
  (e.g. an on-the-fly format converter or a multi-format reference
  loader); the existence of two separate flags is the visible
  workaround.
  - **Pass when:** an end-to-end run that produces at least one
    `_notmerged` sample publishes
    `<basename>_notmerged_table_assigned.tsv` only when
    `--reference_dataset_sintax` is set; the same run without the
    flag publishes only the regular `<basename>_table_assigned.tsv`
    and no shadow artefact; running Part C with
    `--taxonomy_method=sintax` and no `--reference_dataset_sintax`
    aborts at startup with stderr naming
    `reference_dataset_sintax`.

- `[S65]` `params.hash_function` (optional, default `'sha1'`) selects
  the hash vsearch uses to rename amplicons. Accepted values are
  `'sha1'` and `'md5'`; any other value aborts at startup with stderr
  naming `hash_function`. The flag is consumed only by
  `filter_and_convert_to_fasta` ([S01]), which passes the matching
  `--relabel_sha1` / `--relabel_md5` to vsearch — every amplicon name
  downstream is therefore a SHA1 (40 hex chars) or MD5 (32 hex chars)
  digest. No pipeline step assumes a fixed hash-string length:
    1. the `uniq --check-chars=<width>` dedup in
       `extract_expected_error_values` ([S01]) and
       `build_expected_error_file` ([S28]) derives `<width>` from
       `hash_function` (40 for sha1, 32 for md5) rather than
       hard-coding it;
    2. `extract_ee.awk` ([S01]) splits the fasta header on `>;=`
       delimiters and never measures the name; and
    3. the Python occurrence-table builders key amplicons on the
       hash **string** (dictionary lookup), with no length or
       prefix-bucketing assumption.
  Switching `hash_function` changes only the amplicon names; the set
  of clusters, abundances, and the occurrence-table shape are
  otherwise identical.
  - **Pass when:** a run with `--hash_function md5` produces
    per-sample fastas whose headers carry 32-hex-character names and a
    `.qual` file with one row per unique 32-char name (lowest ee per
    name preserved); the default (`sha1`) run keeps 40-char names; and
    a run with an unsupported `--hash_function` value aborts at startup
    with stderr naming `hash_function`.

- `[S66]` `params.majority_assignment` (boolean, default `false`) is an
  **opt-in** step that runs at the very end of the **regular** Part C
  path ([S49] / [S51]), after `update_occurrence_table`. When set, the
  process `compute_majority_assignment` runs
  `bin/majority_assignment.py` on the regular assigned table
  (`<basename>_table_assigned.tsv`, [S51]) plus the stampa-formatted
  `--reference_dataset` ([S47]), and publishes an **independent** table
  `<basename>_taxonomy_stampa_majority.tsv` under `--results_folder`
  (a sibling of [S51]'s output — the assigned table is not modified).
  The new table has three columns,
  `OTU\tamplicon\ttaxonomy_majority`, one row per OTU. For each OTU the
  step takes every reference accession in the `references` column,
  looks up its `|`-separated lineage in the reference, and reports —
  per taxonomic rank — the most frequent name annotated with its
  support as `name (hits/total)` (`total` = number of accessions,
  `hits` = how many carry that name at the rank); ranks are
  `|`-joined. Lineages of unequal depth are padded with `NA`
  (`itertools.zip_longest`); ties are broken by first-seen order
  (`Counter.most_common`). A `No_hit` row is passed through verbatim as
  `No_hit`. The reference is read transparently whether plain, gzip- or
  bzip2-compressed (detected by magic bytes). The step is wired into
  every entry point where the regular `part_C` runs (end-to-end
  A→B→C, Part B standalone with a reference, and Part C standalone)
  and **never** on the shadow path ([S50]). Because the majority vote
  needs the per-hit accessions that only the stampa method populates
  (sintax leaves `references` at the `NA` placeholder, [S50]), the flag
  requires `--taxonomy_method=stampa`: setting `--majority_assignment`
  together with `--taxonomy_method=sintax` aborts at startup with a
  message naming `majority_assignment`.
  - **Pass when:** golden-file characterization tests for
    `bin/majority_assignment.py` reproduce the byte-exact output of the
    legacy `majority_assignment.py` on a fixture covering: (a) a single
    reference accession, (b) several accessions that agree at every
    rank, (c) accessions that disagree at a rank (most-frequent wins,
    ties broken by first-seen order), (d) lineages of unequal depth
    (`NA` fill), and (e) a `No_hit` row; the same output is produced
    from a plain, a `.gz`, and a `.bz2` reference. A Part C run with
    `--majority_assignment true` (stampa) publishes
    `<basename>_taxonomy_stampa_majority.tsv` alongside
    `<basename>_table_assigned.tsv`; the same run without the flag
    publishes no majority table; running with
    `--majority_assignment true --taxonomy_method sintax` aborts at
    startup with stderr naming `majority_assignment`.


## Fetch — download from ENA/SRA

- `[S97]` setting `--accession` selects a standalone fetch mode that
  downloads raw fastq files for one or more accessions: a single
  accession or a comma-separated list (`--accession PRJEB89924,SRP012345`);
  a `nextflow.config` may instead supply a Groovy list
  (`accession = ['PRJEB89924', 'SRP...']`), mirroring `--fastq_folder`
  (`[S10]`). Dispatch is by param presence inside the entry workflow —
  the same mechanism as the other modes (`[S02]`), because Nextflow's
  strict config parser (25.10+) drops the `-entry` option. `--accession`
  is the sixth mutually-exclusive input-mode selector of `[S02]`; the
  fetch mode runs no Part A / B / C.
  - **Pass when:** `nextflow run main.nf --accession A,B` fans out to
    the fetch stages once per accession and invokes no Part A / B / C
    process; combining `--accession` with any other input-mode selector
    aborts at startup (`[S02]`).
- `[S98]` `--accession` accepts only bioproject and study accessions,
  validated at startup before any network access:
    - bioproject — `^PRJ(E|D|N)[A-Z][0-9]+$` (PRJEB / PRJNA / PRJDB)
    - study — `^(E|D|S)RP[0-9]{6,}$` (ERP / DRP / SRP)
  A value matching neither aborts at startup with a message naming the
  offending accession and the two accepted forms. Validation is a
  pure param check (no network), so it fires regardless of
  connectivity.
  - **Pass when:** `--accession PRJEB89924` and `--accession SRP012345`
    are accepted; `--accession SRR123` (a run accession) and
    `--accession bogus` each abort at startup naming the value and the
    accepted forms.
- `[S99]` fetch resolves each accession to its constituent run
  accessions in a **resolve** stage (one task per `--accession`)
  before any fastq download, so the download stage can fan out one
  task per run.
  - **Pass when:** a stubbed resolve stage emitting a fixture run list
    for one accession produces one downstream download task per run.
- `[S100]` fastq files for each accession are published under a
  subfolder named after that accession
  (`<outdir>/<accession>/<run>_{1,2}.fastq.gz`), so runs from distinct
  accessions in a comma-separated list never collide in a shared
  directory.
  - **Pass when:** `--accession A,B` produces `<outdir>/A/…` and
    `<outdir>/B/…`; no run file is written directly under `<outdir>`.
- `[S101]` the **download** stage runs one task per run accession via
  `fastq-dl`, pinned exactly to `fastq-dl=4.0.1` (`[S69]`) as a
  **per-process** `conda` directive — NOT in `environment.yml` — so
  native-conda users who never run `fetch` never resolve it and Wave
  builds a fetch-only image on demand (D10), with `--provider ena` as
  the default provider and `--max-attempts 5` (raising fastq-dl's
  default of 3) so a transient ENA / network hiccup is retried more
  times before the run-task fails. The per-run granularity gives the
  failure contract: a run that fails to download fails only its own
  task; runs that succeeded are published and cached, so a re-run with
  `-resume` retries only the failed runs.
  - **Pass when:** the download process carries a `conda` directive
    pinning `fastq-dl=4.0.1` and passes `--provider ena` and
    `--max-attempts 5`; `environment.yml` does not mention fastq-dl; a
    full Part A/B/C run triggers no fastq-dl resolution; with a stub in
    which one run-task fails, the sibling run-tasks still publish their
    fixtures.


## Reproducibility

- `[S68]` every run records the versions of the external tools it
  relies on (`vsearch`, `swarm`, `cutadapt`, `mumu`) plus the Python
  interpreter into `software_versions.yml`, published under
  `<outdir>/pipeline_info/` ([S71]). The process
  `dump_software_versions` queries each tool's `--version` on the
  active environment (PATH / conda / module / container) and pipes the
  raw output through `bin/collect_versions.py`, which extracts the
  first `MAJOR.MINOR[.PATCH]` token per tool. A tool missing from the
  environment is recorded as `n/a` (not dropped) so the gap is visible
  in the report rather than silently absent. Because `--outdir` always
  resolves (default `results`), the versions file is published on every
  entry point — including a Part A-only run.
  - **Pass when:** a run publishes
    `<outdir>/pipeline_info/software_versions.yml` naming `vsearch`,
    `swarm`, `cutadapt`, `mumu`, and `python`, each with a non-empty
    value (a real `MAJOR.MINOR[.PATCH]` token or `n/a`);
    `bin/collect_versions.py` unit tests pin the extraction for each
    tool's real `--version` output and the missing-tool `n/a` case.
- `[S69]` the conda environment (`environment.yml`) and the CI tool
  setup (`.github/workflows/test.yml`) pin every shared external tool
  to an **exact** version (`name=X.Y.Z`, never a floating `>=` / `<=`),
  and the two files agree. This keeps the byte-exact characterization
  tests (`[S22]`, `[S35]`, `[S39]`, `[S44]`) reproducible: a floating
  pin lets an upstream release silently change tool output and flip the
  golden files. `mumu` is now distributed on bioconda, so it is a
  pinned `bioconda::mumu` conda dependency alongside `vsearch` /
  `swarm` / `cutadapt`; the `conda` profile resolves it instead of a
  from-source build.
  - **Pass when:** repo-level checks assert that every tool dependency
    in `environment.yml` uses an exact `=` pin, that `environment.yml`
    declares a `bioconda::mumu` package, and that the `vsearch` /
    `swarm` / `cutadapt` / `mumu` pins in `environment.yml` and the CI
    workflow are identical.
- `[S77]` the release version is declared in exactly two places —
  `manifest.version` in `nextflow.config` and `version` in
  `CITATION.cff` — and the two must agree, so a release can never ship a
  workflow whose citation metadata reports a different version. (The
  README "Releasing" section is the human-facing instruction to bump
  both, plus `CITATION.cff`'s `date-released`; this is the automated
  guard behind it.)
  - **Pass when:** a repo-level check asserts that `manifest.version`
    (`nextflow.config`) and `version` (`CITATION.cff`) are identical.
- `[S80]` repository metadata stays internally consistent so the
  workflow never ships dead documentation links or an undocumented
  release:
    - `manifest.homePage` (`nextflow.config`) and `url` (`CITATION.cff`)
      reference the **same** GitHub repository (`owner/repo`);
    - `manifest.defaultBranch` is declared, and `nextflow_schema.json`'s
      `$id` URL points at that branch (and the same `owner/repo`), so the
      published schema link resolves instead of 404-ing on a stale
      branch name;
    - a `CHANGELOG.md` exists and carries an entry for the current
      `manifest.version`.
  This is the metadata sibling of `[S77]`'s version-sync guard.
  - **Pass when:** a repo-level check asserts (a) the `manifest.homePage`
    and `CITATION.cff` `url` resolve to the same `owner/repo`; (b)
    `manifest.defaultBranch` is non-empty and the schema `$id` path
    contains `/<owner>/<repo>/<defaultBranch>/`; (c) `CHANGELOG.md`
    exists and contains the `manifest.version` string.


## Input samplesheet

- `[S70]` `--input <samplesheet.csv>` is the primary, validated way to
  declare samples: a header-keyed CSV in one of two profiles, selected
  by the columns present.
    - **fastq profile** (Part A) — columns `sample`, `fastq_1`,
      `fastq_2` (optional), `run` (optional). One row per sample; an
      empty `fastq_2` marks a single-end sample (merging skipped,
      [S21]). `run` is an optional batch label carried through as
      provenance only (no behavioural effect yet).
    - **fasta profile** (Part B standalone) — columns `sample`,
      `fasta`, `qual` (optional), `stats` (optional). `qual` / `stats`
      default to the `<fasta-dir>/<sample>.{qual,stats}` siblings when
      omitted. A `sample` ending in `notmerged` routes the row to the
      shadow Part B path ([S56]).
  Structural validation runs at startup, before any process, in
  `bin/parse_samplesheet.py`, which aborts **naming the offending row
  and column** when: the header is missing a required column or carries
  an unknown one; the profile cannot be inferred (neither `fastq_1` nor
  `fasta` present); a required cell is empty; two rows share a `sample`
  ([S13] / [S14]); or (fastq profile) a `sample` ends in the reserved
  `notmerged` suffix ([S23]). Path cells are `~`-expanded ([S60]);
  relative paths resolve against the launch directory (Nextflow's
  default). The existence of each listed file is enforced downstream by
  Nextflow's `file(..., checkIfExists: true)` when the row is staged (a
  missing file aborts the run, naming the path). Empty fastq / fasta
  inputs are allowed and travel through to the occurrence table ([S09]).
  Files referenced by the samplesheet are declared `file()` inputs, so
  Nextflow stages them and `-resume` sees a change to the samplesheet
  or any listed file.

  **Trust boundary — path cells are deliberately unrestricted.** The
  path columns accept any absolute or relative path, including `..`
  traversal outside the launch directory; this is by design, because a
  user legitimately points the workflow at data anywhere on their own
  filesystem. The consequence is that whoever writes the samplesheet
  (equivalently, whoever populates a scanned `--fastq_folder` /
  `--fasta_folder`) can make the run read and stage any file the
  launching user can read. The samplesheet and the input folders are
  therefore a **trust boundary**: run only samplesheets and folders you
  or a trusted party produced, exactly as you would only execute a
  script you trust. The structural guards ([S93] sample-ID charset,
  [S95] delimiter rejection, [S23] reserved suffix) constrain the
  *sample-ID* and *cell-integrity* surface — they intentionally do **not**
  sandbox the path cells' reach, which the operator owns.

  `--input` is mutually exclusive with the folder-scan inputs
  (`--fastq_folder` / `--fasta_folder`); setting `--input` together
  with either aborts at startup. The folder-scan inputs keep their
  [S10] multi-folder and [S11] / [S12] pattern-table semantics and
  auto-generate the same internal per-sample channel, but do not give
  `-resume` the same guarantee (a file added to a scanned folder is not
  seen unless the folder set changes); `--input` is recommended for
  reproducible runs. Dispatch ([S02]) gains `--input`, with the entry
  point (Part A vs Part B standalone) selected by the active profile.
  - **Pass when:** `bin/parse_samplesheet.py` unit tests cover, per
    profile — a valid sheet → normalized rows; a missing required
    column; an unknown column; an un-inferrable profile; a duplicate
    `sample` (error lists both rows); a reserved `notmerged` sample
    (fastq profile); an empty required cell; a single-end row (empty
    `fastq_2`); and sibling `qual` / `stats` defaulting (fasta
    profile). An end-to-end `--input` (fastq profile) run with
    `--project_name` / `--results_folder` set drives Part A → Part B
    and publishes the Part B occurrence table to `occurrence_table/`
    and the step logs to `logs/part_b/` ([S46] / [S45] /
    D16) under `--results_folder`, exactly like the equivalent
    `--fastq_folder` run; setting `--input` together with
    `--fastq_folder` (or `--fasta_folder`) aborts at startup. Part A's
    per-sample artefacts ([S19]) publish under `<outdir>/per_sample/`
    ([S71]) regardless of whether the samples came from `--input` or a
    folder scan.
- `[S95]` no samplesheet cell may contain a raw TAB, newline, or
  carriage return. `bin/parse_samplesheet.py` emits its normalized rows
  as TSV (columns joined by TAB, rows by newline) for the workflow to
  consume via `splitCsv`, so a cell carrying one of those characters —
  reachable through CSV quoting, e.g. a quoted `"a<TAB>b"` path cell —
  would shift columns or inject a phantom row into the normalized output.
  Such a cell aborts at startup, naming the offending row and column,
  before any process runs. The `sample` column is already constrained to
  the stricter `[S93]` charset; this rule closes the same hole in the
  path (`fastq_1` / `fastq_2` / `fasta` / `qual` / `stats`) and `run`
  columns.
  - **Pass when:** a samplesheet whose `fastq_1` cell contains an
    embedded TAB (or newline) exits non-zero naming the row and column;
    a sheet whose cells are delimiter-free is accepted.


## Output (`--outdir`)

- `[S71]` `--outdir <dir>` is the single root for every published
  artefact, in a fixed layout. **Data files and step logs are kept
  apart** (D15): data lands in the part subdirectory, logs land under a
  dedicated `<outdir>/logs/` tree organised by the producing pipeline
  stage (D16 — superseding D15's data-mirroring sub-layout).
    - `<outdir>/per_sample/` — Part A's per-sample **data** files
      ([S19]): `<sample>.{fas,qual,stats}`.
    - `<outdir>/occurrence_table/` — Part B **data** ([S46] / [S59]:
      `<basename>_table.tsv`, `<basename>_table.fas`) and Part C
      ([S51] / [S61] / [S66]: `<basename>_table_assigned.tsv`,
      `<basename>_taxonomy_<method>.tsv`,
      `<basename>_taxonomy_stampa_majority.tsv`), regular and
      `_notmerged` shadow siblings together.
    - `<outdir>/logs/part_a/per_sample/` — Part A's per-step **logs**
      ([S19]): `_merging` / `_trimming_forward` / `_trimming_reverse` /
      `_dereplicating` / `_clustering` (and their `_notmerged` shadow
      siblings, [S04]).
    - `<outdir>/logs/part_a/<basename>_read_counts.tsv` — the Part A
      read-count summary ([S86]), a sibling of `per_sample/`.
    - `<outdir>/logs/part_b/` — Part B's six step logs ([S45]), regular
      and `_notmerged` shadow siblings together.
    - `<outdir>/logs/part_c/` — Part C's `<basename>_taxonomy.log`
      (sintax `--log`, or the stampa path's gathered per-chunk
      `vsearch.log` slices), regular and `_notmerged` shadow siblings
      together.
    - `<outdir>/pipeline_info/` — `software_versions.yml` ([S68]); a
      sibling of `occurrence_table/`, **not** inside it. Nextflow's own
      execution reports (timeline / report / trace / dag) also live
      here; they are not step logs and live outside the `logs/` tree.
  Only stages that run produce a `logs/part_*` directory.
  Inputs are never written to. `publishDir` materialises the tree on
  first publish (D08). Two sibling helpers resolve each process's
  target — `publish_dir('<subdir>')` → `<outdir>/<subdir>` for data and
  `log_dir('<subdir>')` → `<outdir>/logs/<subdir>` for logs — so
  routing lives in one place rather than scattered across the modules.

  `--outdir` defaults to `results`. `--results_folder` is a
  **deprecated alias**: when `--outdir` is unset and `--results_folder`
  is set, the latter becomes the output root (with a deprecation
  warning); when both are set, `--outdir` wins
  (`effective_outdir(outdir, results_folder)`). `--fastq_folder`
  reverts to an **input-only** parameter — Part A no longer publishes
  into it. This is a breaking change to output locations (manifest
  version bumped; see the README migration note). Part B no longer
  aborts when `--results_folder` is unset ([S26] superseded) because
  `--outdir` always resolves to a value.

  `[S59]`'s closed whitelist re-anchors to
  `<outdir>/occurrence_table/` (the original eight entries — the
  `pipeline_info/` directory moves out to its own sibling under
  `<outdir>`).
  - **Pass when:** `effective_outdir` unit tests pin the precedence
    (`--outdir` > deprecated `--results_folder` > `results`) and `~`
    expansion. An end-to-end run (`--fastq_folder` or `--input`) with
    `--outdir <d>` publishes Part A data files under `<d>/per_sample/`
    and Part A step logs under `<d>/logs/part_a/per_sample/`; the Part B table
    and the Part C assigned table under `<d>/occurrence_table/`, with
    the Part B six step logs under `<d>/logs/part_b/` and the Part C
    `_taxonomy.log` under `<d>/logs/part_c/`; and `software_versions.yml`
    under `<d>/pipeline_info/`; nothing is written to the input folder;
    `<d>/occurrence_table/` matches the `[S59]` whitelist and contains
    no `*.log`; `--results_folder <d>` (no `--outdir`) routes to the
    same layout with a deprecation warning.

  Reconciliation: any earlier bullet or pass-when that names a publish
  location as `params.fastq_folder` / `--results_folder` now resolves
  to this `[S71]` layout — Part A data files under
  `<outdir>/per_sample/`, Part B / Part C tables under
  `<outdir>/occurrence_table/`, all step logs under the stage-organised
  `<outdir>/logs/part_{a,b,c}/` tree (D15/D16), `software_versions.yml`
  under `<outdir>/pipeline_info/` — with `--results_folder` the
  deprecated alias for `--outdir`.
- `[S72]` numeric parameters are range-validated at workflow startup,
  before any process is wired, so an out-of-range value aborts
  immediately with a message naming the parameter rather than surfacing
  as an obscure tool error mid-pipeline (the same fail-fast contract as
  `[S63]`'s `--join_padding_length`). Each accepted range is pinned to
  the vsearch option the parameter feeds:
    - `fastq_encoding` — `33` or `64` (vsearch `--fastq_ascii`)
    - `threads` — integer `1..256`
    - `percentage` — real in `(0, 1]` (cleaving threshold, `[S22]`)
    - `chimera_minsize` — integer `>= 1` (vsearch `--minsize`, `[S34]`)
    - `stripright` — integer `>= 0` (`0` disables the trim, vsearch
      `--fastq_stripright`, `[S24]`)
    - `iddef` — integer `0..4` (vsearch `--iddef`)
    - `stampa_chunk_size` — integer `>= 0` (`0` = no-split sentinel,
      `[S49]`)
    - `stampa_maxrejects` — integer `>= 0` (`0` = no-limit sentinel,
      `[S49]`)
    - `stampa_id` — real in `[0, 1]` (vsearch `--id`, `[S49]`)
    - `sintax_cutoff` — real in `[0, 1]` (vsearch `--sintax_cutoff`,
      `[S50]`)
    - `sintax_randseed` — integer `>= 0` (`0` = random-data-source
      sentinel, vsearch `--randseed`, `[S91]`)
  - **Pass when:** for each parameter above, a representative
    out-of-range value aborts the run before any process executes with
    stderr naming the offending parameter; the default configuration
    passes validation.
- `[S73]` when a reference dataset path is supplied, its first FASTA
  header is sniffed at workflow startup and must match the format the
  flag declares, so a swapped or mis-formatted reference aborts before
  any process runs instead of producing empty / wrong assignments:
    - `--reference_dataset` must be **stampa-formatted** — a
      space-separated lineage follows the id (`>id <lineage>`);
    - `--reference_dataset_sintax` must be **sintax-formatted** — the
      header carries a `tax=` annotation (`>id;tax=d:Dom,p:Phy,...;`).
  Only the first header line is read (a declared input, like the
  `--input` samplesheet header in `[S70]`); the rest is left for
  vsearch. Plain and gzip (`.gz`) references are sniffed; a bzip2
  (`.bz2`) reference is skipped with a warning (pure-Groovy startup has
  no bzip2 decompressor). A path that is unset, or set but not present
  on the launch filesystem, is not sniffed — presence is enforced by the
  mode-specific `[S47]` assert and `file()` staging. The check is
  format-only and runs for whichever flags are set, independent of the
  selected run mode.
  - **Pass when:** a sintax-formatted file passed to
    `--reference_dataset` aborts at startup naming `reference_dataset`; a
    stampa-formatted file passed to `--reference_dataset_sintax` aborts
    naming `reference_dataset_sintax`; a `.gz` reference is decompressed
    and sniffed the same way; correctly-formatted references pass.
- `[S74]` when primer trimming runs (`--no_trimming` false, `[S20]`),
  `--forward_primer` and `--reverse_primer` are validated at startup to
  be IUPAC nucleotide strings — the codes `A C G T U R Y S W K M B D H V
  N` plus `I` (inosine), upper- or lower-case, at least 3 characters —
  the same alphabet `reverse_complement.sh` understands. A primer that is
  too short or carries a non-IUPAC character aborts at startup naming the
  parameter, instead of producing a confusing cutadapt error mid-run.
  Primer values reach cutadapt through quoted shell interpolation so a
  stray character cannot reshape the command.
  - **Pass when:** the primer check accepts valid IUPAC primers
    (including lowercase, ambiguity codes, and `I`) and rejects
    too-short / non-IUPAC values; an end-to-end Part A run with a
    malformed `--forward_primer` aborts at startup naming
    `forward_primer`.
- `[S94]` the startup header sniffs — the `--input` samplesheet profile
  probe (`[S70]`) and the reference-dataset format check (`[S73]`) —
  read the first line through a **length-bounded** reader that pulls at
  most a fixed cap of characters from the (possibly gzip-decompressed)
  stream before the first newline. A header is a short line, so the cap
  never truncates a legitimate one; it exists so a crafted input whose
  first line is enormous (no newline for gigabytes) — or a gzip
  **decompression bomb** whose first "line" inflates without bound —
  cannot exhaust launcher memory at startup. Bounding the characters
  pulled from a gzip reader also bounds how much is decompressed. A
  single shared `read_bounded_line()` helper backs both call sites so
  they cannot diverge.
  - **Pass when:** `read_bounded_line()` returns at most `max_chars`
    characters, stops at (and drops) the first newline, strips a
    trailing carriage return, and returns `null` on empty input; the
    existing `[S70]` / `[S73]` sniffs keep passing through it.

- `[S75]` a site adapts the workflow to its cluster through a
  **`-c` site-config file**, never by editing `nextflow.config`. A
  documented template (`conf/site.config.example`) carries every knob a
  site is expected to override — slurm `slurm_queue` / `slurm_account` /
  `slurm_queue_size`, the resource ceiling (`max_cpus` / `max_memory` /
  `max_time`) and dataset/reference sizes (`[S07]`), the
  container/conda cache directories, and the environment-module names
  (`[S08]`) — to be copied, edited, and passed with
  `-c my-site.config`, which native Nextflow merges over the pipeline
  defaults and the active profiles. In addition,
  `params.slurm_clusterOptions` (default `null`) is a free-form
  passthrough appended verbatim to every sbatch submit line (QoS,
  constraints, reservations) — the slurm profile's `clusterOptions`
  closure combines it with `--account=<slurm_account>`, yielding `null`
  (directive omitted) when neither is set. This is the one
  cluster-specific submission knob that cannot be expressed through
  `--slurm_account` / `--slurm_queue` alone.
  - Config **resolution** is checked automatically by
    `tests/check-site-config.sh` (the shipped template parses and
    resolves; its overrides win over the slurm-profile defaults when
    layered with `-c`; the `clusterOptions` closure is wired to
    `params.slurm_clusterOptions`). Actual job submission stays a manual
    cluster smoke test (`[S07]`), so — per `[S00]` — the runtime
    evaluation of the `clusterOptions` closure is not unit-tested.
  - **Pass when:** `nextflow -c conf/site.config.example config -profile
    slurm` resolves and the template's `slurm_clusterOptions` /
    `singularity.cacheDir` overrides appear in the resolved config;
    `-profile slurm` (no `-c`) exposes `params.slurm_clusterOptions =
    null` and a `clusterOptions` closure that reads it.

- `[S76]` a bundled `demo` profile runs the whole pipeline
  (Part A → B → C) out of the box with **no fixture generation and no
  required flags**: `nextflow run main.nf -profile demo` (or
  `-profile demo,<engine>` to also validate a container setup in one
  command). The profile points `fastq_folder`, `forward_primer`,
  `reverse_primer`, `project_name`, and `reference_dataset` at the
  committed demo dataset under `assets/demo/` — a tiny synthetic
  mergeable paired-end sample plus a stampa-formatted reference,
  produced deterministically by `assets/demo/make_demo.sh`. Unlike the
  `tests/data/` fixtures (`.gitignore`d, regenerated per run), the demo
  files are committed as plain text so a fresh checkout runs with no
  setup. The reads are not biologically meaningful — the profile exists
  so a user, or a new cluster / container environment, can confirm the
  tools and the wiring work end to end before committing real data.
  - **Pass when:** `tests/check-demo-profile.sh` confirms the committed
    demo assets exist and that `nextflow config -profile demo` resolves
    `fastq_folder` / `forward_primer` / `reverse_primer` /
    `project_name` / `reference_dataset` to them (and
    `-profile demo,singularity` composes). The end-to-end run behaviour
    is covered by the Part A→B→C tests in `tests/main.nf.test`; an
    actual `-profile demo` run is the manual smoke test.


## Dependencies

In addition to the upstream-tested tools used by Part A
(`vsearch`, `cutadapt`, `swarm`), Part B's "mumu (ex-lulu)"
post-processing step requires
[`mumu >= 1.1.1`](https://github.com/frederic-mahe/mumu) on
`PATH`. mumu is a daughter-cluster filter that detects unreliable
OTUs by comparing per-sample co-occurrence with a parent OTU's
distribution.


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

The two sides of the brace token must differ (`{1,2}`, `{R1,R2}` —
not `{1,1}`): they are the R1/R2 discriminator, so equal sides would
make the derived R2 name identical to the R1 name and pair a file
with itself. A pattern with equal sides is rejected (`[S67]`).


## Execution robustness

- `[S81]` every piped process script runs under `set -euo pipefail`, so
  a failure in **any** stage of a pipe fails the task rather than being
  masked by the exit status of the last stage alone. Without it, a crash
  in the left-hand command of a pipe — e.g. the `vsearch --fastx_filter`
  feeding `vsearch --uchime_denovo` in `chimera_detection` (`[S34]`) — is
  silently swallowed and a truncated artefact flows downstream into the
  occurrence table. Two classes of pipe are deliberately exempted so
  `pipefail` does not turn a correct outcome into a spurious failure:
    - pipes whose non-zero exit is a **legitimate empty result** stay
      explicitly guarded. In particular the `grep "^>"` header scans in
      `fake_taxonomic_assignment`, `fake_taxonomic_assignment2`, and
      `build_distribution_file` exit non-zero on a header-less input
      (`grep` returns 1 when nothing matches), so they pre-create the
      output and trail the pipe with `|| true` — an empty input yields an
      empty output and a successful task (`[S09]` / `[S33]`).
    - minimum-size computations that close a pipe early
      (`... | sort -n | head -n 1` in `chimera_detection_post_cleave`,
      `[S37]`) are written SIGPIPE-safe (the producer is allowed to run
      to completion) so they do not abort under `pipefail` once the
      sorted output exceeds the OS pipe buffer on a large dataset.
  - **Pass when:** a process whose first pipe stage is forced to fail
    (`chimera_detection` on a non-FASTA `representatives` input) exits
    non-zero instead of publishing a truncated `.uchime`;
    `fake_taxonomic_assignment` on a header-less input still succeeds
    with an empty `.results`; and `chimera_detection_post_cleave` on a
    cleaved input large enough to overflow the pipe buffer completes
    without a SIGPIPE abort.
- `[S82]` `cleanup` defaults to `false` in `nextflow.config`, so a
  successful run **retains** its per-task `work/` directories. This keeps
  `-resume` working across separate invocations — the large-dataset
  workflow, where a follow-up run (adding Part C, re-running after a
  downstream tweak) should reuse cached tasks rather than recompute — and
  leaves a "succeeded but wrong" run inspectable for post-mortem
  debugging. The cost is that `work/` accumulates and is cleaned by hand
  (the README documents `rm -rf work/`); a site that prefers
  auto-reclaim on success sets `cleanup = true` in a `-c site.config`
  ([S75]). The `test` profile keeps `cleanup = false` explicitly so
  nf-test can read task outputs after a run.
  - **Pass when:** `nextflow config main.nf -flat` resolves
    `cleanup = false` with no profile, and also under `-profile test`.
- `[S83]` an air-gapped site — one whose compute (and possibly login)
  nodes cannot reach the network at task start — can run from a
  pre-built container image instead of building one on the fly with
  Seqera Wave ([S08]), via either of two paths, neither needing
  outbound network from the compute nodes:
    - **build once, run offline:** on a connected node, run the pipeline
      once under an engine profile ([S08]) so Wave builds the image and
      caches it in `singularity.cacheDir` / `apptainer.cacheDir` on
      shared scratch; later runs reuse the cached image with no network.
      This uses the existing [S08] profiles unchanged — documentation
      only.
    - **site-supplied image:** the site sets `process.container` to a
      pre-pulled image (a `.sif` path or registry URI) and enables the
      engine in its `-c site.config` ([S75]), composing with the
      executor profile (`-profile slurm`) rather than an engine profile
      ([S08]) — so Wave is never enabled and no build is attempted. A
      plain `-profile slurm` with no such `-c` sets no
      `process.container`, so the default behaviour is unchanged. The
      template in `conf/site.config.example` carries a commented block
      for this.
  Documented in the README "Running with containers" section. Offline
  execution itself is a manual cluster smoke test (as for [S08]); the
  config wiring is checked automatically.
  - **Pass when:** `nextflow config main.nf -flat -profile slurm` with a
    `-c` that sets `process.container` + the engine resolves that
    `process.container` and the engine, with no `wave.enabled = true`; a
    plain `nextflow config main.nf -flat -profile slurm` resolves no
    `process.container`.
- `[S84]` every run writes Nextflow's execution reports — `timeline`,
  `report`, `trace`, and `dag` — into `<outdir>/pipeline_info/`
  (alongside `software_versions.yml`, [S68]). The `trace` and `report`
  are the primary tool for tuning per-step resources on a new cluster:
  they record the requested cpus/memory against the observed
  `peak_rss` / `peak_vmem` / `realtime` / `%cpu` per task, which is how a
  site chooses `--dataset_size_gb` / `--reference_size_gb` ([S79]) or a
  per-step memory override ([S75]). Files use fixed names with
  `overwrite = true`, so a rerun / `-resume` refreshes them in place. The
  output root follows the [S71] precedence (`--outdir`, then the
  `--results_folder` alias, then `results`) when set via `--outdir` /
  `-params-file`. Because the report `file` paths are config values
  resolved at parse time — before profiles are applied, and they cannot
  be closures — a profile that sets `params.outdir` is not seen by them;
  the one shipped such profile, `demo` ([S76]), re-points the four report
  files to its `demo_results/` explicitly so its reports sit beside its
  other outputs.
  - **Pass when:** `nextflow config main.nf -flat` resolves `timeline` /
    `report` / `trace` / `dag` enabled with `file` paths under
    `<outdir>/pipeline_info/` (default `results/`), the `trace` carrying
    `peak_rss` / `peak_vmem`. Report *generation* itself (and following a
    run-time `--outdir`) is upstream Nextflow behaviour, exercised by the
    demo / cluster smoke runs, not unit-tested ([S00]).
- `[S85]` every process that invokes an external tool (or consumes one's
  output) declares a Nextflow `stub:` block that produces its declared
  outputs, so the whole pipeline runs under `nextflow run -stub-run` with
  none of vsearch / swarm / cutadapt / mumu installed. This validates the
  Part A→B→C channel topology (the `.join`s, the regular/shadow branch,
  the stampa scatter/gather) in seconds without tools or real data — the
  fast topology-CI / onboarding smoke check. The three
  input-discovery / samplesheet-validation processes (`discover_inputs`,
  `discover_part_b_fasta`, `validate_samplesheet`) are **exempt**: they
  are pure stdlib-Python glue that must run for real to bootstrap the
  per-sample channel from the filesystem and invoke no bioinformatics
  tool. The `stub:` block sits after the `script:` / `shell:` block (the
  order the strict config parser requires).
  - **Pass when:** a static gate confirms every non-exempt process
    declares a `stub:`; and `nextflow run -profile demo -stub-run`, with
    vsearch / swarm / cutadapt / mumu shadowed by stubs that exit
    non-zero (so any fall-through to a real script fails), completes and
    publishes the placeholder Part A/B/C artefacts (per-sample `.fas`,
    the occurrence table, the `_table_assigned.tsv`, and
    `software_versions.yml`).
- `[S86]` whenever Part A runs (fastq input — end-to-end **or** Part
  A-only; not the fasta-input Part B/C standalone paths), the pipeline
  publishes a per-sample read-count summary table to
  `<outdir>/logs/part_a/<basename>_read_counts.tsv`. The
  basename is the Part B construct `<project>_<N>_samples` when
  `--project_name` is set, otherwise `<N>_samples`, where `N` is the
  number of regular (non-`_notmerged`) samples. The table is a port of
  the legacy genotoul read-tracking summary: a tab-separated file whose
  header is `samples	reads	assembled	F	R	passing`, one row per
  regular sample, followed by a final `Total` row summing each numeric
  column. The five counts per sample are extracted from that sample's
  published Part A logs ([S19]):
    - `reads`     — the `Pairs` count from `<sample>_merging.log`
      (vsearch `--fastq_mergepairs`);
    - `assembled` — the `Merged` count from `<sample>_merging.log`;
    - `F`         — the `Reads with adapters` count from
      `<sample>_trimming_forward.log` (cutadapt forward pass);
    - `R`         — the `Reads with adapters` count from
      `<sample>_trimming_reverse.log` (cutadapt reverse pass);
    - `passing`   — the `Reads written (passing filters)` count from
      `<sample>_trimming_reverse.log`.
  A count whose source log is absent (single-end samples have no
  `_merging.log`; `--no_trimming` runs have no trimming logs) or whose
  source line is missing (cutadapt emits no summary for an empty sample)
  is recorded as `0` — every cell is numeric, never blank. Thousands
  separators in tool output are stripped. Shadow (`_notmerged`) samples
  are excluded; rows are sorted by sample ID for determinism.
  - **Pass when:** a unit test of `bin/build_read_counts.sh` on staged
    fixture logs reproduces the header, the per-sample rows (zeros for
    absent logs/lines), and the `Total` row; and an end-to-end fastq run
    publishes `<outdir>/logs/part_a/<basename>_read_counts.tsv`
    with the expected sample rows and column sums, while a fasta-input
    Part B standalone run publishes no such file.
