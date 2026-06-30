# Open questions / blocked spec items

Each entry below is a decision the workflow can't make for itself. As
long as a question is `open`, the `[Sxx]` IDs it lists cannot have a
green test — they show up as `blocked` in
[`tests/COVERAGE.md`](tests/COVERAGE.md).

When a decision lands, record the resolution inline (do **not** delete
the entry) and update the matching bullets in
[`SPECIFICATIONS.md`](SPECIFICATIONS.md) in the same commit.

| Status legend |
|---------------|
| `open`        | nobody has answered yet                       |
| `proposed`    | a draft answer exists, awaiting confirmation  |
| `resolved`    | answer is final; spec and tests reflect it    |


## D01 — N↔A round-trip for unmerged pairs

**Blocks:** `[S04]` (resolved), `[S05]` (still blocked on marker name)
**Status:** `resolved` — **superseded by `[S04]` / `[S63]` (A-padding redesign, 2026-05-23)**

> **Superseded.** The N→A mask described below was retired by the
> A-padding redesign in `[S04]` / `[S63]`. Shadow Part A now injects
> the padding directly via `vsearch --fastq_join --join_padgap`, so
> there is no mask, no rewrite, and no round-trip to track. The
> resolution body is kept for historical context only —
> `[S04]` / `[S63]` in [`SPECIFICATIONS.md`](SPECIFICATIONS.md) are
> the current contract.

**Resolution (2026-05-18):** the `N→A` rewrite is applied uniformly
to **every** `N` in the fasta sequence lines (no per-position
tracking, no sidecar). The masked fasta is consumed only by `swarm`
and is **not** published; the original (N-containing) fasta is the
one published as `<sampleId>_notmerged.fas`. There is no `A→N` round
trip back, since the `.fas` is preserved through dereplication and
serves as the authoritative shadow-pipeline artefact. The swarm
`.stats` SHA1 IDs match the `.fas` IDs because the SHA1 is computed
in `filter_and_convert_to_fasta` (before the mask step) and the mask
step rewrites only sequence lines, never headers.

A new `[S23]` reserves the `notmerged` suffix so a user-supplied
sample ID can never collide with a shadow-pipeline artefact.

Sub-questions, retained for the record:

- **Which positions get rewritten?** Every `N` (the rewrite is fed
  only to swarm; the published artefact keeps the original `N`s).
- **How is the round-trip tracked?** It isn't — the original is
  preserved upstream, so no tracking is needed.
- **Failure mode** if a sequence acquires unexpected `N`s mid-pipeline:
  vsearch's `--fastq_maxns 8` (shadow) / `--fastq_maxns 0` (regular)
  drops any read with too many `N`s before fasta conversion.


## D02 — Marker for unmerged-pair clusters

**Blocks:** `[S05]`
**Status:** `proposed` — working name is `sampleID_partial`

`SPECIFICATIONS.md` notes the placeholder name. A short marker
suffix (e.g. `+partial`) on the sample column in the occurrence
table would make spreadsheet filtering easier than a separate column.

To decide: keep `_partial`, switch to `+partial`, or add a boolean
`unmerged` column. Decision must be reflected in the occurrence table
schema in [`SPECIFICATIONS.md`](SPECIFICATIONS.md).


## D03 — Collision policy for same-named samples

**Blocks:** `[S13]`, `[S14]`
**Status:** `resolved`

**Resolution (2026-05-19):** **refuse**. Sample IDs **must** be
unique. When `bin/discover_fastq.py` (Part A) or
`bin/discover_fasta.py` (Part B) finds two or more input files that
derive to the same sample ID, the workflow exits non-zero **before
any process runs** and the error message lists every offending file
path so the user can rename / move / remove the duplicate.

Rationale: silent merge (option 1) is destructive and loses
provenance; suffix-disambiguation (option 3) hides the conflict and
makes downstream artefacts (logs, occurrence-table columns) hard to
trace back to the original fastq. Refusing keeps the input → sample
mapping injective and forces the user to resolve the conflict
explicitly.

Sub-decisions:

- **Where the check runs:** Part A — inside `bin/discover_fastq.py`
  (after `check_reserved_suffix`). Part B — inside
  `bin/discover_fasta.py`. Both expose a CLI that exits non-zero on
  duplicates and an importable `check_unique_sample_ids()` helper.
- **Error format:** stderr names every offending path, grouped by
  the colliding sample ID, e.g.

  ```
  error: duplicate sample IDs (each sample ID must be unique):
    A: /run1/A_1.fastq.gz, /run2/A_R1.fastq.gz
  ```

  This format is asserted by the unit tests, so future tweaks must
  stay backward compatible.


## D04 — Part C input mode and output policy

**Blocks:** `[S48]`
**Status:** `resolved` — sub-question 2 (2026-05-19), sub-question 1 (2026-06-30)

**Resolution of sub-question 2 (2026-05-19):** Part C publishes its
updated table as a **sibling file** named
`<basename>_table_assigned.tsv` rather than overwriting Part B's
`<basename>_table.tsv` in place. Rationale: keeping both artefacts in
the results folder preserves the unannotated Part B output for
inspection / debugging and makes the Part B → Part C handoff
non-destructive when the two parts run end-to-end. The `_assigned`
suffix signals that the `identity` / `taxonomy` / `references`
columns have been populated with real assignments. `[S51]` reflects
this output policy and is no longer blocked by D04.

**Resolution of sub-question 1 (2026-06-30):** a dedicated
`--representatives_fasta /path/to/representatives.fas` flag selects
fasta input. It is mutually exclusive with `--occurrence_table` and the
other input-mode selectors ([S02]); `--fasta_folder` is **not** reused
(it is a Part B directory input with a different cardinality and role).
When `--representatives_fasta` is set, Part C skips the occurrence-table
extraction and the join, running the assignment selected by
`--taxonomy_method` directly on the supplied fasta.

**Fasta-input output policy:** with no occurrence table to splice onto,
fasta-input Part C does **not** synthesise a table. Its sole deliverable
is the standalone `<basename>_taxonomy_<method>.tsv` that the assignment
step already publishes (the same artefact the table-input path emits
alongside `_table_assigned.tsv`): the 5-column, headered
`_taxonomy_stampa.tsv` for stampa, or the 4-column `_taxonomy_sintax.tsv`
for sintax. `<basename>` is derived from the fasta filename.
`--majority_assignment` is rejected in this mode (no table to compute a
per-OTU majority on). `[S48]` reflects this and is no longer blocked.

For the record, the original sub-question-1 flag candidates were a
dedicated mutually-exclusive flag (the option taken, named
`--representatives_fasta`) versus reusing `--fasta_folder`; and the
original fasta-output candidates were a legacy `*.results`-shaped TSV, a
synthesised minimal occurrence table, or a hard startup failure — all
superseded by the `_taxonomy_<method>.tsv` policy once the per-method
standalone tables landed.


## D05 — Two-table occurrence output (`--split-occurrence-table`)

**Blocks:** `[S15]` (sub-questions 1–5 resolved, unblocking Part B
implementation; sub-question 6 affects Part C output only and is
deferred)
**Status:** `partial` — sub-questions 1–5 resolved 2026-05-24;
sub-question 6 open, no current `[Sxx]` depends on it

`[S15]` promises a wide single-table form and a sparse two-table form
selected by `--split-occurrence-table`. The original spec wording
left the filenames, replacement semantics, and column projections
unspecified. The scoping pass on 2026-05-24 settled the following:

**Resolution of sub-question 1 — filename convention (2026-05-24):**
`<basename>_clusters.tsv` and `<basename>_occurrences.tsv`. Matches
the `<basename>_*` pattern every other Part B artefact follows
(`[S46]`'s `_table.tsv`, `[S45]`'s step logs, the shadow path's
`_notmerged` token), so a results folder for project `X` stays
grep-friendly per basename.

**Resolution of sub-question 2 — coexist vs replace (2026-05-24):**
**replace**. When the flag is set, only the split pair is published;
`_table.tsv` stays in the Nextflow work directory but does not reach
`--results_folder`. Rationale: the user picked a mode; publishing
both forms doubles the artefact count without giving downstream
consumers a stable contract about which to read.

**Resolution of sub-question 3 — long-format cluster key
(2026-05-24):** the `OTU` column (the renumbered 1..N integer). It
matches what users read off the wide-format table and is the natural
join key against `_clusters.tsv`. The SHA1 `amplicon` column lives
only in `_clusters.tsv` (recoverable via the join when needed).

**Resolution of sub-question 4 — zero-abundance rows in the
long-format (2026-05-24):** dropped. `[S09]` already states this
for empty samples ("not in the two-table mode long-format"); the
same rule applies to zero cells in non-empty samples for
consistency. The long form is sparse by design.

**Resolution of sub-question 5 — schema-section reconciliation
(2026-05-24):** the spec's "Occurrence table schema" section follows
the code (the 13 real metadata columns), not the other way round.
The original section listed four idealised columns (`cluster_id`,
`sequence`, `abundance_total`, `taxonomy`) that never matched the
bash port; the rewrite is honest about the legacy bookkeeping
columns (`cloud`, `length`, `abundance`, `chimera`, `spread`,
`quality`, `identity`, `references`) that downstream consumers
already rely on. Changing the published shape would be a separate
breaking-change decision.

**Sub-question 6 — Part C output shape in split mode (OPEN):**
should `update_occurrence_table` (`[S51]`) emit a split pair
(`<basename>_clusters_assigned.tsv` +
`<basename>_occurrences_assigned.tsv`) when the upstream Part B ran
in split mode, or always emit a single
`<basename>_table_assigned.tsv` regardless of upstream mode? The
working default (and what the current `[S15]` spec wording assumes)
is "Part C always emits a single assigned table", because:

  - Part C's contract is "splice taxonomy onto the per-cluster
    metadata columns" — those columns live in `_clusters.tsv`, so
    the long form is unaffected by Part C and would just be passed
    through unchanged.
  - Users who picked split mode for Part B may want the same layout
    for Part C, but that is a usage assumption, not a structural
    necessity.

This sub-question does not block `[S15]` (which scopes Part B's
output only). It can be resolved once Part B split lands and a real
consumer asks for split assigned output.


## D06 — Input model: validated samplesheet vs in-process folder-glob

**Blocks:** a new input-contract `[Sxx]` (next free is `[S70]`)
**Revises:** `[S10]`, `[S11]`, `[S12]`, `[S27]`; adds `[S70]`
**Status:** `resolved` — option 3 (2026-06-18)

**Resolution (2026-06-18):** option 3, implemented in two slices.
Slice 1: `bin/parse_samplesheet.py` (stdlib, pytest-tested) does
structural validation of the `--input` samplesheet in two profiles
(fastq → Part A, fasta → Part B), inferred from columns; single
`--input`, profile inferred from columns; launchDir-relative paths;
`run` column carried as provenance. Slice 2: wiring —
`validate_samplesheet` process (staged input); `part_A` and `part_b`
source their per-sample channels from `--input` *or* the folder scan
(`discover_*` paths left byte-identical as the fallback);
`samplesheet_profile()` drives the dispatch; mutual exclusion with
`--fastq_folder` / `--fasta_folder` asserted in `validate_params()`.
Full ci suite (47 tests) green.

Coupling found with **D09**: Part A publishes its per-sample artefacts
into `--fastq_folder`, which is null under `--input`. Three Part A
data-file modules (`extract_expected_error_values`, `dereplicate_fasta`,
`list_local_clusters`) published *without* an `enabled:` guard and so
hard-crashed ("publishDir target cannot be null") on a null
`fastq_folder` — a latent inconsistency, since `merge_fastq_pairs` /
`trim_primers` already had the guard. Added the guard to all three:
`--fastq_folder` runs publish as before; `--input` runs skip Part A
per-sample publishing until the unified `--outdir` lands (D09). The
end-to-end `--results_folder` outputs (Part B table, Part C,
`software_versions.yml`) are unaffected. This makes D09 the natural
next step.

The original analysis follows, for the record.

Discovery is performed *inside* processes that declare **no path
inputs** and glob `params.fastq_folder` / `params.fasta_folder` off
the launch filesystem (`modules/local/part_a/discover_inputs.nf`,
`modules/local/part_b/discover_part_b_fasta.nf` — both call the
`bin/discover_*.py` helpers on raw absolute paths). Consequences:

- Nextflow stages and hashes nothing, so `-resume` cannot detect that
  the input set changed, and provenance is opaque.
- On slurm/cloud executors the task may land on a node that cannot see
  the path unless the fastq/fasta dir is on shared storage (the README
  quietly requires exactly this). It would not work on an
  object-store-backed executor at all.

**Question:** adopt an nf-core-style validated samplesheet CSV as the
primary input, keep the folder-glob, or both?

1. **Samplesheet-only.** `--input samplesheet.csv`
   (`sampleId,forward,reverse,run`); discovery becomes channel logic
   (`splitCsv` / `Channel.fromPath`), so inputs are staged and hashed.
   Gains provenance, working `-resume`, cloud/object-store support, and
   an explicit per-run/batch column. Breaking UX change; the carefully
   specified pattern table (`[S11]` / `[S12]`) becomes legacy/optional.
2. **Keep folder-glob, fix staging only.** Move the globbing out of a
   process into a build-time Groovy function. Smaller change, but still
   reads the launch FS and still leaves `-resume` blind to directory
   contents; does not fix cloud / object-store.
3. **Both (recommended).** Samplesheet is the primary, staged,
   hashed input; folder discovery is retained as a convenience that
   *generates* a samplesheet-shaped channel up front, feeding the same
   downstream channels. Keeps the existing UX and the `[S11]` / `[S12]`
   pattern table as the single source of truth for auto-pairing, while
   gaining provenance and resume-correctness.

**Proposed resolution:** option 3.

Sub-questions for the human:
- Final samplesheet columns, and whether a `run`/`batch` column should
  drive per-run chimera/clustering grouping (it does not today).
- Whether duplicate-sample-ID refusal (`[S13]` / `[S14]` / D03) moves
  into samplesheet validation or stays in the discovery helpers.
- Whether `--fastq_folder` / `--fasta_folder` survive as sugar or are
  deprecated with a window.

Until D06 lands the discovery processes stay as-is; the current
`[S10]`–`[S12]` / `[S27]` tests remain valid and green.


## D07 — Shadow Part C standalone toggle: parse-time disk probe vs channel

**Revises:** `[S62]` (replaced the parse-time shadow probe in
`workflow part_c`, `main.nf`)
**Status:** `resolved` — option 2 (2026-06-18)

**Resolution (2026-06-18):** option 2. Replaced the
`shadow_table_path.exists()` parse-time branch with channel logic:
`part_C_shadow` is wired unconditionally whenever
`--reference_dataset_sintax` is set (a param/config check), and the
sibling is routed through
`Channel.fromPath(..., checkIfExists: false).filter { it.exists() }`,
which empties at runtime when the sibling is absent so the branch
self-suppresses. The DAG shape now depends only on the parameter, not
on head-node disk state at parse time. `[S62]` was reworded to match.

Note: `Channel.fromPath` on a literal **absent** path with
`checkIfExists: false` phantom-emits the path (it does *not* yield an
empty channel), so the runtime `.filter { it.exists() }` is required —
verified empirically before wiring. The existing `[S62]` tests (sibling
present + sintax ref; sintax ref set + no sibling; no sibling, no
sintax ref) are unchanged and stay green; the "sintax ref set + no
sibling" case is the regression guard for the new empty-channel path.

The original analysis follows, for the record.

`[S62]` specified that standalone Part C decides whether to run the
shadow path by testing `shadow_table_path.exists()` **at workflow-build
time** (`main.nf:48-55`), and branches the DAG topology on the result.
The DAG shape therefore depends on head-node disk state at parse time:
this defeats `-resume` reproducibility and breaks when the input table
/ results folder is remote.

**Question:** keep the parse-time `.exists()` probe, or express the
toggle as channel logic so topology is data-driven?

1. **Keep `[S62]` as-is.** Simplest; brittle as above.
2. **Channel logic (recommended).** Build the shadow branch
   unconditionally from `Channel.fromPath(sibling, checkIfExists:
   false)`; an empty channel self-suppresses the branch (no work, no
   output). Topology becomes static and the toggle becomes *data
   presence*, not a `java.io.File` call in the workflow body.
3. **Explicit flag.** Replace the implicit sibling-probe with an
   explicit `--shadow_table PATH` (or `--with-shadow`). Less magic, a
   small UX change; composes with option 2.

**Proposed resolution:** option 2, optionally with the explicit
override from option 3. `[S62]` reworded so the sibling is discovered
through a staged channel and the branch self-suppresses on an empty
channel — no `File.exists()` in the workflow body.

The current `[S62]` test stays green until the reword lands.


## D08 — Results-folder creation: build-time `mkdirs()` vs `publishDir`

**Revises:** `[S26]` (removed the `new File(...).mkdirs()` calls from
the three workflow-body sites in `main.nf`)
**Status:** `resolved` — option 2 (2026-06-18)

**Resolution (2026-06-18):** option 2. Dropped the three
`results_dir.mkdirs()` calls from the workflow body; `publishDir`
creates the results folder (and any missing parents) on first publish,
so the workflow performs no filesystem I/O at parse time. `[S26]` was
reworded to match. The `[S26]` integration test
(`tests/main.nf.test`, "Part B standalone …") is unchanged and stays
green — its `file(resultDir).exists()` assertion holds because the run
publishes artefacts into the folder, which `publishDir` materialises.
No new test was added: the change is behaviour-preserving by design (it
swaps the creation mechanism, not the observable contract), and the
existing test is the regression guard. The required-param assert
(`assert params.results_folder`) is untouched.

The original analysis follows, for the record.

`[S26]` mandated that the workflow "creates the folder (and any missing
parent directories) at startup." It was implemented as side-effecting
`java.io.File.mkdirs()` in the workflow body (three sites). Build-time
filesystem I/O runs on the head node only, so it breaks for a remote /
object-store `--results_folder` and is non-idiomatic — Nextflow's
`publishDir` already creates the directory tree itself.

**Question:** keep the explicit startup `mkdirs()`, or let `publishDir`
materialise the folder?

1. **Keep `[S26]` startup `mkdirs()`.** Works on local/shared FS;
   head-node only.
2. **Drop `mkdirs()`; rely on `publishDir` (recommended).** It creates
   the directory on first publish. Removes build-time FS I/O and works
   wherever the executor's publish layer works. Edge case: a run that
   publishes nothing leaves no folder — acceptable, since there is
   nothing to hold.
3. **Creation process.** Materialise the folder inside a tiny process
   so it runs in the task environment (consistent with remote
   executors). Heavier; rarely needed if option 2 suffices.

**Proposed resolution:** option 2; `[S26]` reworded to "`publishDir`
materialises the results folder; the workflow performs no filesystem
I/O at parse time." Coordinate with D09 (`--outdir`).

The current `[S26]` test stays green until the reword lands.


## D09 — Output routing: unified `--outdir` vs publishing into the input folder

**Blocks:** a new output-routing `[Sxx]` (next free after D06's
allocation, i.e. `[S71]`)
**Revises on resolution:** `[S19]` (Part A publishes into
`params.fastq_folder`); touches `[S45]`, `[S46]`, `[S58]`, `[S59]`, and
`[S68]`'s "the Part A-only path has no results folder yet" note
**Status:** `proposed` — **option 2 selected (2026-06-18); spec + code pending**

> **Note (2026-06-18):** the `[S68]` version-capture step publishes
> `pipeline_info/software_versions.yml` into `--results_folder`, which
> the `[S59]` closed whitelist must now permit. Allowing `pipeline_info/`
> in `[S59]` is a down-payment on this decision (the `pipeline_info`
> directory moves under `--outdir` when D09 lands).

Part A publishes its per-sample artefacts **back into**
`--fastq_folder` (`modules/local/part_a/trim_primers.nf:9` and
siblings, `[S19]`), while Part B/C publish to `--results_folder`. There
is no single `--outdir`. Publishing into the input directory mutates
the user's data, mixes inputs with derived files, complicates cleanup
(the README carries a manual "remove them from the fixture dir" step)
and re-runs, and leaves a Part A-only run with nowhere to put
`software_versions.yml` (`[S68]`).

**Question:** introduce a unified `--outdir` (nf-core convention:
inputs immutable, all outputs under `--outdir/<subdir>`) and stop
writing into `--fastq_folder`?

1. **Status quo.** Familiar to current users; the smells above persist.
2. **Unified `--outdir` (recommended).** Per-part subfolders, e.g.
   `<outdir>/part_a`, `<outdir>/part_b`, `<outdir>/pipeline_info`. Part
   A no longer writes into `--fastq_folder`; `--results_folder` becomes
   a back-compat alias. Gives `[S68]` a home on the Part A-only path,
   makes inputs immutable, and unifies cleanup. The `[S59]` whitelist
   anchors to `<outdir>/part_b` instead of `--results_folder`.
3. **Hybrid.** Keep Part A → `fastq_folder` by default, allow
   `--outdir` to override. A half-measure that keeps the
   input-mutation default.

**Proposed resolution:** option 2. This is a breaking change to output
locations, so it needs a version bump and a README/migration note.
Coordinate with D08 (folder creation).

Sub-questions for the human:
- The subfolder layout and names.
- Whether `--results_folder` / the `--fastq_folder` publish location
  survive as deprecated aliases for one release.
- Migration: which existing user scripts break, and whether a
  deprecation window is needed.

Until D09 lands, `[S19]`'s publish-into-`fastq_folder` contract stands
and its tests remain green.


## D10 — Container packaging strategy

**Blocks:** `[S08]` (container profiles)
**Status:** `resolved` — option 1 (Wave, no freeze) (2026-06-19)

**Question:** how should the `docker` / `podman` / `singularity` /
`apptainer` profiles obtain the tool images — given that `mumu` is now
on bioconda (so every dependency is conda-installable), the two adopting
labs run on slurm with **outbound network on the compute nodes**, and
the project values a single pinned source of truth (`[S69]`) without
nf-core membership?

1. **Seqera Wave, no freeze (chosen).** Wave builds the image on the
   fly from `environment.yml` and caches it; the same pinned spec that
   drives `-profile conda` drives the container build. No Dockerfile, no
   registry, no manual image bump. Rebuilds are deterministic from the
   exact conda pins, so tool versions never drift (only the image digest
   may differ). Cost: a soft dependency on the Wave service and outbound
   network at task start — acceptable here (compute nodes are online).
2. **Wave + freeze to a registry (ghcr.io).** Same declarative build,
   but Wave pushes each image into a repo the project owns, for
   permanent self-owned artefacts. Deferred: needs a registry + token;
   can be layered on later (flip `wave.freeze` + `wave.build.repository`)
   without changing the profiles or the workflow.
3. **Self-built Dockerfile on ghcr.io.** Full control and fully
   offline-capable, but the project maintains a Dockerfile + a
   build/push workflow and bumps the image tag with every
   `environment.yml` change. Rejected as more maintenance than the labs
   want.

**Resolution (2026-06-19):** option 1. Each engine profile enables its
engine + `wave.enabled` + `conda.enabled` + `process.conda =
environment.yml`, scoped inside the profile so a plain `nextflow run`
stays bare-PATH. Option 2 is the documented upgrade path if a
permanent image archive is later wanted. See `[S08]` in
[`SPECIFICATIONS.md`](SPECIFICATIONS.md).


## D11 — Shadow pipeline: always-on vs opt-in

**Blocks:** `[S04]`, `[S78]` (and touches `[S05]`, `[S50]`, `[S56]`,
`[S62]`, `[S64]`)
**Status:** `resolved` — opt-in, default off (2026-06-20)

**Question:** the shadow pipeline ([S04]) — recovering unmergeable read
pairs by A-padded joining — fired automatically whenever a pair failed
to merge, and (as seen on the `demo` profile, where 100 % of pairs
merge) `part_B_shadow` / `part_C_shadow` were invoked even with **zero**
not-merged reads, publishing empty `_notmerged` artefacts. The path is
explicitly experimental and its output carries a caveat (the run of
`A`s is artificial padding, not biological sequence). Should it stay
always-on, or become opt-in?

1. **Always-on (status quo).** No flag. But the experimental, caveated
   output is produced unasked, empty `_notmerged` files clutter the
   output even when nothing failed to merge, and the extra
   `part_B_shadow` / `part_C_shadow` work runs needlessly.
2. **Opt-in via `--recover_unmerged` (default `false`) (chosen).** The
   shadow path is off unless the user explicitly asks for it. Default
   runs are simpler (no `_notmerged` artefacts, no wasted shadow work),
   and enabling it is a deliberate, documented acknowledgement of the
   experimental caveat — "hard to misuse". Cost: a breaking change to
   default behaviour for anyone who relied on the shadow output, and a
   migration of the shadow tests to set the flag.

**Resolution (2026-06-20):** option 2. `params.recover_unmerged`
(default `false`) gates the whole shadow path: Part A drops the
not-merged reads instead of routing them into the A-padded join, and the
three entry points only invoke `part_B_shadow` / `part_C_shadow` when it
is set. The `[S23]` reserved-suffix guard stays always-on (independent
of the flag) so a user sample can never collide with shadow naming if
the flag is later enabled. See `[S78]` and the gated clauses on `[S04]`
/ `[S50]` / `[S56]` / `[S62]` / `[S64]` in
[`SPECIFICATIONS.md`](SPECIFICATIONS.md).


## D12 — Strict bash (`pipefail`) for piped process scripts

**Blocks:** `[S81]`
**Revises:** the piped `script:` / `shell:` blocks listed in the audit
below
**Status:** `resolved` — per-script `set -euo pipefail` (option 1),
implemented 2026-06-21 (`[S81]`)

**Resolution (2026-06-21):** option 1, implemented across all 15 piped
process scripts. Every piped `script:` / `shell:` block now carries
`set -euo pipefail`. The two pre-existing lines
(`chimera_detection_post_cleave`, `merge_substring_otus`) stay; nine
scripts (Bucket A) gained the line with no other change; the three
already-`|| true`-guarded scripts (`fake_taxonomic_assignment2`,
`build_distribution_file`, `search_for_terminal_gaps`) gained it on top
of their guards so the invariant is uniform. Bucket B
(`fake_taxonomic_assignment`) gained the `: >` pre-create + `|| true`
guard its sibling already had, so a header-less input still yields an
empty `.results`. Bucket C (`chimera_detection_post_cleave`) had its
`sort -n | head -n 1` minimum swapped to the SIGPIPE-safe
`sort -n | sed -n '1p'`. Covered by `[S81]` (three nf-test cases:
left-pipe failure propagates, header-less Bucket-B guard, large-cleaved
SIGPIPE regression). `-u` was kept (matches the two original scripts;
the audit confirmed no Bucket-A script references an unset var). The
global `process.shell` route (option 2) was rejected: a `#!/bin/bash`
shebang overrides `process.shell`, so it would silently skip the four
shebang-carrying scripts. The original analysis follows, for the record.

Nextflow does not enable `pipefail` for task scripts by default, and the
pipeline sets no global `process.shell`. So in a pipe like
`chimera_detection`'s

```
vsearch --fastx_filter ... --fastaout - | vsearch --uchime_denovo - ...
```

the task's exit status is the **right-hand** command's. If the left
`vsearch` dies (bad input, partial OOM), the task is still marked
successful and a silently truncated `.uchime` table flows into the
occurrence table — exactly the failure mode that goes unnoticed on a
large dataset. There are 15 piped process scripts; only **two**
(`chimera_detection_post_cleave`, `merge_substring_otus`) currently set
`set -euo pipefail`, so protection is inconsistent.

**Question:** make strict-bash uniform via a global `process.shell`, or
per-script `set -euo pipefail`?

1. **Per-script `set -euo pipefail` (chosen).** Add the line to each
   piped script (matching the two that already have it). Explicit,
   visible at the point of use, and — crucially — the only approach that
   actually works uniformly (see the shebang finding below).
2. **Global `process.shell = ['/bin/bash', '-e', '-o', 'pipefail']`.**
   One config line. **Rejected:** the audit found a `#!/bin/bash`
   shebang at the top of a script *overrides* `process.shell` entirely
   (verified empirically — a shebang script ignored a global pipefail),
   so this would silently skip the 4 scripts that carry a shebang
   (`trim_primers`, `list_local_clusters`, and the two already-protected
   ones). A global flag that covers only 11 of 15 scripts is worse than
   an explicit per-script line.

**Audit (2026-06-21).** Empirically established, in this environment:
`swarm`, all five piped `vsearch` subcommands, `cat`, `sort`, `awk` all
exit 0 on empty input; `grep` with no match exits 1; a `#!/bin/bash`
shebang overrides `process.shell`; a SIGPIPE in `cmd | sort -n | head`
aborts under `pipefail` **even inside a `$()` assignment**, but a
failing/SIGPIPE `$()` in a *printf-argument* position does **not** abort.
With those facts, adding `set -euo pipefail` is **not** uniformly "just
add a line" — the 15 scripts split into three buckets:

- **Bucket A — add the line, no other change (9 scripts):**
  `chimera_detection` (the headline win), `rebuild_post_mumu_table`
  (also catches a currently-masked python failure), `global_dereplication`,
  `find_similar_sequences`, `build_expected_error_file`,
  `extract_expected_error_values`, `dump_software_versions` (its
  `$(tool --version | head)` is in printf-argument position, verified to
  survive a missing tool → the `[S68]` `n/a` path is preserved),
  `trim_primers` and `list_local_clusters` (both already carry a
  `#!/bin/bash` shebang but **no** `set` line — add it there).
- **Bucket B — needs a guard *before* the line (1 script):**
  `fake_taxonomic_assignment` pipes `grep "^>" | sed`; `grep` exits 1 on
  a header-less representatives fasta, so the bare `set -euo pipefail`
  would turn a legitimate empty result into a task failure. Mirror the
  `: > file` pre-create + trailing `|| true` that its sibling
  `fake_taxonomic_assignment2` already uses, *then* the line is safe.
- **Bucket C — the line is already present and is itself the hazard
  (1 script):** `chimera_detection_post_cleave` already has
  `set -euo pipefail` **and** a latent SIGPIPE bug —
  `lowest="$(sed ... | sort -n | head -n 1)"`. On the tiny CI fixture the
  sorted output fits the 64 KB pipe buffer so `sort` finishes before
  `head` closes (exit 0); on a real study with enough distinct cleaved
  sizes, `head` closing early SIGPIPEs `sort` (exit 141) and `pipefail`
  aborts the task. This is a **pre-existing bug** (independent of this
  rollout) that bites exactly the large-dataset case the adopting labs
  care about. Fix with a SIGPIPE-safe minimum (`sort -n | sed -n 1p`, or
  a one-pass `awk` min) — *not* another `set` line.

The remaining 4 are already correct: `merge_substring_otus` (already
`set -euo pipefail`; its `head` reads a file, not a pipe);
`fake_taxonomic_assignment2`, `build_distribution_file`,
`search_for_terminal_gaps` (each already guards its pipe with `|| true`).
Minor note: `search_for_terminal_gaps`'s `|| true` is over-broad — it
masks a `vsearch` failure as well as the intended `grep "^H"` no-match;
narrowing it to the `grep` is an optional cleanup.

So the user's instinct — per-script `set -euo pipefail` on the
safe scripts — is the right direction *and* the right approach, but it
is a one-liner only for Bucket A (9 of 15). Bucket B needs a guard in
the same change; Bucket C is a separate SIGPIPE fix that should ship
alongside. `-u` is safe for every Bucket-A script: each uses only
Nextflow `!{...}` / `${...}` interpolation or locally-assigned bash vars
(no unset references), matching the two scripts that already carry
`set -euo pipefail`.

**Proposed resolution:** option 1, in one change covering all three
buckets. New `[Sxx]`: "every piped process script runs under
`set -euo pipefail`; a failure in any stage of a pipe fails the task,
and the documented empty-result pipes (`grep` no-match) stay guarded."
Tests: (a) a process whose first pipe stage is forced to fail exits
non-zero instead of publishing a truncated artefact; (b) the
`fake_taxonomic_assignment` Bucket-B guard keeps a header-less input
succeeding with an empty `.results`; (c) a `chimera_detection_post_cleave`
fixture large enough to exceed the pipe buffer no longer aborts (the
SIGPIPE regression guard). Touching the existing `script:` blocks and
their comments needs the usual authorization.

Until D12 lands, the two per-process `set` lines stand (Bucket C's bug
included) and the other piped scripts remain unprotected.


## D13 — Default run-directory retention (`cleanup`) vs `-resume`

**Blocks:** `[S82]`
**Revises:** `cleanup` in `nextflow.config`; the README cleanup/`-resume`
notes
**Status:** `resolved` — option 2 (default off), implemented 2026-06-21
(`[S82]`)

**Resolution (2026-06-21):** option 2. `cleanup` flipped from `true` to
`false` in `nextflow.config`, so a successful run retains its `work/`
directories and `-resume` works across separate invocations. The option-3
`--cleanup` param knob was **dropped**: a top-level `cleanup` directive
reads `params` at config-parse time (it does not defer like a process
closure), so a `--cleanup`/`-c` param override would not reliably reach
it, and `nextflow config` cannot validate a param-driven value anyway.
Sites that want auto-reclaim set `cleanup = true` directly in their
`-c site.config` ([S75]) — the established override mechanism. Guarded by
`[S82]` (`tests/check-cleanup-default.sh`: `nextflow config` resolves
`cleanup = false` with no profile and under `-profile test`). The
original analysis follows, for the record.

`nextflow.config` sets `cleanup = true` at top level (the `test` profile
overrides it to `false` so nf-test can read work files). `cleanup`
deletes the per-task `work/` directories on **successful** completion.
For the target use case — multi-hour/day runs on very large datasets —
this has two costs the README already warns about but the default works
against:

- `-resume` across separate invocations needs `work/` to survive; with
  `cleanup = true` a successful run leaves nothing to resume from, so a
  follow-up run (e.g. adding Part C, or re-running after a downstream
  tweak) recomputes everything.
- a run that "succeeded but produced something wrong" has no work
  directories left for post-mortem inspection.

A footgun-by-default sits awkwardly against the project's "hard to
misuse" goal, and nf-core ships `cleanup` **off** by default for exactly
these reasons.

**Question:** keep `cleanup = true`, flip it off, or make it a knob?

1. **Status quo (`cleanup = true`).** Tidy by default; defeats `-resume`
   and forensics, as above. Documented but easy to get bitten by.
2. **Default `cleanup = false` (recommended).** Matches nf-core; serves
   the large-dataset / `-resume` workflow out of the box; post-mortem
   debugging stays possible. Cost: `work/` accumulates on scratch and
   must be cleaned by hand (the README already documents `rm -rf work/`),
   so it trades silent disk reclamation for resume-correctness.
3. **Expose `--cleanup` (param, default off).** Same default as option 2
   but lets a user opt back into auto-clean for throwaway runs. Cheap;
   composes with option 2.

**Proposed resolution:** option 2, optionally with the option-3 knob.
This changes a default, so it warrants a CHANGELOG entry and a one-line
README note ("`work/` is no longer auto-deleted; clean it by hand, or
set `cleanup = true` in a `-c` override / `--cleanup` for throwaway
runs"). The `test` profile's explicit `cleanup = false` is unaffected.
New `[Sxx]`: "`cleanup` defaults to `false`; a successful run retains
`work/` so `-resume` works across invocations." Testable via
`nextflow config` resolving `cleanup = false` with no profile (mirrors
the `tests/check-*.sh` config-resolution style).

Until D13 lands, `cleanup = true` stands.


## D14 — Offline / air-gapped container path

**Revises:** `[S08]` (extends it with `[S83]`); extends D10
**Status:** `resolved` — options 1 + 3, implemented 2026-06-21 (`[S83]`)

**Resolution (2026-06-21):** options 1 + 3, both realised **without new
code in the engine profiles**. (1) Documented the build-once-online,
run-offline-from-cache recipe in the README — it uses the existing [S08]
Wave profiles unchanged. (3) The "site-supplied image" override is
realised through the existing `-c site.config` mechanism ([S75]) rather
than a dedicated `--container` flag: a site enables the engine and sets
`process.container` in its `-c` and composes it with `-profile slurm`
(the executor profile), so Wave is never turned on. A dedicated
`--container` param was **rejected**: a lone flag could not also disable
Wave (the engine profiles hard-enable it), and a param read inside a
profile block evaluates at config-parse time so it would not reliably
reflect a CLI/`-c` override anyway (the same eval-order limitation that
sank D13's `--cleanup` knob). The `-c` path sidesteps both issues — it
sets the engine, the image, and (by not using an engine profile) leaves
Wave off, all in one place. Covered by `[S83]`
(`tests/check-offline-container.sh`: a site `-c` resolves
`process.container` + the engine with no `wave.enabled`; a plain
`-profile slurm` resolves no container). `conf/site.config.example` and
the README "Air-gapped clusters" section document both paths. Option 2
(freeze to a project-owned registry) remains the documented upgrade path
from D10 if a site cannot build its own image. The original analysis
follows, for the record.

D10 chose Wave-from-`environment.yml` with no freeze and no registry,
explicitly on the assumption that "the two adopting labs run on slurm
with **outbound network on the compute nodes**." For the broader
"portability to other HPC setups" goal, that assumption does not hold
everywhere: many clusters air-gap their compute nodes (and sometimes the
login node), so Wave — which resolves and builds the image at task start
— cannot reach the network from where tasks run. Such a site is left
with `-profile modules` as the only option, since `conda` may also be
blocked. There is currently no way to point the workflow at a
pre-built image.

**Question:** what offline story should the container profiles support?

1. **Document a "build once online, run offline" recipe (recommended,
   cheap).** On a connected node, run the pipeline (or `nextflow inspect`
   / a Wave pull) once so the Wave-built image lands in
   `singularity.cacheDir` / `apptainer.cacheDir` on shared scratch; then
   compute nodes reuse the cached `.sif` with no further network. Pure
   documentation + a note in `conf/site.config.example`.
2. **Freeze to a project-owned registry (D10 option 2).** Flip
   `wave.freeze` + `wave.build.repository` so Wave pushes a permanent
   image to ghcr.io. Self-owned artefacts, but needs a registry + token
   and a release step — heavier; deferrable.
3. **Add an optional `--container <uri>` override (recommended,
   cheap).** Let a site set `process.container` to a pre-pulled `.sif`
   path or a registry URI, bypassing Wave entirely. One param + a
   `withName`-free `process.container = params.container ?: null`,
   scoped so it is null (Wave path unchanged) unless set.

**Proposed resolution:** options 1 + 3 together — they are both low-cost
and cover the realistic air-gapped cases without committing to running a
registry (option 2 stays the documented upgrade path from D10). New
`[Sxx]`: "an air-gapped site can run from a pre-built image via the
container cache (build-once recipe) or a `--container` override, with no
network from compute nodes." Wiring is testable by `nextflow config`
(the override resolves to `process.container`); execution stays a manual
cluster smoke test, as for the rest of `[S08]`.

Sub-questions for the human:
- Is `--container` a single image for all tools (matches the single
  `environment.yml` image Wave builds today), or per-process? Single is
  almost certainly right here.
- Should the project commit to publishing a frozen image per release
  (option 2) for sites that cannot build their own, or leave that to the
  site? This is the only part that adds ongoing maintenance.

Until D14 lands, the Wave-online assumption from D10 stands and
air-gapped sites use `-profile modules`.


## D15 — Group step logs under a dedicated `<outdir>/logs/` tree

**Blocks:** no new `[Sxx]` (refines the existing layout authority `[S71]`)
**Revises on resolution:** `[S19]`, `[S04]`, `[S45]`, `[S59]`, `[S71]`, and
Part C's `_taxonomy.log` publish location (`[S61]`/`[S50]` sintax path)
**Status:** `superseded by [D16]` — the dedicated `logs/` tree stands, but
its internal `per_sample` / `occurrence_table` (data-mirroring) sub-layout
was replaced by the stage-based `part_a` / `part_b` / `part_c` layout
(2026-06-22, before any release). Originally: `resolved` — option 2
selected (2026-06-22).

`[S71]` routes every published artefact under `<outdir>`, but the
per-step **log** files sit interleaved with the data files they
describe: Part A logs land in `<outdir>/per_sample/` next to
`<sample>.{fas,qual,stats}` (`[S19]`), and the Part B six step logs
(`[S45]`) plus Part C's `_taxonomy.log` land in
`<outdir>/occurrence_table/` next to the tables and FASTA. A user who
wants only the results, or only the logs, has to filter by extension
across two directories.

**Question:** keep logs interleaved with data, or collect them under a
single dedicated directory?

1. **Status quo.** Logs stay beside the data they describe; no change.
2. **Dedicated `<outdir>/logs/` tree mirroring the data subdirs
   (recommended).** Step logs move to `<outdir>/logs/per_sample/`
   (Part A) and `<outdir>/logs/occurrence_table/` (Part B/C). Data
   files stay where they are. A second helper `log_dir('<subdir>')`
   (sibling of `publish_dir`) keeps routing in one place. Nextflow's
   `pipeline_info/` reports are unaffected — they are not step logs.
3. **Flat `<outdir>/logs/`.** All logs in one directory, no
   sub-structure. Simpler tree but mixes per-sample and global logs;
   filenames already disambiguate them.

**Resolution:** option 2 (2026-06-22, confirmed by the developer). The
existing `per_sample` / `occurrence_table` subdir names are reused under
`logs/` to minimise spec and test churn. This is a breaking change to
output locations for anyone scripting against the old log paths; it is
grouped with the other `[S71]` breaking changes in the migration note.

Consequences for the spec:
- `[S71]` gains the `<outdir>/logs/{per_sample,occurrence_table}/`
  layer and the `log_dir('<subdir>')` helper.
- `[S19]` / `[S04]`: Part A (and shadow) step logs move to
  `<outdir>/logs/per_sample/`; the data files (`.fas` / `.qual` /
  `.stats`) stay in `<outdir>/per_sample/`.
- `[S45]`: the six Part B step logs move to
  `<outdir>/logs/occurrence_table/`.
- `[S59]`: the closed `<outdir>/occurrence_table/` whitelist drops the
  six step logs (now data-only: `_table.tsv`, `_table.fas`, and the
  Part C assigned tables); a complementary guarantee is that
  `<outdir>/occurrence_table/` contains **no** `*.log`.
- Part C's `_taxonomy.log` moves to `<outdir>/logs/occurrence_table/`.


## D16 — Organise `<outdir>/logs/` by pipeline stage, not by data directory

**Blocks:** no new `[Sxx]` (refines `[S71]`; supersedes [D15]'s sub-layout)
**Revises on resolution:** `[S19]`, `[S04]`, `[S45]`, `[S59]`, `[S71]`,
`[S86]`, and Part C's `_taxonomy.log` location
**Status:** `resolved` — option 2 selected (2026-06-22); spec + code in
this commit series

[D15] created the dedicated `<outdir>/logs/` tree but named its
sub-directories after the **data** directories they mirror —
`logs/per_sample/` (Part A) and `logs/occurrence_table/` (Part B + Part
C). That mirroring has three frictions: `logs/occurrence_table/` holds
no occurrence table and commingles Part B with Part C; the `[S86]`
read-count summary — a Part A artefact — was published into
`logs/occurrence_table/`, so a Part A-only run created an
`occurrence_table` log directory with no occurrence table; and there is
no single place to see "everything one stage did".

**Question:** keep the data-mirroring sub-layout ([D15]) or regroup the
logs by the pipeline stage that produced them?

1. **Status quo ([D15]).** `logs/per_sample/` + `logs/occurrence_table/`.
2. **Stage-based layout (recommended).** Regroup as:
   - `logs/part_a/per_sample/` — Part A per-sample step logs ([S19],
     including `_notmerged` shadow siblings [S04]);
   - `logs/part_a/<basename>_read_counts.tsv` — the Part A-wide
     read-count summary ([S86]), a sibling of `per_sample/` (it is
     project-wide, not per-sample);
   - `logs/part_b/` — Part B's six step logs ([S45]);
   - `logs/part_c/` — Part C's `<basename>_taxonomy.log` (the sintax
     path's vsearch `--log`, or the stampa path's per-chunk
     `vsearch.log` slices gathered into one file).
   Only stages that actually run produce a directory, so the tree is
   truthful about what executed.
3. **Flat `logs/`** with stage-prefixed filenames. Rejected — loses the
   per-sample grouping and reads worse than sub-directories.

**Resolution:** option 2, snake_case (`part_a` / `part_b` / `part_c`) to
match the repository's `modules/local/part_*` filesystem convention.

Trade-off accepted: the logs tree no longer mirrors the **data** tree
(`per_sample/`, `occurrence_table/` are unchanged — renaming those is
`[S71]`'s much larger blast radius and out of scope), so finding a
result's build log now goes by stage (table → Part B → `logs/part_b/`)
rather than by matching leaf name. The data-dir/log-dir asymmetry is the
price of the more honest, stage-centric grouping. Done now, before any
release, so no published output path is broken twice.


## D17 — Which hard-coded threshold values stay fixed (not exposed)

**Blocks:** no new `[Sxx]` (records the boundary of the parameter
surface added for `[S17]`, `[S35]`, `[S42]`, `[S88]`–`[S90]`)
**Status:** `resolved` — fixed-constant list confirmed (2026-06-24)

Exposing the tunable thresholds (`[S88]`–`[S90]`, plus `[S17]`/`[S35]`/
`[S42]` promotions) raised the question of where to stop. The following
values were reviewed and deliberately **left hard-coded** because they
are structural to the method rather than knobs a user would tune:

1. **swarm `--differences 1`** (`global_clustering`,
   `list_local_clusters`). `d = 1` *is* the ASV definition this pipeline
   implements; changing it changes the method, not a threshold.
2. **vsearch `--id 1.0`** in `search_for_terminal_gaps`. The step detects
   sequences identical modulo terminal gaps — only `1.0` is meaningful.
3. **vsearch `--fastq_maxns 0`** in `filter_and_convert_to_fasta`. "Any
   ambiguous base drops the read" is intrinsic to exact ASVs (`[S65]`).
4. **vsearch sentinels** `--maxaccepts 0`, `--fasta_width 0`,
   `--rowlen 0`, and `--differences`/`--iddef` defaults that are not
   abundance/identity thresholds. These select "no limit" / "unwrapped"
   behaviour, not a tunable cutoff.

**Resolution:** keep 1–4 fixed. If a future use-case needs one of them,
re-open this decision and promote it the same way (`params` default +
schema range + `[Sxx]` + test) rather than editing the module in place.
