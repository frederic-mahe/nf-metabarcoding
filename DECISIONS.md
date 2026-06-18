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
**Status:** `partial` — sub-question 2 resolved, sub-question 1 still open

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

Sub-question 1 (the CLI flag toggling between table-input and
fasta-input modes) is still open and continues to block `[S48]`.

Two sub-questions still need a human answer before Part C can move
beyond skeleton:

1. **Which CLI flag toggles between table-input and fasta-input
   modes?** Two candidates:
    - `--occurrence_table /path/to/table.tsv` to consume Part B's
      `_table.tsv` and extract a fasta on the fly (via
      `extract_fasta_sequences_from_occurrence_table`); else
      `--fasta_input /path/to/representatives.fas` for a
      stand-alone fasta input. The two flags are mutually
      exclusive.
    - reuse the existing `--fasta_folder` for the fasta-input
      case; introduce a single `--occurrence_table` flag for the
      table-input case.

2. **How should Part C publish its result when the input is a
   fasta file?** With no occurrence table to splice back onto, the
   options are:
    - emit a stand-alone TSV with just the taxonomy columns
      (`amplicon\tabundance\tidentity\ttaxonomy\treferences`) —
      same shape as the legacy `*.results` file;
    - synthesise a minimal occurrence table from the fasta
      headers (no per-sample columns) so the output shape is
      identical to the table-input case;
    - fail at startup with a message asking the user to provide
      an occurrence table.

   A separate but related question: when Part B and Part C run
   end-to-end, should Part C **overwrite** Part B's `_table.tsv`
   in place, or publish a sibling file (e.g.
   `<basename>_taxonomy.tsv`) so the unannotated Part B output is
   preserved alongside the annotated one?

Until D04 lands, the Part C `[Sxx]` tests stay tagged `pending` and
the workflow stub exposes no user-visible CLI surface for the
ambiguous parts.


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
**Revises on resolution:** `[S10]`, `[S11]`, `[S12]`, `[S27]` and the
discovery processes `discover_inputs` / `discover_part_b_fasta`
**Status:** `proposed` — **option 3 selected (2026-06-18); spec + code pending**

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
