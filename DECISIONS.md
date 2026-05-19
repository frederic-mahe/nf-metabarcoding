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
**Status:** `resolved`

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

**Blocks:** `[S48]`, `[S51]`
**Status:** `open`

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
