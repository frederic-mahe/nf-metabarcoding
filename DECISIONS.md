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

**Blocks:** `[S04]`, `[S05]`
**Status:** `open`

When paired-end reads fail to merge, [`SPECIFICATIONS.md`](SPECIFICATIONS.md)
calls for joining R1 and R2 with `N`s, then converting `N→A` before
feeding the sequence to swarm (which does not handle `N`), then back
to `N` for the occurrence table. Open sub-questions:

- **Which positions get rewritten?** Every `N`, or only the
  artificial joiner `N`s? Real reads may legitimately contain `N`s.
- **How is the round-trip tracked?** Per-sequence sidecar (cluster ID
  → list of N positions), or do we forbid downstream `A`s at known
  join positions?
- **Failure mode** if a sequence acquires unexpected `N`s mid-pipeline.


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

**Blocks:** `[S14]` (the warning bullet `[S13]` is testable today)
**Status:** `open`

Three plausible answers:

1. **Merge** — concatenate reads from same-named inputs (fast,
   destructive)
2. **Refuse** — abort with a clear error (safest, but breaks
   discovery if two folders happen to overlap)
3. **Suffix-disambiguate** — append a hash or folder name to one of
   the collisions (verbose, never silently merges)

Sub-decisions follow from the choice (e.g. if (3), what suffix
algorithm; if (1), how is provenance recorded).
