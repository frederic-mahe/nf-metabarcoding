# Migration plan — `nextflow_schema.json` + `nf-schema` validation

Status: **proposal, not yet approved**. This document is a plan only;
no validation behaviour has changed. It touches specified behaviour
(`[S57]`, `[S72]`, `[S74]`, `[S58]`, `[S61]`, `[S65]`, `[S63]`) and the
test strategy for several `[Sxx]`, so per `CLAUDE.md` it needs a human
review and a `DECISIONS.md` entry (proposed `D11`, see Phase 0) before
any code lands.


## 1. Goal

Replace the hand-rolled Groovy parameter-validation layer with the
standard nf-core `nf-schema` plugin driven by a `nextflow_schema.json`
file, **for the parts that fit a JSON-schema model** — types, enums,
numeric ranges, string patterns, and `--help` generation. Keep in
Groovy the validation that genuinely cannot be expressed as a static
per-parameter schema (file-content sniffing, tilde expansion,
run-mode-dependent requirements).

Why: a `nextflow_schema.json` is the de-facto nf-core interface
contract — it gives a single declarative source for types/defaults,
auto-generated `--help`, IDE/Seqera Platform/`nf-core` tooling support,
and removes ~200 lines of bespoke validation Groovy
(`numeric_param_spec` / `check_numeric_param` / `validate_numeric_params`
+ several enum asserts + the `usage()` heredoc).

Non-goal: this is **not** a rewrite of the run-mode dispatch in
`main.nf`, nor of `bin/parse_samplesheet.py`. Those stay (see §4, §6).


## 2. Current state (what exists today)

All validation is pure Groovy, no plugins:

- `nextflow.config` `params {}` block — central defaults, heavily
  commented (protected by `CLAUDE.md`); also profile-injected params
  inside `slurm` / `modules` / `demo` / `local` profiles.
- `modules/local/functions.nf`:
  - `usage()` — hand-tuned grouped `--help` heredoc (`[S57]`).
  - `numeric_param_spec()` / `check_numeric_param()` /
    `validate_numeric_params()` — numeric ranges (`[S72]`).
  - `check_primer_format()` — IUPAC regex (`[S74]`).
  - `check_reference_format()` / `reference_first_header()` — FASTA
    header content sniff, gzip-aware (`[S73]`).
  - `normalize_path()` / `lookup_user_home()` — `~` expansion (`[S60]`).
  - `validate_params()` — enum asserts (`publish_mode` `[S58]`,
    `taxonomy_method` `[S61]`, `hash_function` `[S65]`),
    `join_padding_length` positivity (`[S63]`), `majority_assignment`
    × `sintax` incompatibility (`[S66]`), mode-selector mutual
    exclusivity (`[S02]`/`[S70]`), `--results_folder` deprecation warning.
  - `samplesheet_profile()` — header sniff to route fastq/fasta (`[S70]`).
  - `resource_size_warnings()` — slurm size-hint warnings (`[S79]`).
- `main.nf` workflow body — `--help` short-circuit, calls
  `validate_params()`, then **mode-specific** inline asserts
  (`fastq_folder` required, `forward_primer`/`reverse_primer` required
  unless `--no_trimming`, `reference_dataset` per `taxonomy_method`,
  `project_name` when Part B, `occurrence_table` for Part C).
- `modules/local/validate_samplesheet.nf` → `bin/parse_samplesheet.py`
  — structural + semantic samplesheet checks (`[S70]`).

Tests:
- Unit tests of the Groovy helpers:
  `tests/functions/check_numeric_param.nf.test`,
  `check_primer_format.nf.test`, `check_reference_format.nf.test`,
  `normalize_path.nf.test`, `resource_size_warnings.nf.test`,
  `effective_outdir.nf.test`.
- `tests/main.nf.test` — ~31 cases assert on **exact error-message
  text** (`"must be one of"`, param names) and on `--help` markers
  (`"Usage:"`, specific flag strings).
- `tests/COVERAGE.md` maps every `[Sxx]` → test file; `coverage-gate.sh`
  enforces the link.


## 3. What `nf-schema` can replace (in scope)

`nf-schema` (nf-core plugin) provides `validateParameters()`,
`paramsHelp()`, `paramsSummaryLog()`, and `samplesheetToList()`. The
JSON-schema vocabulary covers:

| Current Groovy | `nextflow_schema.json` mechanism | Spec |
|---|---|---|
| `usage()` heredoc | `paramsHelp()` from per-property `description` / `help_text`, grouped by `$defs` | `[S57]` |
| numeric ranges | `type: integer/number` + `minimum`/`maximum`/`exclusiveMinimum` | `[S72]` |
| `fastq_encoding ∈ {33,64}` | `enum: [33, 64]` | `[S72]` |
| `publish_mode` enum | `enum: [copy, copyNoFollow, link, move, rellink, symlink]` | `[S58]` |
| `taxonomy_method` enum | `enum: [stampa, sintax]` | `[S61]` |
| `hash_function` enum | `enum: [sha1, md5]` | `[S65]` |
| `join_padding_length ≥ 1` | `type: integer, minimum: 1` | `[S63]` |
| primer IUPAC regex | `pattern: "^[ACGTURYSWKMBDHVNIacgturyswkmbdhvni]{3,}$"` (validity only; *requiredness* stays Groovy — it is mode-dependent) | `[S74]` |
| path params | `format: file-path` / `directory-path` (see §4 caveat about `exists` + tilde) | — |

Defaults: `nextflow_schema.json` can also *carry* defaults, but we will
**keep `nextflow.config`'s `params {}` block as the source of defaults**
(its comments are protected by `CLAUDE.md` and carry `[Sxx]` rationale).
The schema mirrors them; a new sync test (§5) guards drift, in the
spirit of the existing `manifest.version == CITATION.cff` test (`[S77]`).


## 4. What must stay in Groovy (out of scope for the schema)

These cannot be expressed as a static per-parameter JSON schema, or are
clearer/safer left as-is:

1. **`check_reference_format()` (`[S73]`)** — reads and decompresses the
   reference FASTA to inspect its first header (stampa vs sintax). Pure
   content inspection; no schema equivalent. **Keep.**
2. **`normalize_path()` / `~` expansion (`[S60]`)** — schema validation
   does not tilde-expand. See the ordering caveat below. **Keep.**
3. **Mode-specific requiredness** — `forward_primer`/`reverse_primer`
   required *unless* `--no_trimming`; `reference_dataset` vs
   `reference_dataset_sintax` chosen by `taxonomy_method`;
   `project_name` required when Part B runs; `occurrence_table` for
   Part C standalone; `fastq_folder` required only without `--input`.
   These depend on the **resolved run mode**, which is computed in the
   `main.nf` dispatch. JSON-schema `if/then/else` / `dependentRequired`
   can express *some* of this but not the dynamic dispatch, and it
   would split one decision across two files. **Keep in Groovy.**
4. **Mode-selector mutual exclusivity (`[S02]`/`[S70]`)** — expressible
   as schema `oneOf`, but the current assert produces a precise message
   naming exactly which selectors were set. **Keep in Groovy** (low
   value to move, higher risk of a worse message).
5. **`majority_assignment` × `sintax` incompatibility (`[S66]`)** —
   a cross-param rule; could be schema `allOf/if/then` but is clearer
   inline. **Keep in Groovy.**
6. **`resource_size_warnings()` (`[S79]`)** — runtime warning keyed on
   the active `workflow.profile`, not a param constraint. **Keep.**
7. **`--results_folder` deprecation warning (`[S71]`)** — a warning, not
   a rejection. **Keep** (or model as schema `deprecated`, but the
   custom message is more useful).
8. **Samplesheet semantics (`[S70]`)** — `bin/parse_samplesheet.py`
   does single-end inference, reserved `_notmerged` suffix rejection,
   sibling `qual`/`stats` defaulting, and emits a normalized TSV.
   `nf-schema`'s `samplesheetToList` + a `schema_input.json` validates
   *structure* only and cannot do the inference/defaulting. **Defer**
   to a separate, optional later effort (§6 Phase 4); out of scope here.


### Critical caveats / gotchas

- **`validateParameters()` validates *every* param against the schema**
  and warns (or fails, depending on config) on params not declared in
  it. The schema must therefore also declare **all profile-injected
  params**: `slurm_queue`, `slurm_account`, `slurm_clusterOptions`,
  `slurm_queue_size`, `dataset_size_gb`, `reference_size_gb`,
  `max_cpus`, `max_memory`, `max_time`, `module_vsearch`/`_swarm`/
  `_cutadapt`/`_mumu`, `forward_primer`, `reverse_primer` — or they
  must be added to the plugin's ignore list. Missing one turns every
  run under that profile into a warning storm.
- **`max_memory` / `max_time`** are Nextflow `MemoryUnit`/`Duration`
  literals (`128.GB`, `240.h`). Model as `type: string` with a
  `pattern` (e.g. `^\d+(\.\d+)?\.?(K|M|G|T)?B$`), not a number.
- **Ordering vs tilde expansion (`[S60]`)**: do **not** mark
  tilde-accepting path params with `format: file-path-exists` /
  `exists: true`. `validateParameters()` runs at workflow start and
  would stat the *unexpanded* `~/...` path and wrongly fail. Use
  `format: file-path` (shape only) and let the existing
  `file(normalize_path(...), checkIfExists: true)` at the use site do
  existence. Keep `validateParameters()` ordering compatible with the
  `--help` short-circuit (call help first, validate after).
- **Error-message text changes.** `nf-schema` emits its own
  multi-line validation report; it will **not** match the current
  `"--<name> must be one of ..."` strings. Every `tests/main.nf.test`
  case that asserts on those strings (≈31 cases reference message text;
  `--help` markers too) must be updated. Those are human-reviewed tests
  → explicit authorization required before editing (`CLAUDE.md` rule 2).
- **Plugin fetch on air-gapped HPC.** `nf-schema` is downloaded to
  `$NXF_PLUGINS_DIR` on first use; offline clusters must pre-seed it.
  The project targets slurm/HPC, so document this (and pin the version).
- **Loss of fast unit tests.** `check_numeric_param` / `check_primer_format`
  are unit-tested directly today (`tests/functions/*.nf.test`).
  Schema validation only runs inside a workflow, so those checks move
  to slower `tests/main.nf.test` integration cases. Net: fewer, slower
  tests for the same `[Sxx]`. Call this out for the human reviewer —
  it may be a reason to keep `[S72]`/`[S74]` in Groovy after all.


## 5. Proposed end state — file changes

- **`nextflow.config`**
  - add `plugins { id 'nf-schema@<pinned>' }` (confirm latest 2.x
    compatible with `nextflowVersion >= 25.04.0`).
  - add a `validation { ... }` block: parametrise help, and set the
    ignore list / monochrome / lenient options as needed.
  - keep the `params {}` defaults block unchanged (source of defaults).
- **`nextflow_schema.json`** (new, repo root) — `$defs` groups mirroring
  `usage()` sections (Entry-point / Part A / Part B / Part C / Runtime /
  Slurm / Modules), every param typed with `description` + `help_text`,
  enums/ranges/patterns per §3, all profile params declared (§4 caveat).
- **`modules/local/functions.nf`**
  - remove `numeric_param_spec` / `check_numeric_param` /
    `validate_numeric_params` (`[S72]`), the `usage()` heredoc (`[S57]`),
    `check_primer_format` (`[S74]`), and the enum/positivity asserts in
    `validate_params` (`[S58]`/`[S61]`/`[S65]`/`[S63]`) **only after**
    their schema equivalents are proven green.
  - **keep** `check_reference_format`, `normalize_path`,
    `samplesheet_profile`, `resource_size_warnings`, the mode-selector
    and `majority_assignment` asserts, the deprecation warning.
- **`main.nf`**
  - replace the `--help` short-circuit body with `paramsHelp()`
    (`[S57]`); optionally add `paramsSummaryLog()` to the run banner.
  - call `validateParameters()` (after the help short-circuit) in place
    of the removed-from-`validate_params` checks; keep the
    Groovy-resident checks (§4) and all mode-specific inline asserts.
- **Tests**
  - new `tests/main.nf.test` cases (or updated existing ones, with
    authorization) asserting the schema-driven rejections + new `--help`.
  - new sync test (Python or nf-test): every key in
    `nextflow_schema.json` ↔ `params {}` defaults agree (drift guard,
    cf. `[S77]`).
  - retarget/remove `tests/functions/check_numeric_param.nf.test` and
    `check_primer_format.nf.test` (human-authorized).
- **`tests/COVERAGE.md`** — update the test pointers for `[S57]`,
  `[S72]`, `[S74]`, `[S58]`, `[S61]`, `[S65]`, `[S63]` in the same
  commit as the code (coverage-gate requirement).
- **Docs** — `README.md` / `SPECIFICATIONS.md` notes where the help
  text or message format is now schema-generated.


## 6. Phased execution (each phase ends green: `nf-test --tag ci`,
`coverage-gate.sh`, `shellcheck`/`flake8` on touched files)

- **Phase 0 — decide.** Open `DECISIONS.md` `D11`: "adopt nf-schema?"
  capturing the two judgement calls — (a) accept schema-generated
  `--help`/message text (rewrites human-reviewed test assertions);
  (b) accept moving `[S72]`/`[S74]` from fast unit tests to integration
  tests. **Stop for human review** (TDD rule: no tests against
  undefined behaviour; `CLAUDE.md` golden rule).
- **Phase 1 — additive, no removals.** Add the plugin +
  `nextflow_schema.json` mirroring today's params (all profile params
  declared). Wire `validateParameters()` and `paramsHelp()` *alongside*
  the existing Groovy (do not delete anything yet). Add the schema↔params
  sync test. Verify the full suite still passes unchanged. This proves
  the plugin loads, the schema covers every param, and offline plugin
  fetch is understood — with zero behaviour change.
- **Phase 2 — switch `--help` (`[S57]`).** Replace `usage()` with
  `paramsHelp()`. Update the `--help` assertions in `tests/main.nf.test`
  (authorized). Update COVERAGE. Remove `usage()`.
- **Phase 3 — retire schema-covered Groovy validators.** Remove the
  numeric/enum/pattern/positivity checks (`[S72]`/`[S74]`/`[S58]`/
  `[S61]`/`[S65]`/`[S63]`) now that the schema enforces them. Update the
  affected `tests/main.nf.test` message assertions and retarget the two
  `tests/functions/*` unit tests (authorized). Update COVERAGE.
- **Phase 4 — (optional, separate proposal) samplesheet schema.**
  Evaluate `samplesheetToList` + `schema_input.json` for the *structural*
  layer of `[S70]`, keeping `parse_samplesheet.py` for inference/
  defaulting. Likely deferred — listed for completeness, not part of
  this migration.

Each phase is a small, self-contained commit (or a few), on this
`tmp_*` branch, with `Co-Authored-By: Florian Filloux`. Human merges to
`dev`.


## 7. Risks summary

- Test churn on human-reviewed assertions (message text + `--help`) —
  needs authorization; mitigated by phasing (Phase 2/3 isolate it).
- Coarser test granularity for `[S72]`/`[S74]` (unit → integration).
- Plugin availability on air-gapped compute nodes — pin + pre-seed.
- Schema/params default drift — mitigated by the sync test.
- `validateParameters()` strictness vs profile-injected params and the
  tilde-expansion ordering — both addressed in §4 and validated in
  Phase 1 before any removal.

If the human reviewer judges the test-granularity loss (§4 caveat) or
the `--help` text churn too costly, a smaller alternative is: adopt the
schema for `--help` generation and Platform/tooling interop **only**
(Phases 1–2) and keep all value validation in Groovy (skip Phase 3).
