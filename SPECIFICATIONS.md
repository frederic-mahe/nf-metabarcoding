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
- `[S02]` each part can be run separately, or all at once
  - **Pass when:** `--only=A|B|C` (or equivalent flag) runs that part
    in isolation; default runs all three
- `[S03]` fastq files can be paired-end or single-end, compressed
  (`.gz`, `.bz2`) or not
  - **Pass when:** the same input data in any of these forms produces
    byte-identical Part A outputs
- `[S04]` when processing paired-end fastq files, reads that can be
  merged are processed normally; reads that cannot be merged are
  collected via `vsearch --fastq_mergepairs --fastqout_notmerged_fwd`
  and `--fastqout_notmerged_rev` (reads stay in sync) and routed
  through a parallel **shadow** Part A pipeline:
    1. join R1/R2 with `vsearch --fastq_join` (8-N padding — the
       tool default);
    2. trim primers (when `--no_trimming` is false) and convert to
       fasta with `vsearch --fastq_filter --fastq_maxns 8` so the
       join-padding `N`s survive the filter;
    3. dereplicate and extract ee values per the normal pipeline; the
       published `.fas` retains the `N`s;
    4. replace every `N` with `A` **inline** just before swarm (swarm
       rejects `N`s). The masked fasta is **not** published; only
       sequence lines are rewritten so SHA1 headers (computed in step
       2) are unchanged and the published `.stats` IDs match the
       `.fas` IDs.
  - shadow-pipeline sample IDs are `<sampleId>_notmerged`; published
    artefacts are `<sampleId>_notmerged.{fas,qual,stats}` plus the
    same per-step logs as the normal pipeline with the `_notmerged`
    prefix (`<sampleId>_notmerged_{merging,trimming,dereplicating,clustering}.log`).
  - **Pass when:** running Part A on a paired-end fixture whose reads
    cannot overlap produces non-empty `<sampleId>_notmerged.{fas,qual,stats}`
    in `params.fastq_folder`, the published `.fas` contains the 8-`N`
    join padding, and every `_notmerged_<step>.log` from `[S19]` is
    present and non-empty.
- `[S05]` unmerged-pair clusters appear in the occurrence table with a
  per-sample marker (working name: `sampleID_partial`)
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) — final marker name


## Interfaces between parts

Each part is a directory contract so parts can be tested in
isolation.

- **Part A → Part B**
  - `<sampleId>.fas` — dereplicated FASTA, vsearch
    `--sizeout`-style abundance annotation, `--fasta_width 0`
  - `<sampleId>.qual` — TSV with columns `sha1\tee\tlength`, sorted
    by length then SHA1
  - `<sampleId>.stats` — swarm `--statistics-file` output, filtered
    to clusters with > 2 reads
- **Part B → Part C**
  - one occurrence table (see schema below)
- **Part C output**
  - the occurrence table with taxonomy column updated (no rows
    added or removed)


## Occurrence table schema

- **single-table mode** (default)
  - columns: `cluster_id`, `sequence`, `abundance_total`, one column
    per sample (`<sampleId>` — integer read count), `taxonomy` (Part
    C, empty until then)
  - rows: one per cluster; empty samples (`[S08]`) contribute a
    zero-filled column
- **two-table mode** (`--split-occurrence-table`)
  - `occurrences.tsv` (long format): `clusterId`, `sampleId`,
    `abundance`
  - `clusters.tsv`: `clusterId`, `sequence`, `abundance_total`,
    `taxonomy`


## Workflow requirements

- `[S06]` configurable via a config file or command-line parameters
  - **Pass when:** the same run can be reproduced from either a
    `nextflow.config`/`params.yml` or from `--key value` flags
- `[S07]` runs locally or on HPC with slurm
  - Not part of automated CI; covered by manual smoke tests on the
    cluster. Profiles live in `nextflow.config`.
- `[S08]` runs local or containerized applications (vsearch, swarm,
  cutadapt)
  - Not part of automated CI; covered by a `-profile docker` /
    `-profile singularity` manual smoke test.
- `[S09]` empty input samples must travel through and appear in the
  occurrence table (but not in the two-table mode long-format)
  - **Pass when:** an empty fastq pair produces a row in the
    occurrence table with `abundance = 0` for that sample (and no
    error)
- `[S10]` accepts a directory, or a list of directories (absolute or
  relative paths). `--fastq_folder` accepts a single path or a
  comma-separated list (`--fastq_folder a,b,c`); a `nextflow.config`
  may instead supply a Groovy list (`fastq_folder = ['a', 'b']`).
  - **Pass when:** every fastq file in every listed folder is
    discovered; identical inputs spread across two folders produce
    identical artefacts to the single-folder run.
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
- `[S15]` users can choose between a single occurrence table or a
  two-part table (long-format + per-cluster metadata)
  - **Pass when:** `--split-occurrence-table` toggles between the two
    schemas defined above
- `[S16]` expects demultiplexed fastq files
  - **Pass when:** documented in README; demultiplexing is out of
    scope (could be added later as a subworkflow)
- `[S17]` per-cluster minimum-read threshold defaults to > 2 reads
  - **Pass when:** `list_local_clusters` emits no row with
    `reads <= 2` on a fixture that contains a singleton cluster
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
- `[S19]` every Part A per-sample artefact is published into
  `params.fastq_folder` (the same folder(s) the input fastq files
  came from):
    - data files: `<sampleId>.fas`, `<sampleId>.qual`,
      `<sampleId>.stats`
    - per-step log files:
        - merging       → `<sampleId>_merging.log`
        - trimming      → `<sampleId>_trimming.log` (only when the
          trimming step runs — see `[S20]`)
        - dereplicating → `<sampleId>_dereplicating.log`
        - clustering    → `<sampleId>_clustering.log`
  - **Pass when:** running Part A on any sample produces all three
    data files and all four log files in `params.fastq_folder`,
    each non-empty (three logs when `--no_trimming` is set: no
    `_trimming.log`).
- `[S20]` `--no_trimming` toggle (default: `false`) skips the primer
  trimming step. The toggle and the primer parameters are mutually
  exclusive:
  - when `--no_trimming` is `true`, `forward_primer` and
    `reverse_primer` must be empty (not set on the CLI or in the
    config). The minimum-length and max-N filters that the trimming
    step used to enforce are taken over by
    `filter_and_convert_to_fasta` (vsearch `--fastq_minlen 32`,
    `--fastq_maxns <caller-supplied>`).
  - when `--no_trimming` is `false` (the default),
    `forward_primer` and `reverse_primer` are required (see
    `[S18]`).
  - **Pass when:** (a) running with `--no_trimming true` and empty
    primers succeeds and does **not** publish a `_trimming.log`;
    (b) running with `--no_trimming true` together with a non-empty
    `forward_primer` or `reverse_primer` exits non-zero and the
    error names the conflicting parameter.
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
  emits an augmented `(stats, swarms, representatives)` triple. The
  fastidious flag used by the upstream swarm run is propagated to
  this step (default: `--fastidious`) so the representatives output
  filename matches the swarm parameters.
  - **Pass when:** golden-file characterization tests for
    `bin/cluster_cleaver.py` reproduce the byte-exact output of the
    legacy `tmp/OTU_cleaver.py` on a fixture covering: (a) a cluster
    that splits on a sub-seed, (b) a cluster that does not split,
    (c) a candidate sub-seed below the per-sample-presence threshold,
    (d) a sort tie within a sub-cluster.
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
- `[S26]` Part B requires `--results_folder` (no default). The
  value is an absolute or relative path to the folder where Part B
  publishes every artefact. The workflow creates the folder (and
  any missing parent directories) at startup if it does not exist;
  an existing folder is reused as-is. Like `[S25]`, this parameter
  is **only** required for Part B.
  - **Pass when:** running Part B with `--results_folder` set to a
    non-existent absolute or relative path succeeds, creates the
    folder, and every Part B artefact lands inside it.
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
  `vsearch --derep_fulllength` with `--sizein --sizeout --fasta_width
  0`. The output is `<project_name>_<N>_samples.fas`; the vsearch
  log is published as `<project_name>_<N>_samples.log`.
  - **Pass when:** running on two `.fas` fixtures that share a
    sequence yields one record per distinct sequence; the
    `;size=N` value is the sum of the per-sample sizes.
- `[S32]` Part B's `global_clustering` runs swarm on the
  globally-dereplicated fasta from `[S31]` with `--differences 1
  --fastidious --usearch-abundance` and the four output flags
  `--internal-structure`, `--output-file`, `--statistics-file`,
  `--seeds`. Output filenames follow the
  `<project_name>_<N>_samples_1f.{swarms,stats,struct}` and
  `<project_name>_<N>_samples_1f_representatives.fas` scheme; the
  run log is `<project_name>_<N>_samples_1f.log`.
  - **Pass when:** the four output files plus the log exist and are
    non-empty for the documented fixture.
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
- `[S37]` Part B's `chimera_detection2` re-runs uchime_denovo on
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
  kept iff `chimera == "N" && ee/length <= 0.0002 && (abundance >= 3
  || spread >= 2)`. Sample columns appear in sorted order; empty
  samples (`[S09]`) contribute a zero-filled column.
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
  legacy lulu-recommended parameters (`--id 0.84 --iddef 1
  --maxaccepts 0 --query_cov 0.9 --maxhits 10
  --userfields query+target+id`). The output is the userout stream
  with the `;size=N;` annotation stripped from every column.
  - **Pass when:** the output is a 3-column TSV
    (`query\ttarget\tid`); no `;size=` annotations remain.
- `[S43]` Part B's `run_mumu` invokes the `mumu` binary
  (`>= 1.1.1`) with `--otu_table`, `--match_list`,
  `--new_otu_table`, and `--log`. The cleaned-up intermediate
  inputs (`_reduced.table`, `.match_list`) are not kept.
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
- `[S45]` Part B publishes one step-level log file per major step,
  named with the project-wide basename `<project>_<N>_samples`
  (the same construct as the final occurrence table, see `[S46]`)
  and a step suffix. The six step logs are:
    - `<basename>_dereplication.log` — `global_dereplication`
      vsearch log (`[S31]`)
    - `<basename>_clustering.log` — `global_clustering` swarm log
      (`[S32]`)
    - `<basename>_chimera_detection.log` — the canonical chimera
      run's stderr (`chimera_detection2`'s output, `[S37]`);
      `chimera_detection`'s pre-cleave stderr is kept internal
    - `<basename>_cleaving.log` — `cleaving`'s stderr from
      `bin/cluster_cleaver.py` (`[S22]`)
    - `<basename>_superstring_clustering.log` — combined stderr of
      `search_for_terminal_gaps` + `merge_substring_otus`
      (`[S38]`, `[S39]`)
    - `<basename>_post_clustering_curation.log` — `run_mumu`'s
      `--log` output (`[S43]`)
  - **Pass when:** running Part B on the documented fixture
    publishes the six step logs above into `params.results_folder`;
    each file is non-empty.
- `[S46]` Part B publishes its final occurrence table as
  `<project>_<N>_samples_table.tsv` (after `rebuild_post_mumu_table`
  and the size=0→1 awk hotfix). Intermediate tables produced along
  the way (`*.OTU.filtered.cleaved.table`, `*.nosubstringOTUs.table`,
  `*_raw_mumu.table`, `*.mumu.table`) are no longer published —
  only the final `_table.tsv` lands in `params.results_folder`.
  - **Pass when:** running Part B on the documented fixture
    publishes `<project>_<N>_samples_table.tsv` to
    `params.results_folder`; the file is non-empty and starts with
    the OTU-table header row; none of the intermediate
    `.OTU.filtered.cleaved.table`, `.nosubstringOTUs.table`,
    `_raw_mumu.table`, or `.mumu.table` files are present in the
    results folder.


## Part C — taxonomic assignment

Part C re-implements the legacy
[`stampa`](./tmp/stampa/) method on top of nextflow. Inputs: either
the Part B occurrence table (`[S46]`) or a standalone fasta file
of representative sequences; plus a reference dataset (fasta,
optionally compressed). Output: an occurrence table whose
`taxonomy` / `identity` / `references` columns have been updated
from placeholder values to real taxonomic assignments.

- `[S47]` Part C requires `--reference_dataset` (no default). The
  value is an absolute or relative path to a fasta file (plain
  `.fasta`/`.fas`, gzip-compressed `.gz`, or bzip2-compressed
  `.bz2`). The workflow aborts at startup if Part C is asked to
  run without `--reference_dataset`. The parameter is **only**
  required for Part C; Parts A and B do not consume it.
  - **Pass when:** running Part C without `--reference_dataset`
    exits non-zero and stderr names the missing parameter;
    supplying it via `-params-file` or `--key value` lets the
    workflow proceed.
- `[S48]` Part C accepts either an occurrence table (the
  `[S46]` `_table.tsv`) or a fasta file as its primary input.
  When given an occurrence table, the process
  `extract_fasta_sequences_from_occurrence_table` reads column 4
  (amplicon ID), column 2 (abundance), and column 10 (sequence)
  and emits a fasta of representatives with header
  `<amplicon>;size=<abundance>;`. When given a fasta, the
  extraction step is skipped.
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) D04 — which CLI
    flag toggles between the two modes (`--occurrence_table` vs
    `--fasta_input`), and how the updated table is reconstructed
    when the input is a fasta (no occurrence table to splice back
    onto).
- `[S49]` Part C's primary taxonomic-assignment path is a port of
  the legacy `stampa.sh` / `stampa_merge.py` pipeline:
    1. one `vsearch --usearch_global` self-search against the
       reference dataset with `--top_hits_only --output_no_hits
       --maxaccepts 0 --maxrejects 0 --notrunclabels --rowlen 0`
       and `--userfields query+id<iddef>+target`. The whole
       representatives fasta goes through a single vsearch
       invocation (no slurm array split — nextflow handles
       parallelism by process).
    2. `bin/stampa_merge.py` parses the userout, computes the
       last-common-ancestor taxonomy across top hits per
       amplicon, and emits a TSV with columns
       `amplicon\tabundance\tidentity\ttaxonomy\treferences`
       (same shape as `fake_taxonomic_assignment`'s `[S33]`
       output).
  - **Pass when:** **(skeleton phase)** the process and helper
    exist, their CLI parses without error, and a smoke test
    against a tiny reference fixture produces a TSV with the
    documented shape. Real-world LCA correctness is pinned by
    porting the legacy `stampa_merge.py` unit tests in a follow-up.
- `[S50]` Part C's shadow path uses `vsearch --sintax` against
  the **same** reference dataset and emits an alternative
  taxonomy TSV with the same column shape as the primary path
  (`[S49]`). Sintax-specific columns (bootstrap confidences) are
  collapsed into the `references` field; users opt in via
  `--taxonomy_method sintax` (default `stampa`).
  - **Pass when:** **(skeleton phase)** the sintax process exists
    and runs against a tiny reference fixture without error; the
    taxonomy column it emits is documented.
- `[S51]` Part C's `update_occurrence_table` splices the
  taxonomic assignment back onto the `[S46]` occurrence table by
  amplicon ID, overwriting the `identity`, `taxonomy`, and
  `references` columns. Rows are neither added nor removed. The
  output filename is `<project>_<N>_samples_table.tsv` (same
  construct as `[S46]`); when Part B and Part C run end-to-end,
  Part C's output overwrites Part B's `_table.tsv` in the same
  results folder.
  - **Blocked by:** [`DECISIONS.md`](DECISIONS.md) D04 — whether
    Part C overwrites in place or publishes a sibling file
    (e.g. `<basename>_taxonomy.tsv`).


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
