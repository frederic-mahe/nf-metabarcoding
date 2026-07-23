# Spec coverage

Each row maps an `[Sxx]` ID from
[`../SPECIFICATIONS.md`](../SPECIFICATIONS.md) to the test(s) that
cover it. **Status legend:**

- `done`    — test exists and passes
- `red`     — test exists and fails (the code is not yet there — this is the TDD red phase)
- `TODO`    — no test yet
- `n/a`     — meta or not testable in automated CI (e.g. slurm/HPC)
- `blocked` — spec needs clarification before a test can be written; see [`../DECISIONS.md`](../DECISIONS.md)

When a spec changes, update this table in the same commit. The
coverage gate (`bash tests/coverage-gate.sh`) checks that every
`[Sxx]` in SPECIFICATIONS appears here, and that every
`// COVERAGE: [Sxx]` in `tests/` references a known ID.


## Workflow structure

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S00]`| TDD discipline (meta)                                                      | `tests/coverage-gate.sh`                        | done    | —          |
| `[S01]`| three-part workflow (fastq→fasta, fasta→occurrence, taxonomic assignment)  | `tests/main.nf.test`, `tests/bin/reverse_complement.bats` | done    | —          |
| `[S02]`| each part can be run separately or all at once; the five input-mode selectors are mutually exclusive | `tests/main.nf.test`                            | done    | —          |
| `[S03]`| paired-end or single-end, compressed (gz/bz2) or uncompressed input        | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/main.nf.test` | done | — |
| `[S04]`| unmerged paired reads → shadow pipeline (A-padded join, no mask round-trip)| `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/join_notmerged.nf.test`, `tests/main.nf.test` | done   | D15, D16   |
| `[S05]`| unmerged clusters appear in occurrence table with per-sample marker        | —                                               | blocked | D01, D02   |


## Workflow requirements

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S06]`| read config file *or* command-line parameters                              | —                                               | n/a     | Nextflow's intrinsic `params` indirection; required-param branch covered by `[S18]` |
| `[S07]`| runs locally or on HPC (slurm)                                             | `tests/check-slurm-config.sh` (profile wiring: executor + resourceLimits + closures resolve, composes with conda, job arrays off by default); `tests/check-local-params.sh` (slurm-only account/size params stay defined on a no-profile local run); execution = manual cluster smoke test | done    | per-process resource tiers, dataset-scaled memory, retry escalation, `submitRateLimit` 50/1min, opt-in `--slurm_array_size` |
| `[S08]`| runs local *or* containerized applications (docker/podman/singularity/apptainer via Wave) | `tests/check-container-profiles.sh` (profile wiring); execution = manual cluster smoke test | done    | —          |
| `[S09]`| empty input samples travel through and appear in the occurrence table      | `tests/main.nf.test`, `tests/python/test_build_filtered_contingency_table.py` | done    | —          |
| `[S10]`| accept a directory or a list of directories (absolute or relative)         | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S11]`| auto-discover fastq files using common name patterns                       | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S12]`| auto-deduce sample names from fastq file names                             | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S13]`| abort if two or more samples share a derived ID; list duplicate file paths | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py`, `tests/main.nf.test` | done | — |
| `[S14]`| collision policy for same-named samples: refuse (sample IDs must be unique) | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py` | done    | —          |
| `[S15]`| export single occurrence table *or* two-part (long + metadata) table       | `tests/main.nf.test`                            | TODO    | —          |
| `[S16]`| expect demultiplexed fastq files                                           | —                                               | n/a     | —          |
| `[S17]`| per-cluster minimum-read threshold `--min_cluster_size` (default > 2 reads) | `tests/processes/part_a/list_local_clusters.nf.test`, `tests/bin/filter_swarm_stats.bats`, `tests/main.nf.test`   | done    | —          |
| `[S18]`| required params (forward/reverse_primer, fastq_folder) must be supplied    | `tests/main.nf.test`                            | done    | —          |
| `[S19]`| Part A steps emit per-sample `<sampleId>_<step>.log` files (data → `per_sample/`, logs → `logs/part_a/per_sample/`) | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/trim_primers.nf.test`, `tests/processes/part_a/dereplicate_fasta.nf.test`, `tests/processes/part_a/list_local_clusters.nf.test`, `tests/main.nf.test` | done | D15, D16 |
| `[S20]`| `--no_trimming` toggle skips primer trimming; mutually exclusive w/ primers| `tests/main.nf.test`                            | done    | —          |
| `[S21]`| unpaired fastq files skip the merging step                                 | `tests/main.nf.test`                            | done    | —          |
| `[S22]`| Part B re-cleaves global swarm clusters using per-sample sub-seed presence | `tests/python/test_cluster_cleaver.py`          | done    | —          |
| `[S23]`| `notmerged` reserved suffix — sample IDs ending in `notmerged` are rejected | `tests/python/test_reserved_keyword.py`, `tests/main.nf.test` | done | —          |
| `[S93]`| sample IDs restricted to safe charset `[A-Za-z0-9._-]` (shared validator) | `tests/python/test_sample_id.py` | done | — |
| `[S24]`| shadow pipeline 3'-strip via `vsearch --fastq_stripright` (default 30)     | `tests/processes/part_a/strip_reads.nf.test`    | done    | —          |
| `[S25]`| Part B requires `--project_name` (no default)                              | `tests/main.nf.test`                            | done    | —          |
| `[S26]`| Part B requires `--results_folder` (no default); auto-created if missing   | `tests/main.nf.test`                            | done    | —          |
| `[S27]`| Part B fasta channel excludes `_notmerged.fas`; sample IDs must be unique  | `tests/python/test_discover_fasta.py`, `tests/bin/discover_fasta.bats`, `tests/main.nf.test` | done | — |
| `[S28]`| Part B `build_expected_error_file` — merge per-sample `.qual` into project-wide file | `tests/processes/part_b/build_expected_error_file.nf.test` | done | — |
| `[S29]`| Part B `build_distribution_file` — sequence ↔ sample mapping from `.fas` headers | `tests/processes/part_b/build_distribution_file.nf.test` | done | — |
| `[S30]`| Part B `list_all_cluster_seeds_of_size_greater_than_2` — concatenate per-sample `.stats` | `tests/processes/part_b/list_all_cluster_seeds_of_size_greater_than_2.nf.test` | done | — |
| `[S31]`| Part B `global_dereplication` — vsearch --derep_fulllength across every input `.fas` | `tests/processes/part_b/global_dereplication.nf.test` | done | — |
| `[S32]`| Part B `global_clustering` — swarm `--fastidious` on the global fasta      | `tests/processes/part_b/global_clustering.nf.test` | done | — |
| `[S33]`| Part B `fake_taxonomic_assignment` — placeholder taxonomy TSV              | `tests/processes/part_b/fake_taxonomic_assignment.nf.test` | done | — |
| `[S34]`| Part B `chimera_detection` — vsearch uchime_denovo on representatives     | `tests/processes/part_b/chimera_detection.nf.test` | done | — |
| `[S35]`| Part B `build_occurrence_table` — merge swarm/uchime/quality/stampa/distribution into the filtered occurrence table | `tests/python/test_build_filtered_contingency_table.py`, `tests/processes/part_b/build_occurrence_table.nf.test` | done | — |
| `[S36]`| Part B `fake_taxonomic_assignment2` — placeholder taxonomy TSV for cleaved representatives | `tests/processes/part_b/fake_taxonomic_assignment2.nf.test` | done | — |
| `[S37]`| Part B `chimera_detection_post_cleave` — uchime_denovo on pre-cleave + cleaved reps with dynamic --minsize | `tests/processes/part_b/chimera_detection_post_cleave.nf.test` | done | — |
| `[S38]`| Part B `search_for_terminal_gaps` — vsearch self-cluster (id=1.0) for sub/super-string OTUs | `tests/processes/part_b/search_for_terminal_gaps.nf.test` | done | — |
| `[S39]`| Part B `merge_substring_otus` — merge pupil OTUs onto masters; sort; assert reads conserved | `tests/python/test_merge_substring_otus.py`, `tests/processes/part_b/merge_substring_otus.nf.test` | done | — |
| `[S40]`| Part B `extract_otu_fasta` / `extract_mumu_fasta` — FASTA from OTU table cols (4, 2, 10); the latter skips `$2==0` | `tests/processes/part_b/extract_otu_fasta.nf.test`, `tests/processes/part_b/extract_mumu_fasta.nf.test` | done | — |
| `[S41]`| Part B `trim_metadata_for_mumu` — keep amplicon + sample cols (cut -f 4,14-)         | `tests/processes/part_b/trim_metadata_for_mumu.nf.test` | done | — |
| `[S42]`| Part B `find_similar_sequences` — vsearch --usearch_global self-search; strip ;size= | `tests/processes/part_b/find_similar_sequences.nf.test` | done | — |
| `[S43]`| Part B `run_mumu` — mumu binary (>=1.1.1) post-clustering filter; `--minimum_relative_cooccurrence` coupled to cleaving `--percentage` (`[S22]`) as `1 - percentage` (default cleaving `0.05` → `0.95`) | `tests/processes/part_b/run_mumu.nf.test` | done | — |
| `[S44]`| Part B `rebuild_post_mumu_table` — splice old metadata onto mumu rows; renumber      | `tests/python/test_rebuild_table_after_mumu.py`, `tests/processes/part_b/rebuild_post_mumu_table.nf.test` | done | — |
| `[S45]`| Part B publishes six step-level log files under `logs/part_b/<basename>_<step>.log` | `tests/main.nf.test`, `tests/processes/part_b/global_dereplication.nf.test`, `tests/processes/part_b/global_clustering.nf.test`, `tests/processes/part_b/chimera_detection_post_cleave.nf.test`, `tests/processes/part_b/cleaving.nf.test`, `tests/processes/part_b/merge_substring_otus.nf.test`, `tests/processes/part_b/run_mumu.nf.test` | done | D15, D16 |
| `[S46]`| Part B publishes final occurrence table as `<basename>_table.tsv`                    | `tests/main.nf.test`, `tests/processes/part_b/rebuild_post_mumu_table.nf.test` | done   | — |
| `[S102]`| Optional post-mumu re-clustering gated on `--recluster_id` (master switch); `--recluster_iddef` inert without it | `tests/main.nf.test`, `tests/python/test_schema_params_sync.py` | done | D20 |
| `[S103]`| Part B `recluster_search` — `vsearch --cluster_size` on the post-mumu FASTA → `^H` connexions | `tests/processes/part_b/recluster_search.nf.test` | done | D20 |
| `[S104]`| Part B `recluster_merge` — `bin/recluster_otu_table.py` folds members onto centroids; renumber `1..N`; `cloud` stays `NA`; reads conserved | `tests/python/test_recluster_otu_table.py`, `tests/processes/part_b/recluster_merge.nf.test` | done | D20 |
| `[S105]`| Recluster replaces the emitted table + fasta, feeds Part C, logs `_reclustering.log`; OFF → byte-identical | `tests/main.nf.test` | done | D20 |
| `[S47]`| Part C requires `--reference_dataset` (no default)                                   | `tests/main.nf.test`                            | done    | —          |
| `[S48]`| Part C accepts an occurrence table (`--occurrence_table`, join → `_table_assigned.tsv`) or a representatives fasta (`--representatives_fasta`, no join → standalone `_taxonomy_<method>.tsv`); the two are mutually exclusive and fasta input rejects `--majority_assignment` | `tests/processes/part_c/extract_fasta_sequences_from_occurrence_table.nf.test`, `tests/subworkflows/part_c.nf.test`, `tests/main.nf.test` | done   | —        |
| `[S49]`| Part C stampa primary path — splitFasta + per-chunk `vsearch --usearch_global` + `bin/stampa_merge.py` + gather (`collectFile`) + `sort_taxonomy` (`LC_ALL=C sort -k2,2nr -k1,1d`); published `<basename>_taxonomy_stampa.tsv` carries a header row; per-chunk `vsearch.log` gathered into a published `logs/part_c/<basename>_taxonomy.log` | `tests/processes/part_c/assign_taxonomy_stampa.nf.test`, `tests/processes/part_c/sort_taxonomy.nf.test`, `tests/subworkflows/part_c.nf.test`, `tests/main.nf.test`, `tests/python/test_stampa_merge.py` | done | — |
| `[S50]`| Part C sintax shadow path — `part_C_shadow` runs `vsearch --sintax` on `<basename>_notmerged_table.tsv` and publishes `<basename>_notmerged_table_assigned.tsv` (its 4-col `_taxonomy_sintax.tsv` is never published) | `tests/processes/part_c/assign_taxonomy_sintax.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S91]`| `params.sintax_randseed` (default `0`, vsearch's random-data-source sentinel) threads through to `vsearch --sintax --randseed`; shared by the regular and shadow sintax paths; integer `>= 0` (`[S72]`) | `tests/processes/part_c/assign_taxonomy_sintax.nf.test` | done   | —          |
| `[S51]`| Part C `update_occurrence_table` — splice taxonomy back onto the occurrence table (tolerates a header row in the assignments TSV) | `tests/processes/part_c/update_occurrence_table.nf.test`, `tests/python/test_update_occurrence_table.py` | done   | D04        |
| `[S52]`| **Retired** (was: Part A U/u → T/t normalisation; dropped with the A-padding redesign — see `[S04]`, `[S63]`) | —                                | retired | —          |
| `[S53]`| **Retired** (was: `discover_fasta.py` rejects U/u; dropped with the A-padding redesign — see `[S04]`, `[S63]`) | —                                | retired | —          |
| `[S54]`| every vsearch fastq-emitting module preserves the canonical 4-line layout            | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/strip_reads.nf.test`, `tests/processes/part_a/join_notmerged.nf.test` | done | — |
| `[S55]`| every vsearch fasta-emitting module preserves the single-line-sequence layout        | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`, `tests/processes/part_a/dereplicate_fasta.nf.test`, `tests/processes/part_b/global_dereplication.nf.test` | done | — |
| `[S56]`| shadow Part B workflow — runs Part B as-is on A-padded shadow inputs; publishes `_notmerged` artefacts | `tests/main.nf.test` | done   | —          |
| `[S106]`| fastq-reading vsearch processes accept the full representable quality range (`--fastq_qmax = 126 - offset`, and the same as `--fastq_qmaxout` on merge) under either encoding, for PacBio HiFi-grade input | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/strip_reads.nf.test`, `tests/processes/part_a/filter_and_convert_to_fasta.nf.test` | done | — |
| `[S107]`| Part B publishes a per-sample read/cluster tracking summary `<basename>_read_counts.tsv` to `logs/part_b/` (reads_in/reads_kept + clusters per curation stage, empty samples as zero rows, Total row, optional recluster column; shadow path mirrored) | `tests/python/test_build_part_b_read_counts.py`, `tests/processes/part_b/summarize_part_b_read_counts.nf.test` | done | — |
| `[S57]`| `--help` prints a usage block describing all modes/params and exits without running any process | `tests/main.nf.test` | done | — |
| `[S58]`| `params.publish_mode` threads through every `publishDir` directive; invalid values abort at startup | `tests/main.nf.test` | done   | —          |
| `[S59]`| `occurrence_table/` whitelist: Part B table + post-mumu fasta; Part C `_table_assigned`, `_taxonomy_<method>` (sintax standalone only on explicit `--taxonomy_method=sintax`, never shadow), `_taxonomy_stampa_majority` (logs → `logs/part_b/` + `logs/part_c/`) | `tests/main.nf.test`                            | done    | D15, D16   |
| `[S60]`| path-typed params normalise leading `~` / `~user` at workflow startup (reference_dataset, occurrence_table, fastq_folder, fasta_folder, results_folder) | `tests/functions/normalize_path.nf.test`, `tests/main.nf.test` | done | — |
| `[S61]`| `params.taxonomy_method` validated at startup: accepts `stampa` (default) and `sintax`; invalid value aborts before any process runs; regular sintax publishes a 4-col `<basename>_taxonomy_sintax.tsv` (header `query/taxonomy/strand/cutoff_taxonomy`) | `tests/main.nf.test` | done | — |
| `[S62]`| Part C standalone probes for `<basename>_notmerged_table.tsv` sibling; runs `part_C_shadow` iff present, no-op otherwise | `tests/main.nf.test` | done | — |
| `[S63]`| `params.join_padding_length` (default 8): A-padding length used by `vsearch --fastq_join` in the shadow Part A path | `tests/processes/part_a/join_notmerged.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S64]`| `params.reference_dataset_sintax`: optional sintax-formatted reference; gates `part_C_shadow`; required by regular Part C when `--taxonomy_method=sintax` | `tests/main.nf.test` | done   | —          |
| `[S65]`| `params.hash_function` (default `sha1`, accepts `md5`): selects vsearch `--relabel_*`; `.qual` dedup width derived from the hash; invalid value aborts at startup | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`, `tests/processes/part_a/extract_expected_error_values.nf.test`, `tests/processes/part_b/build_expected_error_file.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S66]`| `params.majority_assignment` (opt-in, default `false`): regular Part C's final majority-rule step; publishes `<basename>_taxonomy_stampa_majority.tsv`; requires `--taxonomy_method=stampa` (sintax combo aborts at startup) | `tests/python/test_majority_assignment.py`, `tests/processes/part_c/compute_majority_assignment.nf.test`, `tests/main.nf.test` | done | — |
| `[S67]`| reject a `--fastq_pattern` whose `{<r1>,<r2>}` brace token has identical sides (e.g. `{1,1}`); equal sides would derive an R2 name equal to R1 and pair a file with itself | `tests/python/test_discover_fastq.py` | done | — |
| `[S96]`| cap `--fastq_pattern` at 3 `*` wildcards (ReDoS guard); each `*` is a greedy `.*`, an unbounded count enables catastrophic backtracking on scanned filenames | `tests/python/test_fastq_pattern_redos.py` | done | — |
| `[S68]`| every run records external-tool + Python versions into `pipeline_info/software_versions.yml` (missing tool → `n/a`); published on every entry point, incl. Part A-only | `tests/python/test_collect_versions.py`, `tests/processes/dump_software_versions.nf.test`, `tests/main.nf.test` | done | — |
| `[S69]`| exact, agreeing version pins in `environment.yml` and CI, incl. `bioconda::mumu` | `tests/python/test_reproducible_pins.py` | done | — |
| `[S77]`| `manifest.version` (`nextflow.config`) and `version` (`CITATION.cff`) must agree | `tests/python/test_version_sync.py` | done | — |
| `[S80]`| repo metadata consistent: `manifest.homePage` ↔ `CITATION.cff` url same `owner/repo`; `manifest.defaultBranch` set and schema `$id` points at it; `CHANGELOG.md` documents `manifest.version` | `tests/python/test_metadata_consistency.py` | done | — |
| `[S81]`| every piped process script runs under `set -euo pipefail` (left-pipe failure fails the task); documented empty-result pipes (`grep` no-match) stay `\|\| true`-guarded; SIGPIPE-safe min in post-cleave | `tests/processes/part_b/chimera_detection.nf.test`, `tests/processes/part_b/fake_taxonomic_assignment.nf.test`, `tests/processes/part_b/chimera_detection_post_cleave.nf.test` | done | — |
| `[S82]`| `cleanup` defaults to `false` (retain `work/` so `-resume` works across invocations; auto-clean is opt-in via `-c`) | `tests/check-cleanup-default.sh` (config resolution) | done | — |
| `[S83]`| air-gapped container path: build-once-then-offline (Wave cache) or a site-supplied `process.container` via `-c` (composes with `-profile slurm`, no Wave); documented in README + `conf/site.config.example` | `tests/check-offline-container.sh` (config resolution); offline execution = manual cluster smoke test | done | — |
| `[S84]`| every run writes Nextflow execution reports (timeline/report/trace/dag) to `<outdir>/pipeline_info/` for per-step resource tuning; demo profile re-points them to `demo_results/` | `tests/check-execution-reports.sh` (config resolution); generation = manual demo / cluster smoke run | done | — |
| `[S85]`| every tool-invoking process declares a `stub:`; `-profile demo -stub-run` runs the full Part A→B→C topology with no tool installed (discovery processes exempt) | `tests/check-stub-run.sh` (static gate + tool-free stub-run) | done | — |
| `[S86]`| Part A runs publish a per-sample read-count summary `<basename>_read_counts.tsv` to `logs/part_a/` (reads/assembled/F/R/passing from the merging + trimming logs, zeros for absent, Total row) | `tests/bin/build_read_counts.bats`, `tests/main.nf.test` | done | D16 |
| `[S78]`| `--recover_unmerged` (default false) gates the whole shadow path ([S04]); off by default produces no `_notmerged` artefacts / no `part_B_shadow` | `tests/main.nf.test` | done | — |
| `[S79]`| slurm runs warn at startup when `--dataset_size_gb` / `--reference_size_gb` are unset (fixed-fallback memory may OOM large runs); silent off-slurm | `tests/functions/resource_size_warnings.nf.test` | done | — |
| `[S70]`| `--input` samplesheet (fastq / fasta profiles): structural validation in `bin/parse_samplesheet.py`; folder-scan fallback, mutually exclusive | `tests/python/test_parse_samplesheet.py`, `tests/main.nf.test` | done | — |
| `[S95]`| reject any samplesheet cell containing a raw TAB / newline / CR (would corrupt the normalized TSV); aborts naming row + column | `tests/python/test_samplesheet_delimiters.py` | done | — |
| `[S71]`| `--outdir` single output root (data → `per_sample`/`occurrence_table`, logs → `logs/part_{a,b,c}`, `pipeline_info`); `--results_folder` deprecated alias; `--fastq_folder` input-only | `tests/functions/effective_outdir.nf.test`, `tests/main.nf.test` | done | D15, D16 |
| `[S72]`| numeric params range-validated at startup (fastq_encoding, threads, percentage, chimera_minsize, stripright, iddef, stampa_chunk_size, stampa_maxrejects, stampa_id, sintax_cutoff, primer_error_rate, primer_overlap_fraction, fastq_minlen, min_cluster_size, similar_id, similar_query_cov, similar_maxhits, max_ee, min_abundance, min_spread); out-of-range aborts naming the param. Declared in `nextflow_schema.json`, enforced by nf-schema `validateParameters()` | `tests/main.nf.test`, `tests/python/test_schema_params_sync.py` | done | — |
| `[S73]`| reference dataset format sniffed at startup: `--reference_dataset` must be stampa-shaped (`>id <lineage>`), `--reference_dataset_sintax` sintax-shaped (`>id;tax=...;`); plain + gzip read, bz2 skipped; mismatch aborts naming the flag | `tests/functions/check_reference_format.nf.test`, `tests/main.nf.test` | done | — |
| `[S74]`| primers validated at startup when trimming runs: IUPAC codes (A C G T U R Y S W K M B D H V N) + I, any case, >=3 nt; malformed aborts naming the param; quoted shell interpolation into cutadapt | `tests/functions/check_primer_format.nf.test`, `tests/main.nf.test` | done | — |
| `[S94]`| startup header sniffs ([S70]/[S73]) read the first line length-bounded (shared `read_bounded_line()`); guards against a huge first line / gzip decompression bomb OOM | `tests/functions/read_bounded_line.nf.test`, `tests/functions/check_reference_format.nf.test` | done | — |
| `[S75]`| site config via `-c` template (`conf/site.config.example`); `--slurm_clusterOptions` passthrough combined with `--account` | `tests/check-site-config.sh` (config resolution); submission = manual cluster smoke test | done | — |
| `[S76]`| self-contained `demo` profile (committed `assets/demo/` dataset) runs Part A→B→C with no flags; `-profile demo,<engine>` composes | `tests/check-demo-profile.sh` (assets + config resolution); run covered by Part A→B→C `tests/main.nf.test` | done | — |
| `[S87]`| vendored institutional cluster profiles (`conf/clusters/<name>.config`, shipped: abims/genotoul/ifb_core/meso/saga) each imply slurm via `conf/slurm.config`, clamp `resourceLimits` to the largest node, route partition/account, compose with an engine; vendored not fetched ([S83]); genotoul defaults job arrays on (`process.array=50`) + module-loaded container engine | `tests/check-cluster-profiles.sh` (profile wiring); submission = manual cluster smoke test | done | — |
| `[S92]`| a cluster profile may require `--slurm_account` via `params.require_slurm_account=true` (default false); when on and the account is unset/blank the run aborts at startup naming the flag; the `abims` profile sets it true | `tests/functions/slurm_account_requirement_error.nf.test`, `tests/check-cluster-profiles.sh` (abims resolves the flag) | done | — |
| `[S88]`| `--primer_error_rate` (default 0.1) → cutadapt `--error-rate` in both `trim_primers` passes | `tests/processes/part_a/trim_primers.nf.test` | done | — |
| `[S89]`| `--primer_overlap_fraction` (default 0.667) → cutadapt `--overlap` = int(primer_length × fraction) in `trim_primers` | `tests/processes/part_a/trim_primers.nf.test` | done | — |
| `[S90]`| `--fastq_minlen` (default 32) → vsearch `--fastq_minlen` in `filter_and_convert_to_fasta` | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test` | done | — |


## Fetch — download from ENA/SRA

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S97]`| `--accession` selects standalone fetch mode (param dispatch, not `-entry`); single or comma-separated list; sixth mutually-exclusive selector ([S02]); runs no Part A/B/C | `tests/main.nf.test` | done | D19 |
| `[S98]`| `--accession` accepts only bioproject (`^PRJ(E\|D\|N)[A-Z][0-9]+$`) and study (`^(E\|D\|S)RP[0-9]{6,}$`); anything else aborts at startup | `tests/functions/check_accession_format.nf.test`, `tests/main.nf.test` | done | D19 |
| `[S99]`| resolve stage maps each accession → its run accessions before download (one task per accession) | `tests/processes/fetch/resolve_runs.nf.test`, `tests/main.nf.test` | done | D19 |
| `[S100]`| fastq published under a per-accession subfolder `<outdir>/<accession>/` | `tests/processes/fetch/download_run.nf.test`, `tests/main.nf.test` | done | D19 |
| `[S101]`| per-run download via `fastq-dl=4.0.1` (per-process `conda`, not in `environment.yml`), `--provider ena`; failed run fails only its task | `tests/processes/fetch/download_run.nf.test`, `tests/python/test_reproducible_pins.py`, `tests/main.nf.test` | done | D19 |


## Per-process tests

| Process in `main.nf`            | Test file                                              | Covers       | Status |
|---------------------------------|--------------------------------------------------------|--------------|--------|
| `merge_fastq_pairs`             | `tests/processes/part_a/merge_fastq_pairs.nf.test`            | S01, S03, S04, S54, S106 | done   |
| `trim_primers`                  | `tests/processes/part_a/trim_primers.nf.test`                 | S01, S19, S88, S89 | done   |
| `filter_and_convert_to_fasta`   | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`  | S01, S55, S65, S90, S106 | done  |
| `extract_expected_error_values` | `tests/processes/part_a/extract_expected_error_values.nf.test`| S01, S65     | done   |
| `dereplicate_fasta`             | `tests/processes/part_a/dereplicate_fasta.nf.test`            | S01, S19, S55| done   |
| `list_local_clusters`           | `tests/processes/part_a/list_local_clusters.nf.test`          | S17, S19     | done   |
| `join_notmerged`                | `tests/processes/part_a/join_notmerged.nf.test`               | S04, S19, S54, S63 | done   |
| `strip_reads`                   | `tests/processes/part_a/strip_reads.nf.test`                  | S24, S54, S106 | done   |
| `build_expected_error_file`     | `tests/processes/part_b/build_expected_error_file.nf.test`    | S28, S65     | done   |
| `build_distribution_file`       | `tests/processes/part_b/build_distribution_file.nf.test`      | S29          | done   |
| `list_all_cluster_seeds_of_size_greater_than_2` | `tests/processes/part_b/list_all_cluster_seeds_of_size_greater_than_2.nf.test` | S30 | done |
| `global_dereplication`          | `tests/processes/part_b/global_dereplication.nf.test`         | S31, S55     | done   |
| `global_clustering`             | `tests/processes/part_b/global_clustering.nf.test`            | S32, S45     | done   |
| `fake_taxonomic_assignment`     | `tests/processes/part_b/fake_taxonomic_assignment.nf.test`    | S33          | done   |
| `chimera_detection`             | `tests/processes/part_b/chimera_detection.nf.test`            | S34          | done   |
| `cleaving`                      | `tests/processes/part_b/cleaving.nf.test`                     | S22          | done   |
| `build_occurrence_table`        | `tests/processes/part_b/build_occurrence_table.nf.test`       | S35          | done   |
| `fake_taxonomic_assignment2`    | `tests/processes/part_b/fake_taxonomic_assignment2.nf.test`   | S36          | done   |
| `chimera_detection_post_cleave`            | `tests/processes/part_b/chimera_detection_post_cleave.nf.test`           | S37          | done   |
| `search_for_terminal_gaps`      | `tests/processes/part_b/search_for_terminal_gaps.nf.test`     | S38          | done   |
| `merge_substring_otus`          | `tests/processes/part_b/merge_substring_otus.nf.test`         | S39          | done   |
| `extract_otu_fasta`             | `tests/processes/part_b/extract_otu_fasta.nf.test`            | S40          | done   |
| `extract_mumu_fasta`            | `tests/processes/part_b/extract_mumu_fasta.nf.test`           | S40          | done   |
| `trim_metadata_for_mumu`        | `tests/processes/part_b/trim_metadata_for_mumu.nf.test`       | S41          | done   |
| `find_similar_sequences`        | `tests/processes/part_b/find_similar_sequences.nf.test`       | S42          | done   |
| `run_mumu`                      | `tests/processes/part_b/run_mumu.nf.test`                     | S43          | done   |
| `rebuild_post_mumu_table`       | `tests/processes/part_b/rebuild_post_mumu_table.nf.test`      | S44, S46     | done   |
| `recluster_search`              | `tests/processes/part_b/recluster_search.nf.test`            | S103         | done   |
| `recluster_merge`               | `tests/processes/part_b/recluster_merge.nf.test`             | S104, S105   | done   |
| `summarize_part_b_read_counts`  | `tests/processes/part_b/summarize_part_b_read_counts.nf.test` | S107        | done   |
| `extract_fasta_sequences_from_occurrence_table` | `tests/processes/part_c/extract_fasta_sequences_from_occurrence_table.nf.test` | S48 | done   |
| `assign_taxonomy_stampa`        | `tests/processes/part_c/assign_taxonomy_stampa.nf.test`       | S49          | done   |
| `sort_taxonomy`                 | `tests/processes/part_c/sort_taxonomy.nf.test`               | S49          | done   |
| `assign_taxonomy_sintax`        | `tests/processes/part_c/assign_taxonomy_sintax.nf.test`       | S50          | done   |
| `update_occurrence_table`       | `tests/processes/part_c/update_occurrence_table.nf.test`      | S51          | done   |
| `compute_majority_assignment`   | `tests/processes/part_c/compute_majority_assignment.nf.test`  | S66          | done   |
| `dump_software_versions`        | `tests/processes/dump_software_versions.nf.test`              | S68          | done   |
| `resolve_runs`                  | `tests/processes/fetch/resolve_runs.nf.test`                 | S99          | done   |
| `download_run`                  | `tests/processes/fetch/download_run.nf.test`                 | S100, S101   | done   |
