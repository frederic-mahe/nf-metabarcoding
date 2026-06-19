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
| `[S02]`| each part can be run separately or all at once                             | `tests/main.nf.test`                            | done    | —          |
| `[S03]`| paired-end or single-end, compressed (gz/bz2) or uncompressed input        | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/main.nf.test` | done | — |
| `[S04]`| unmerged paired reads → shadow pipeline (A-padded join, no mask round-trip)| `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/join_notmerged.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S05]`| unmerged clusters appear in occurrence table with per-sample marker        | —                                               | blocked | D01, D02   |


## Workflow requirements

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S06]`| read config file *or* command-line parameters                              | —                                               | n/a     | Nextflow's intrinsic `params` indirection; required-param branch covered by `[S18]` |
| `[S07]`| runs locally or on HPC (slurm)                                             | manual smoke test; `slurm` profile in `nextflow.config` | n/a     | per-process resource tiers, dataset-scaled memory, retry escalation |
| `[S08]`| runs local *or* containerized applications (docker/podman/singularity/apptainer via Wave) | `tests/check-container-profiles.sh` (profile wiring); execution = manual cluster smoke test | done    | —          |
| `[S09]`| empty input samples travel through and appear in the occurrence table      | `tests/main.nf.test`, `tests/python/test_build_filtered_contingency_table.py` | done    | —          |
| `[S10]`| accept a directory or a list of directories (absolute or relative)         | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S11]`| auto-discover fastq files using common name patterns                       | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S12]`| auto-deduce sample names from fastq file names                             | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S13]`| abort if two or more samples share a derived ID; list duplicate file paths | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py`, `tests/main.nf.test` | done | — |
| `[S14]`| collision policy for same-named samples: refuse (sample IDs must be unique) | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py` | done    | —          |
| `[S15]`| export single occurrence table *or* two-part (long + metadata) table       | `tests/main.nf.test`                            | TODO    | —          |
| `[S16]`| expect demultiplexed fastq files                                           | —                                               | n/a     | —          |
| `[S17]`| per-cluster minimum-read threshold (> 2 reads)                             | `tests/processes/part_a/list_local_clusters.nf.test`, `tests/main.nf.test`   | done    | —          |
| `[S18]`| required params (forward/reverse_primer, fastq_folder) must be supplied    | `tests/main.nf.test`                            | done    | —          |
| `[S19]`| Part A steps emit per-sample `<sampleId>_<step>.log` files                 | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/trim_primers.nf.test`, `tests/processes/part_a/dereplicate_fasta.nf.test`, `tests/processes/part_a/list_local_clusters.nf.test` | done | — |
| `[S20]`| `--no_trimming` toggle skips primer trimming; mutually exclusive w/ primers| `tests/main.nf.test`                            | done    | —          |
| `[S21]`| unpaired fastq files skip the merging step                                 | `tests/main.nf.test`                            | done    | —          |
| `[S22]`| Part B re-cleaves global swarm clusters using per-sample sub-seed presence | `tests/python/test_cluster_cleaver.py`          | done    | —          |
| `[S23]`| `notmerged` reserved suffix — sample IDs ending in `notmerged` are rejected | `tests/python/test_reserved_keyword.py`, `tests/main.nf.test` | done | — |
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
| `[S43]`| Part B `run_mumu` — mumu binary (>=1.1.1) post-clustering filter                     | `tests/processes/part_b/run_mumu.nf.test` | done | — |
| `[S44]`| Part B `rebuild_post_mumu_table` — splice old metadata onto mumu rows; renumber      | `tests/python/test_rebuild_table_after_mumu.py`, `tests/processes/part_b/rebuild_post_mumu_table.nf.test` | done | — |
| `[S45]`| Part B publishes six step-level log files under `<basename>_<step>.log`              | `tests/main.nf.test`, `tests/processes/part_b/global_dereplication.nf.test`, `tests/processes/part_b/global_clustering.nf.test`, `tests/processes/part_b/chimera_detection_post_cleave.nf.test`, `tests/processes/part_b/cleaving.nf.test`, `tests/processes/part_b/merge_substring_otus.nf.test`, `tests/processes/part_b/run_mumu.nf.test` | done | — |
| `[S46]`| Part B publishes final occurrence table as `<basename>_table.tsv`                    | `tests/main.nf.test`, `tests/processes/part_b/rebuild_post_mumu_table.nf.test` | done   | — |
| `[S47]`| Part C requires `--reference_dataset` (no default)                                   | `tests/main.nf.test`                            | done    | —          |
| `[S48]`| Part C accepts either an occurrence table or a fasta file (fasta-input branch blocked)| `tests/processes/part_c/extract_fasta_sequences_from_occurrence_table.nf.test` | done   | D04        |
| `[S49]`| Part C stampa primary path — splitFasta + per-chunk `vsearch --usearch_global` + `bin/stampa_merge.py` + `collectFile(sort:)` | `tests/processes/part_c/assign_taxonomy_stampa.nf.test`, `tests/main.nf.test`, `tests/python/test_stampa_merge.py` | done | — |
| `[S50]`| Part C sintax shadow path — `part_C_shadow` runs `vsearch --sintax` on `<basename>_notmerged_table.tsv` and publishes `<basename>_notmerged_table_assigned.tsv` | `tests/processes/part_c/assign_taxonomy_sintax.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S51]`| Part C `update_occurrence_table` — splice taxonomy back onto the occurrence table   | `tests/processes/part_c/update_occurrence_table.nf.test` | done   | D04        |
| `[S52]`| **Retired** (was: Part A U/u → T/t normalisation; dropped with the A-padding redesign — see `[S04]`, `[S63]`) | —                                | retired | —          |
| `[S53]`| **Retired** (was: `discover_fasta.py` rejects U/u; dropped with the A-padding redesign — see `[S04]`, `[S63]`) | —                                | retired | —          |
| `[S54]`| every vsearch fastq-emitting module preserves the canonical 4-line layout            | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/strip_reads.nf.test`, `tests/processes/part_a/join_notmerged.nf.test` | done | — |
| `[S55]`| every vsearch fasta-emitting module preserves the single-line-sequence layout        | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`, `tests/processes/part_a/dereplicate_fasta.nf.test`, `tests/processes/part_b/global_dereplication.nf.test` | done | — |
| `[S56]`| shadow Part B workflow — runs Part B as-is on A-padded shadow inputs; publishes `_notmerged` artefacts | `tests/main.nf.test` | done   | —          |
| `[S57]`| `--help` prints a usage block describing all modes/params and exits without running any process | `tests/main.nf.test` | done | — |
| `[S58]`| `params.publish_mode` threads through every `publishDir` directive; invalid values abort at startup | `tests/main.nf.test` | done   | —          |
| `[S59]`| Part B `--results_folder` whitelist: table + post-mumu fasta + six step logs only                  | `tests/main.nf.test`                            | done    | —          |
| `[S60]`| path-typed params normalise leading `~` / `~user` at workflow startup (reference_dataset, occurrence_table, fastq_folder, fasta_folder, results_folder) | `tests/functions/normalize_path.nf.test`, `tests/main.nf.test` | done | — |
| `[S61]`| `params.taxonomy_method` validated at startup: accepts `stampa` (default) and `sintax`; invalid value aborts before any process runs | `tests/main.nf.test` | done | — |
| `[S62]`| Part C standalone probes for `<basename>_notmerged_table.tsv` sibling; runs `part_C_shadow` iff present, no-op otherwise | `tests/main.nf.test` | done | — |
| `[S63]`| `params.join_padding_length` (default 8): A-padding length used by `vsearch --fastq_join` in the shadow Part A path | `tests/processes/part_a/join_notmerged.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S64]`| `params.reference_dataset_sintax`: optional sintax-formatted reference; gates `part_C_shadow`; required by regular Part C when `--taxonomy_method=sintax` | `tests/main.nf.test` | done   | —          |
| `[S65]`| `params.hash_function` (default `sha1`, accepts `md5`): selects vsearch `--relabel_*`; `.qual` dedup width derived from the hash; invalid value aborts at startup | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`, `tests/processes/part_a/extract_expected_error_values.nf.test`, `tests/processes/part_b/build_expected_error_file.nf.test`, `tests/main.nf.test` | done   | —          |
| `[S66]`| `params.majority_assignment` (opt-in, default `false`): regular Part C's final majority-rule step; publishes `<basename>_table_assigned_majority.tsv`; requires `--taxonomy_method=stampa` (sintax combo aborts at startup) | `tests/python/test_majority_assignment.py`, `tests/processes/part_c/compute_majority_assignment.nf.test`, `tests/main.nf.test` | done | — |
| `[S67]`| reject a `--fastq_pattern` whose `{<r1>,<r2>}` brace token has identical sides (e.g. `{1,1}`); equal sides would derive an R2 name equal to R1 and pair a file with itself | `tests/python/test_discover_fastq.py` | done | — |
| `[S68]`| every run records external-tool + Python versions into `pipeline_info/software_versions.yml` (missing tool → `n/a`) | `tests/python/test_collect_versions.py`, `tests/processes/dump_software_versions.nf.test` | done | — |
| `[S69]`| exact, agreeing version pins in `environment.yml` and CI, incl. `bioconda::mumu` | `tests/python/test_reproducible_pins.py` | done | — |
| `[S70]`| `--input` samplesheet (fastq / fasta profiles): structural validation in `bin/parse_samplesheet.py`; folder-scan fallback, mutually exclusive | `tests/python/test_parse_samplesheet.py`, `tests/main.nf.test` | done | — |
| `[S71]`| `--outdir` single output root (`per_sample`/`occurrence_table`/`pipeline_info`); `--results_folder` deprecated alias; `--fastq_folder` input-only | `tests/functions/effective_outdir.nf.test`, `tests/main.nf.test` | done | — |
| `[S72]`| numeric params range-validated at startup (fastq_encoding, threads, percentage, chimera_minsize, stripright, iddef, stampa_chunk_size, stampa_maxrejects, stampa_id, sintax_cutoff); out-of-range aborts naming the param | `tests/functions/check_numeric_param.nf.test`, `tests/main.nf.test` | done | — |
| `[S73]`| reference dataset format sniffed at startup: `--reference_dataset` must be stampa-shaped (`>id <lineage>`), `--reference_dataset_sintax` sintax-shaped (`>id;tax=...;`); plain + gzip read, bz2 skipped; mismatch aborts naming the flag | `tests/functions/check_reference_format.nf.test`, `tests/main.nf.test` | done | — |
| `[S74]`| primers validated at startup when trimming runs: IUPAC codes (A C G T U R Y S W K M B D H V N) + I, any case, >=3 nt; malformed aborts naming the param; quoted shell interpolation into cutadapt | `tests/functions/check_primer_format.nf.test`, `tests/main.nf.test` | done | — |
| `[S75]`| site config via `-c` template (`conf/site.config.example`); `--slurm_clusterOptions` passthrough combined with `--account` | `tests/check-site-config.sh` (config resolution); submission = manual cluster smoke test | done | — |


## Per-process tests

| Process in `main.nf`            | Test file                                              | Covers       | Status |
|---------------------------------|--------------------------------------------------------|--------------|--------|
| `merge_fastq_pairs`             | `tests/processes/part_a/merge_fastq_pairs.nf.test`            | S01, S03, S04, S54 | done   |
| `trim_primers`                  | `tests/processes/part_a/trim_primers.nf.test`                 | S01, S19     | done   |
| `filter_and_convert_to_fasta`   | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`  | S01, S55, S65 | done  |
| `extract_expected_error_values` | `tests/processes/part_a/extract_expected_error_values.nf.test`| S01, S65     | done   |
| `dereplicate_fasta`             | `tests/processes/part_a/dereplicate_fasta.nf.test`            | S01, S19, S55| done   |
| `list_local_clusters`           | `tests/processes/part_a/list_local_clusters.nf.test`          | S17, S19     | done   |
| `join_notmerged`                | `tests/processes/part_a/join_notmerged.nf.test`               | S04, S19, S54, S63 | done   |
| `strip_reads`                   | `tests/processes/part_a/strip_reads.nf.test`                  | S24, S54     | done   |
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
| `extract_fasta_sequences_from_occurrence_table` | `tests/processes/part_c/extract_fasta_sequences_from_occurrence_table.nf.test` | S48 | done   |
| `assign_taxonomy_stampa`        | `tests/processes/part_c/assign_taxonomy_stampa.nf.test`       | S49          | done   |
| `assign_taxonomy_sintax`        | `tests/processes/part_c/assign_taxonomy_sintax.nf.test`       | S50          | done   |
| `update_occurrence_table`       | `tests/processes/part_c/update_occurrence_table.nf.test`      | S51          | done   |
| `compute_majority_assignment`   | `tests/processes/part_c/compute_majority_assignment.nf.test`  | S66          | done   |
| `dump_software_versions`        | `tests/processes/dump_software_versions.nf.test`              | S68          | done   |
