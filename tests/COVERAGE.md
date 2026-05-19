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
| `[S01]`| three-part workflow (fastq→fasta, fasta→occurrence, taxonomic assignment)  | `tests/main.nf.test`                            | red     | —          |
| `[S02]`| each part can be run separately or all at once                             | `tests/main.nf.test`                            | TODO    | —          |
| `[S03]`| paired-end or single-end, compressed (gz/bz2) or uncompressed input        | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/main.nf.test` | red | — |
| `[S04]`| unmerged paired reads → shadow pipeline (N-join, N→A mask before swarm)    | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/join_notmerged.nf.test`, `tests/processes/part_a/mask_ns_for_swarm.nf.test`, `tests/main.nf.test` | done | — |
| `[S05]`| unmerged clusters appear in occurrence table with per-sample marker        | —                                               | blocked | D01, D02   |


## Workflow requirements

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S06]`| read config file *or* command-line parameters                              | `tests/main.nf.test`                            | TODO    | —          |
| `[S07]`| runs locally or on HPC (slurm)                                             | —                                               | n/a     | —          |
| `[S08]`| runs local *or* containerized applications                                 | —                                               | n/a     | —          |
| `[S09]`| empty input samples travel through and appear in the occurrence table      | `tests/main.nf.test`                            | red     | —          |
| `[S10]`| accept a directory or a list of directories (absolute or relative)         | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S11]`| auto-discover fastq files using common name patterns                       | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S12]`| auto-deduce sample names from fastq file names                             | `tests/python/test_discover_fastq.py`, `tests/bin/discover_fastq.bats`, `tests/main.nf.test` | done | —          |
| `[S13]`| abort if two or more samples share a derived ID; list duplicate file paths | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py`, `tests/main.nf.test` | done | — |
| `[S14]`| collision policy for same-named samples: refuse (sample IDs must be unique) | `tests/python/test_discover_fastq.py`, `tests/python/test_discover_fasta.py` | done    | —          |
| `[S15]`| export single occurrence table *or* two-part (long + metadata) table       | `tests/main.nf.test`                            | TODO    | —          |
| `[S16]`| expect demultiplexed fastq files                                           | —                                               | n/a     | —          |
| `[S17]`| per-cluster minimum-read threshold (> 2 reads)                             | `tests/processes/part_a/list_local_clusters.nf.test`   | red     | —          |
| `[S18]`| required params (forward/reverse_primer, fastq_folder) must be supplied    | `tests/main.nf.test`                            | done    | —          |
| `[S19]`| Part A steps emit per-sample `<sampleId>_<step>.log` files                 | `tests/processes/part_a/merge_fastq_pairs.nf.test`, `tests/processes/part_a/trim_primers.nf.test`, `tests/processes/part_a/dereplicate_fasta.nf.test`, `tests/processes/part_a/list_local_clusters.nf.test` | done | — |
| `[S20]`| `--no_trimming` toggle skips primer trimming; mutually exclusive w/ primers| `tests/main.nf.test`                            | red     | —          |
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


## Per-process tests

| Process in `main.nf`            | Test file                                              | Covers       | Status |
|---------------------------------|--------------------------------------------------------|--------------|--------|
| `merge_fastq_pairs`             | `tests/processes/part_a/merge_fastq_pairs.nf.test`            | S01, S03, S04| red    |
| `trim_primers`                  | `tests/processes/part_a/trim_primers.nf.test`                 | S01          | red    |
| `filter_and_convert_to_fasta`   | `tests/processes/part_a/filter_and_convert_to_fasta.nf.test`  | S01          | red    |
| `extract_expected_error_values` | `tests/processes/part_a/extract_expected_error_values.nf.test`| S01          | red    |
| `dereplicate_fasta`             | `tests/processes/part_a/dereplicate_fasta.nf.test`            | S01          | red    |
| `list_local_clusters`           | `tests/processes/part_a/list_local_clusters.nf.test`          | S17          | red    |
| `join_notmerged`                | `tests/processes/part_a/join_notmerged.nf.test`               | S04, S19     | red    |
| `mask_ns_for_swarm`             | `tests/processes/part_a/mask_ns_for_swarm.nf.test`            | S04          | red    |
| `strip_reads`                   | `tests/processes/part_a/strip_reads.nf.test`                  | S24          | red    |
| `build_expected_error_file`     | `tests/processes/part_b/build_expected_error_file.nf.test`    | S28          | red    |
| `build_distribution_file`       | `tests/processes/part_b/build_distribution_file.nf.test`      | S29          | red    |
| `list_all_cluster_seeds_of_size_greater_than_2` | `tests/processes/part_b/list_all_cluster_seeds_of_size_greater_than_2.nf.test` | S30 | red |
| `global_dereplication`          | `tests/processes/part_b/global_dereplication.nf.test`         | S25, S31     | red    |
| `global_clustering`             | `tests/processes/part_b/global_clustering.nf.test`            | S32          | red    |
