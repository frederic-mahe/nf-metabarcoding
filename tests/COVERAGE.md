# Spec coverage

Each row maps a bullet from [`../SPECIFICATIONS.md`](../SPECIFICATIONS.md)
to the test(s) that cover it. **Status legend:**

- `done`    — test exists and passes
- `red`     — test exists and fails (the code is not yet there — this is the TDD red phase)
- `TODO`    — no test yet
- `n/a`     — meta or not testable
- `blocked` — spec needs clarification before a test can be written

When a spec changes, update this table in the same commit.

## Workflow structure

| # | Spec bullet                                                                 | Test file                                       | Status  |
|---|-----------------------------------------------------------------------------|-------------------------------------------------|---------|
| 1 | three-part workflow (fastq→fasta, fasta→occurrence, taxonomic assignment)   | `tests/main.nf.test`                            | red     |
| 2 | each part can be run separately or all at once                              | `tests/main.nf.test`                            | TODO    |
| 3 | paired-end or single-end input                                              | `tests/processes/merge_fastq_pairs.nf.test`     | red     |
| 4 | compressed (gz/bz2) or uncompressed input                                   | `tests/main.nf.test`                            | TODO    |
| 5 | unmerged paired reads → parallel pipeline (N-join, N↔A round-trip)          | `tests/processes/merge_fastq_pairs.nf.test`     | blocked |
| 6 | unmerged clusters appear as `sampleID_partial` (placeholder)                | —                                               | blocked |

## Workflow requirements

| #  | Spec bullet                                                          | Test file                                       | Status  |
|----|----------------------------------------------------------------------|-------------------------------------------------|---------|
| 7  | read config file *or* command-line parameters                        | `tests/main.nf.test`                            | TODO    |
| 8  | run locally or on HPC (slurm)                                        | —                                               | n/a     |
| 9  | run local *or* containerized applications                            | —                                               | n/a     |
| 10 | empty input samples travel through and appear in the occurrence table| `tests/main.nf.test`                            | red     |
| 11 | accept a directory or a list of directories (absolute or relative)   | `tests/main.nf.test`                            | TODO    |
| 12 | auto-discover fastq files using common name patterns                 | `tests/main.nf.test`                            | TODO    |
| 13 | auto-deduce sample names from fastq file names                       | `tests/main.nf.test`                            | TODO    |
| 14 | warn if two or more samples share a name                             | `tests/main.nf.test`                            | TODO    |
| 15 | merge samples with the same name                                     | —                                               | blocked |
| 16 | export single occurrence table *or* two-part (long + metadata) table | `tests/main.nf.test`                            | TODO    |
| 17 | expect demultiplexed fastq files                                     | —                                               | n/a     |

## Per-process tests

| Process in `main.nf`            | Test file                                              | Status |
|---------------------------------|--------------------------------------------------------|--------|
| `merge_fastq_pairs`             | `tests/processes/merge_fastq_pairs.nf.test`            | red    |
| `trim_primers`                  | `tests/processes/trim_primers.nf.test`                 | red    |
| `convert_fastq_to_fasta`        | `tests/processes/convert_fastq_to_fasta.nf.test`       | red    |
| `extract_expected_error_values` | `tests/processes/extract_expected_error_values.nf.test`| red    |
| `dereplicate_fasta`             | `tests/processes/dereplicate_fasta.nf.test`            | red    |
| `list_local_clusters`           | `tests/processes/list_local_clusters.nf.test`          | red    |

## Blocked items needing spec clarification

- **#5/#6** — N↔A conversion: which positions get rewritten, how is the
  round-trip tracked, what marker replaces `_partial`? (SPECIFICATIONS lines 27–31)
- **#15** — Should same-named samples be merged, refused, or kept distinct
  with a suffix? (SPECIFICATIONS line 49)
