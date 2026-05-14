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
| `[S03]`| paired-end or single-end, compressed (gz/bz2) or uncompressed input        | `tests/processes/merge_fastq_pairs.nf.test`, `tests/main.nf.test` | red | — |
| `[S04]`| unmerged paired reads → parallel pipeline (N-join, N↔A round-trip)         | `tests/processes/merge_fastq_pairs.nf.test`     | blocked | D01        |
| `[S05]`| unmerged clusters appear in occurrence table with per-sample marker        | —                                               | blocked | D01, D02   |


## Workflow requirements

| Spec   | Bullet                                                                     | Test file                                       | Status  | Blocked by |
|--------|----------------------------------------------------------------------------|-------------------------------------------------|---------|------------|
| `[S06]`| read config file *or* command-line parameters                              | `tests/main.nf.test`                            | TODO    | —          |
| `[S07]`| runs locally or on HPC (slurm)                                             | —                                               | n/a     | —          |
| `[S08]`| runs local *or* containerized applications                                 | —                                               | n/a     | —          |
| `[S09]`| empty input samples travel through and appear in the occurrence table      | `tests/main.nf.test`                            | red     | —          |
| `[S10]`| accept a directory or a list of directories (absolute or relative)         | `tests/main.nf.test`                            | TODO    | —          |
| `[S11]`| auto-discover fastq files using common name patterns                       | `tests/main.nf.test`                            | TODO    | —          |
| `[S12]`| auto-deduce sample names from fastq file names                             | `tests/main.nf.test`                            | TODO    | —          |
| `[S13]`| warn if two or more samples share a name                                   | `tests/main.nf.test`                            | TODO    | —          |
| `[S14]`| collision policy for same-named samples                                    | —                                               | blocked | D03        |
| `[S15]`| export single occurrence table *or* two-part (long + metadata) table       | `tests/main.nf.test`                            | TODO    | —          |
| `[S16]`| expect demultiplexed fastq files                                           | —                                               | n/a     | —          |
| `[S17]`| per-cluster minimum-read threshold (> 2 reads)                             | `tests/processes/list_local_clusters.nf.test`   | red     | —          |


## Per-process tests

| Process in `main.nf`            | Test file                                              | Covers       | Status |
|---------------------------------|--------------------------------------------------------|--------------|--------|
| `merge_fastq_pairs`             | `tests/processes/merge_fastq_pairs.nf.test`            | S01, S03, S04| red    |
| `trim_primers`                  | `tests/processes/trim_primers.nf.test`                 | S01          | red    |
| `convert_fastq_to_fasta`        | `tests/processes/convert_fastq_to_fasta.nf.test`       | S01          | red    |
| `extract_expected_error_values` | `tests/processes/extract_expected_error_values.nf.test`| S01          | red    |
| `dereplicate_fasta`             | `tests/processes/dereplicate_fasta.nf.test`            | S01          | red    |
| `list_local_clusters`           | `tests/processes/list_local_clusters.nf.test`          | S17          | red    |
