# Test fixtures

These fixtures are **generated**, not committed — run
`bash tests/data/generate.sh` from the repo root to (re)create them.
The `tests/data/*.fastq*` paths are listed in `.gitignore`.

Fixtures are deliberately tiny: a handful of reads is enough to
exercise each behaviour. Real-world correctness of `vsearch`,
`cutadapt`, and `swarm` is the responsibility of those projects, not
this workflow.

## Inventory

| Fixture                                                   | Used by spec bullet         | What it tests                                              |
|-----------------------------------------------------------|-----------------------------|------------------------------------------------------------|
| `paired_merge_ok_1.fastq.gz` / `paired_merge_ok_2.fastq.gz` | #1, #3                    | the happy path: paired reads that overlap and merge        |
| `paired_merge_fail_1.fastq.gz` / `paired_merge_fail_2.fastq.gz` | #5                    | reads that *cannot* merge (parallel pipeline)              |
| `single_end.fastq.gz`                                     | #3                          | single-end input                                           |
| `uncompressed_1.fastq` / `uncompressed_2.fastq`           | #4                          | uncompressed input                                         |
| `empty_1.fastq.gz` / `empty_2.fastq.gz`                   | #10                         | empty input samples must travel through                    |

Numbers refer to rows in [`../COVERAGE.md`](../COVERAGE.md).

## Sequence design

Reads are synthetic and deterministic. Each read contains the forward
primer (`CCAGCASCYGCGGTAATTCC`) at the 5' end and the reverse-complement
of the reverse primer (`ACTTTCGTTCTTGATYRA`) at the 3' end of the
amplicon. Quality strings are constant `I` (Phred 40) — quality-aware
behaviour should be tested with a dedicated fixture.
