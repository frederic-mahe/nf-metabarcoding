# Specifications

## tests

- test-driven development (write tests first, each specification must
  be covered by tests)
- tests are minimalistic input files (as small as possible, just
  enough to test a specification)
- pipeline relies on vsearch, cutadapt, and swarm. These are
  well-tested and do not need to be tested specifically


## workflow structure

- three parts:
    1) from fastq files to dereplicated fasta files (merge reads, trim
    primers, extract quality, local-clustering with swarm). See
    `../fred-metabarcoding-pipeline/` for details.
    2) from dereplicated fasta files to an occurrence table (vsearch,
    swarm, and python scripts
    3) taxonomic assignment: update the occurrence table
- each part can be run separately, or all at once
- fastq files can be paired-end, or single, compressed or not
- when processing paired-end fastq files, reads that can be merged are
  processed normaly. Reads that cannot be merged follow a parallel
  pipeline (reads are joined with Ns, Ns are converted to As when
  passed to swarm. The exact process of replacing Ns with As and back
  needs to be defined, to avoid the risk of replacing legit As). These
  appear in the occurrence table as sampleID_partial (`_partial` is a
  placeholder, maybe there is a better way to mark clusters of
  non-merged reads)



## workflow requirements

- read a config file or use command-line parameters
- run locally or on HPC with slurm
- empty input samples must travel through and appear in the occurrence
  table
- accept a directory, or a list of directories (absolute or relative
  paths)
- search listed directories automatically (see list of common fastq
  name patterns below)
- automatically deduce sample names from fastq file names
- warn if two or more samples have the same names
- merge samples with the same names?? (not clear what is the best
  option)
- allow users to export a single occurrence table, or a two-part table
  (occurrences in long-format, per-cluster metadata in another table)
- expect demultiplexed fastq files (demultiplexing could be dealt with
  by a subworkflow)

List of common name pattern for paired-end fastq files (see
`../fred-metabarcoding-pipeline/` for a more up-to-date version):
- `_L001_R{1,2}_001.fastq.gz` (MiSeq)
- `_L001_R{1,2}.fastq.gz`
- `_[1-9]_{1,2}.fastq.gz`
- `_[1-9]_{1,2}_.*.fastq.gz`
- `_L001_.*_R{1,2}.fastq.bz2` (variant of `_L001_R{1,2}.fastq.gz`)
- `_L005_R{1,2}.fastq.gz` (lane can be greater than 1!)
- `_L001_R{1,2}_002.fastq.bz2` (the last segment is supposed to always be 001!)
- `_R{1,2}.fastq.gz`
- `_{1,2}.fastq.gz`

- note: `fastq` might be `fq`, and there could be no compression
  (`.(fastq|fq)(.gz|.bz2)?`)

