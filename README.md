# fred-metabarcoding-pipeline-nextflow

Fred's metabarcoding pipeline with Nextflow

This is a first attempt at converting to
[Nextflow](https://www.nextflow.io/) the pipeline I use for my own
metabarcoding projects.

Current status:
- [x] draft of the first section of the pipeline (process individual samples),
- [ ] collect info from log files,
- [ ] add checks:
  - [ ] user-set parameters (empty strings, unrealistic values?),
  - [ ] dependencies and minimal versions (vsearch, cutadapt, swarm, bash >= 4),
- [ ] adapt to work with slurm,
- [ ] automatically deduce the fastq file naming pattern,
- [ ] automatically deduce compression (gz, bz2) or the lack-of,
- [ ] extend to multiplexed datasets,
- [ ] draft of the second section of the pipeline (work at the study scale),
- [ ] deduce reference database from primers

List of common name pattern for paired-end fastq files:
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

List of dependencies that could trigger a complete or partial re-run:
- new versions of vsearch, cutadapt, and swarm,
- new versions of external python scripts,
- new version of the reference database,
- new set of fastq files
