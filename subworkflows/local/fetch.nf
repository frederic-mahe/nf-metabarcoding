include { resolve_runs } from '../../modules/local/fetch/resolve_runs.nf'
include { download_run } from '../../modules/local/fetch/download_run.nf'


workflow fetch_accessions {
    // [S97]/[S99]/[S101]: two-stage fetch. `resolve_runs` maps each
    // accession to its run accessions (one task per accession); the run
    // list is expanded into one (accession, run) pair per line, and
    // `download_run` fetches each run (one task per run). The per-run
    // granularity is what gives the [S101] failure contract: a failed
    // run fails only its own task, so -resume retries just the missing
    // runs while the succeeded ones stay cached.
    take:
    accessions   // channel of accession strings (one per --accession entry)

    main:
    resolve_runs(accessions)

    // [S99]→[S101]: expand <accession>_runs.tsv (one run accession per
    // line) into a (accession, run) channel — one emission per run — so
    // download_run fans out one task per run. Blank lines are dropped so
    // an accession that resolves to no runs simply yields no download.
    def per_run = resolve_runs.out.runs
        .flatMap { accession, runs_file ->
            runs_file.readLines()
                .collect { it.trim() }
                .findAll { it }
                .collect { run -> tuple(accession, run) }
        }

    download_run(per_run)

    emit:
    reads = download_run.out.reads
}
