#!/usr/bin/env nextflow
include { downloadRuns } from './modules/downloadRuns.nf'
include { qualityControl } from './modules/qualityControl.nf'
include { krakenClassify; makePatternFiles; filterByTaxid } from './modules/taxidFiltering.nf'
include { makeGB; runBreseq; summarizeBreseq; summarizeSources } from './modules/findVariants.nf'
include { variantStats } from './modules/statisticalTests.nf'

params {
    // always use params file to pass arguments; otherwise ints, floats, and Booleans may be interpreted as strings
    run_list: Path
    predownload_outputs: Path
    read_length_threshold: Integer
    
    fastp_additional_short: String
    fastp_additional_long: String
    genome_length: Float
    coverage_threshold: Integer

    k2db: Path
    taxid: Integer

    reference_gb: Path
    target_genes: String
    target_type: String
    buffer_upstream: Integer
    buffer_downstream: Integer

    breseq_additional: String
    filter_intergenic: String
    filter_synonymous: String
    category_colname: String
    subcategory_colname: String
    terms_colname: String

    p_adjust_method: String
    report_all: String
}

workflow {
    main:
    downloadRuns(params.run_list, params.predownload_outputs, params.read_length_threshold)
    def fqdump = downloadRuns.out.fqdump
    def runids_dir = downloadRuns.out.runids_dir

    qualityControl(fqdump, params.fastp_additional_short, params.fastp_additional_long, runids_dir, params.predownload_outputs, params.genome_length, params.coverage_threshold)
    def fastp = qualityControl.out.fastp
    def runids_postqc_dir = qualityControl.out.runids_postqc_dir

    krakenClassify(fastp, runids_postqc_dir, params.k2db)
    def kraken2 = krakenClassify.out.kraken2

    makePatternFiles(params.taxid)
    def patternsdir = makePatternFiles.out.patternsdir

    filterByTaxid(kraken2, patternsdir)
    def grepq = filterByTaxid.out.grepq

    makeGB(params.reference_gb, params.target_genes, params.target_type, params.buffer_upstream, params.buffer_downstream)
    def gb_for_breseq = makeGB.out.gb_for_breseq

    // runBreseq(grepq, runids_postqc_dir, gb_for_breseq, params.breseq_additional, params.predownload_outputs, params.filter_intergenic, params.filter_synonymous)
    // def breseq_tables = runBreseq.out.breseq_tables
    // def breseq_htmls = runBreseq.out.breseq_htmls

    runBreseq(grepq, runids_postqc_dir, gb_for_breseq, params.breseq_additional)
    def breseq_htmls = runBreseq.out.breseq_htmls

    summarizeBreseq(breseq_htmls, params.predownload_outputs, params.filter_intergenic, params.filter_synonymous)
    def breseq_tables = summarizeBreseq.out.breseq_tables

    // summarizeSources(breseq_tables, params.category_colname, params.subcategory_colname)

    variantStats(breseq_tables, params.predownload_outputs, params.terms_colname, params.p_adjust_method, params.report_all)

    publish:
    kraken_reports = krakenClassify.out.kraken2_reports
    breseq_t = breseq_tables
    breseq_h = breseq_htmls
    // source_summary = summarizeSources.out.mutation_frequencies
    stats_dir = variantStats.out.stats
}


output {
    kraken_reports {
        path { "./" }
        mode "copy"
    }
    breseq_t {
        path { "./" }
        mode "copy"
    }
    breseq_h {
        path { "./" }
        mode "copy"
    }
    // source_summary {
    //     path {"./"}
    //     mode "copy"
    // }
    stats_dir {
        path {"./"}
        mode "copy"
    }
}