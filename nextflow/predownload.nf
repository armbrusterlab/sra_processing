#!/usr/bin/env nextflow
include { getMetadata } from './modules/getMetadata.nf'
include { checkTaxidCoverage } from './modules/checkCoverage.nf'
include { predictEnvSources } from './modules/predictSources.nf'

params {
    // always use params file to pass arguments; otherwise ints, floats, and Booleans may be interpreted as strings
    sra_query: String
    genome_length: Float
    taxid: Integer
    coverage_threshold: Integer
    model_dir: Path
}

workflow {
    main:
    getMetadata(params.sra_query)
    def taxdir = getMetadata.out.taxonomy_analysis

    checkTaxidCoverage(taxdir, params.genome_length, params.taxid, params.coverage_threshold)

    predictEnvSources(getMetadata.out.metadata_pysradb, params.model_dir)

    publish:
    metadata_esearch = getMetadata.out.metadata_esearch
    metadata_pysradb = predictEnvSources.out.metadata_pysradb // it's overwritten to include model predictions
    taxonomy_analysis = taxdir

    taxid = checkTaxidCoverage.out.taxid
    taxid_passed = checkTaxidCoverage.out.taxid_passed
    taxid_passed_list = checkTaxidCoverage.out.taxid_passed_list
}


output {
    metadata_esearch {
        path { "metadata" }
        mode 'copy'
    }
    metadata_pysradb {
        path { "metadata" }
        mode 'copy'
    }
    taxonomy_analysis {
        path { "metadata" }
        mode 'copy'
    }

    taxid {
        path { "metadata/coverage_check" }
        mode 'copy'
    }
    taxid_passed {
        path { "metadata/coverage_check" }
        mode 'copy'
    }
    taxid_passed_list {
        path { "metadata/coverage_check" }
        mode 'copy'
    }
}