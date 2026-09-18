/*
 * Returns runs with sufficient coverage for the given taxid.
 */
process checkTaxidCoverage {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path taxdir
    path predictions // table of predicted environments; no longer part of pysradb table
    val genome_length
    val taxid
    val coverage_threshold

    output:
    path "coverage_taxid_*.tsv", emit: taxid
    path "coverage_taxid_*_passedWithMetadata.tsv", emit: taxid_passed
    path "coverage_pass.txt", emit: taxid_passed_list

    script:
    """
    projDir="${workflow.projectDir}"
    python "\$projDir/../scripts/screen_sra_coverage.py" --taxdir "${taxdir}" --predictions "${predictions}" --genome_length ${genome_length} --taxid ${taxid} --coverage_threshold ${coverage_threshold}
    """
}