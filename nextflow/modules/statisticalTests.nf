/*
 * Run statistical tests on breseq output.
 */
process variantStats {
    conda "${workflow.projectDir}/envs/stats_env.yml"

    input:
    path breseq_tables
    path predownload_outputs
    val terms_colname
    val p_adjust_method
    val report_all
    val adjust_separately

    output:
    path "stats/", emit: stats

    script:
    """
    projDir="${workflow.projectDir}"

    mutants="${breseq_tables}/mutations.tsv" # at this time, stat tests are only run on the mutations table
    predictions="${predownload_outputs}/environment_predictions.tsv"
    
    # each output table is also saved as a UTF-16LE version which displays without garbled text in Excel
    Rscript "\$projDir/../scripts/mutant_environment_tests.R" "\$mutants" "\$predictions" "stats/" "${terms_colname}" "${p_adjust_method}" "${report_all}" "${adjust_separately}"
    """
}