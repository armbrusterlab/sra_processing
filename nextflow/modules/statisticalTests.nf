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

    output:
    path "stats/", emit: stats

    script:
    """
    projDir="${workflow.projectDir}"

    mutants="${breseq_tables}/mutations.tsv" # at this time, only the mutations table is used for stats
    predictions="${predownload_outputs}/environment_predictions.tsv"
    
    Rscript "\$projDir/../scripts/mutant_environment_tests.R" "\$mutants" "\$predictions" "stats/" "${terms_colname}" "${p_adjust_method}" "${report_all}"
    """
}