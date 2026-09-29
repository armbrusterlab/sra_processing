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

    mutants="${breseq_tables}/mutations.tsv" # at this time, stat tests are only run on the mutations table
    predictions="${predownload_outputs}/environment_predictions.tsv"
    
    # TSVs are saved with UTF-16 encoding
    Rscript "\$projDir/../scripts/mutant_environment_tests.R" "\$mutants" "\$predictions" "stats/" "${terms_colname}" "${p_adjust_method}" "${report_all}"
    """
}