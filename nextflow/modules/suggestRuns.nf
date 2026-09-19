/*
 * Suggest runs to download, assuming the user's goal is to statistically compare mutants between categories.
 */
process suggestRuns {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path taxid_passed
    val terms_colname

    output:
    path "suggested_runs/", emit: suggested

    script:
    """
    projDir="${workflow.projectDir}"
    
    python "\$projDir/../scripts/suggest_sra_runs.py" -f ${taxid_passed} -o "suggested_runs" -t "${terms_colname}"
    """
}