/*
 * Predicts isolation source category and subcategory from pysradb metadata.
 */
process predictEnvSources {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path metadata
    path model_dir

    output:
    path "environment_predictions.tsv", emit: predictions

    script:
    """
    projDir="${workflow.projectDir}"
    python "\$projDir/../scripts/predict_environmental_source.py" -m "${metadata}" -d "${model_dir}" -o "environment_predictions.tsv"
    """
}