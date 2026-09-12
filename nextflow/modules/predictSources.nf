/*
 * Predicts isolation source category and subcategory from pysradb metadata.
 */
process predictEnvSources {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path metadata, stageAs: 'metadata.tsv'
    path model_dir

    output:
    path "metadata_pysradb.tsv", emit: metadata_pysradb

    script:
    """
    projDir="${workflow.projectDir}"
    python "\$projDir/../scripts/predict_environmental_source.py" -m "metadata.tsv" -d "${model_dir}" -o "metadata_pysradb.tsv"
    """
}