/*
 * Retrieves metadata for the SRA query specified.
 */
process getMetadata {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    val sra_query

    output:
    path "metadata_esearch.csv", emit: metadata_esearch
    path "metadata_pysradb.tsv", emit: metadata_pysradb
    path "taxonomy_analysis/", emit: taxonomy_analysis

    script:
    """
    projDir="${workflow.projectDir}"
    bash "\$projDir/../scripts/get_sra_metadata.sh" "${sra_query}"
    """
}