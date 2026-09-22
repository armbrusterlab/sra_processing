/*
 * Creates a GB file containing only target genes to use in the downstream breseq process
 */
process makeGB {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path reference_gb // user may download from NCBI or supply their own
    val target_genes // space-separated string of target genes
    val target_type // e.g. "locus_tag" if the reference genes are provided as locus tags from the reference genome
    val buffer_upstream // in bases
    val buffer_downstream // in bases

    output:
    path "for_breseq.gb", emit: gb_for_breseq

    script:
    """
    projDir="${workflow.projectDir}"

    python "\$projDir/../scripts/extract_gene_from_gb_v2.py" -i ${reference_gb} -o "for_breseq.gb" -g ${target_genes} -t ${target_type} -a "unknown" -u ${buffer_upstream} -d ${buffer_downstream}
    """
}

/*
 * Runs breseq to find variants of target genes among downloaded SRA data
 */
process runBreseq {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path grepq
    path runids_postqc_dir
    path gb_for_breseq
    val breseq_additional

    output:
    path "breseq_export/", emit: breseq_htmls

    script:
    """
    export projDir="${workflow.projectDir}"

    # run breseq
    numCores=${task.cpus}
    export ref="${gb_for_breseq}"
    export out="breseq/" 
    export additional="${breseq_additional}"

    start=\$(date +"%Y-%m-%d %H:%M:%S") # human-readable time

    # any singletons are processed with the paired IDs
    cat "${runids_postqc_dir}/ids_paired_postQC.txt" "${runids_postqc_dir}/ids_single_long_postQC.txt" "${runids_postqc_dir}/ids_single_short_postQC.txt" | parallel -j "\$numCores" '
        id={}
        echo "Processing \$id"
        files=\$(find -L ${grepq} -type f -name "\$id*")
        echo \$files
        required="-r \$ref -o \${out}/\${id} \$files"
        foo=\$(echo \$files | awk "{print \\\$1}")
        type=\$(basename \$(dirname -- \$(dirname -- \$foo)))
        echo "Read length: \$type"
        if [[ \$type == "long" ]]; then 
            required="-x \$required"
        fi
        ADDITIONAL=\$additional REQUIRED=\$required bash "\$projDir/../scripts/breseq_command.sh"
        '

    end=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "Start time: \$start"
    echo "End time: \$end"

    # move breseq output htmls to a single dir
    exportdir="breseq_export/"
    mkdir -p \$exportdir
    for d in \$(find breseq/ -mindepth 2 -type d -name "output"); do 
        temp=\${d%/*}
        run=\${temp##*/}
        mv \$d "\$exportdir/\${run}"
    done
    """
}
/*
 * The breseq summarizer has been separated into its own process so that modifications won't invalidate the cache for actually running breseq.
 */
process summarizeBreseq {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path breseq_tables
    path predownload_outputs
    val filter_intergenic
    val filter_synonymous

    output:
    path "breseq_summary_tables/", emit: breseq_tables

    script:
    """
    projDir="${workflow.projectDir}"

    find -L ${breseq_tables} -mindepth 1 -maxdepth 1 -type d | sort > "run_outputs.txt" # -L option so that find follows symlinks

    additional_flags=""
    if [[ "${filter_intergenic}" == "True" ]]; then
        additional_flags+="-i "
    fi
    if [[ "${filter_synonymous}" == "True" ]]; then
        additional_flags+="-s"
    fi
    echo "Additional flags: \$additional_flags"

    python "\$projDir/../scripts/summarize_breseq.py" "run_outputs.txt" "breseq_summary_tables" \$additional_flags -n -1

    # Removed the join_breseq_metadata.py call since the output file has a lot of redundant data, but the prediction data (formerly from pysradb table) will be used in downstream analysis
    """
}
/*
 * For each mutant, counts the associated environments. This process is no longer being used.
 */
process summarizeSources {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path breseq_tables
    val category_colname
    val subcategory_colname

    output:
    path "mutation_frequencies.tsv", emit: mutation_frequencies

    script:
    """
    projDir="${workflow.projectDir}"

    mutant_summary="${breseq_tables}/breseq_summary_withMetadata.tsv"
    python "\$projDir/../scripts/summarize_mutants_envsource.py" -f \$mutant_summary -c "${category_colname}" -s "${subcategory_colname}" # default outname is "mutation_frequencies.tsv"
    """
}