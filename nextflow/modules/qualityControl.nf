/*
 * Use fastp/fastplong for quality control based on read type
 */
process qualityControl {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path fqdump // from downloadRuns process
    val fastp_additional_short
    val fastp_additional_long

    path runids_dir
    path predownload_outputs // from part 1: the predownload pipeline
    val genome_length
    val coverage_threshold

    output:
    path "fastp/", emit: fastp
    path "run_ids_postqc/", emit: runids_postqc_dir

    script:
    """
    projDir="${workflow.projectDir}"

    # short reads
    if [[ -z "${fqdump}/paired/" ]]; then
        echo "No paired-end reads."
        mkdir -p "fastp/paired/reports"
    else
        python "\$projDir/../scripts/fastp_parallel.py" -i "${fqdump}/paired/" -o fastp/paired/ -r fastp/paired/reports -1 "_1" -2 "_2" -a "${fastp_additional_short}"
    fi

    if [[ -z "${fqdump}/single/short/" ]]; then
        echo "No short single-end reads."
        mkdir -p "fastp/single/short/reports"
    else
        python "\$projDir/../scripts/fastp_parallel.py" -i "${fqdump}/single/short/" -o fastp/single/short/ -r fastp/single/short/reports -a "${fastp_additional_short}"
    fi

    # long reads
    if [[ -z "${fqdump}/single/long/" ]]; then
        echo "No long single-end reads."
        mkdir -p "fastp/single/long/reports"
    else
        python "\$projDir/../scripts/fastplong_parallel.py" -i "${fqdump}/single/long/" -o fastp/single/long/ -r fastp/single/long/reports -a "${fastp_additional_long}"
    fi

    mkdir "run_ids_postqc/"
    python "\$projDir/../scripts/refilter_after_QC.py" -q "fastp/" -i ${runids_dir} -m \$(find "${predownload_outputs}/coverage_check/" -name "*_passedWithMetadata.tsv") -L ${genome_length} -c ${coverage_threshold}

    # filter the singletons file by intersecting with the filtered PE ids list
    comm -12 <(cat "${runids_dir}/ids_singletons.txt" | sort) <(cat "run_ids_postqc/ids_paired_postQC.txt" | sort) > "run_ids_postqc/ids_singletons_postQC.txt"
    """
}