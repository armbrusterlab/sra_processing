/*
 * Use kraken2 to classify reads by taxid, then filter reads to those matching a specified taxid or its children. 
 */
process krakenClassify {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    // from qualityControl process
    path fastp 
    path runids_postqc_dir

    path k2db // user must prepare a kraken2 database in advance

    output:
    path "kraken2/", emit: kraken2
    path "kraken2_reports/", emit: kraken2_reports

    script:
    """
    numCores=\$(( \$(nproc) / 4 ))

    echo "Processing paired reads."
    start=\$(date +"%Y-%m-%d %H:%M:%S")

    mkdir -p "kraken2/paired/"
    while IFS= read -r id; do
        files=\$(find ${fastp}/paired/ -mindepth 1 -maxdepth 1 -type f -name "\$id*.clean.fastq.gz" | sort)
        echo "Running kraken2 on \$id"
        krakenout="kraken2/paired/\${id}"
        mkdir -p \$krakenout
        reportout="kraken2_reports/paired"
        mkdir -p \$reportout
        k2 classify --db ${k2db} --use-daemon --threads \$numCores --paired --classified-out \$krakenout/\${id}#.fq --output "-" --report \${reportout}/\${id}.txt \$files # the "#" is replaced by _1 and _2
    done < "${runids_postqc_dir}/ids_paired_postQC.txt"

    end=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Processing single reads, short"
    start2=\$(date +"%Y-%m-%d %H:%M:%S")

    mkdir -p "kraken2/single/short/"
    while IFS= read -r id; do
        files=\$(find ${fastp}/single/short/ -mindepth 1 -maxdepth 1 -type f -name "\$id*.clean.fastq.gz" | sort)
        echo "Running kraken2 on \$id"
        krakenout="kraken2/single/short/\${id}"
        mkdir -p \$krakenout
        reportout="kraken2_reports/single/short"
        mkdir -p \$reportout
        k2 classify --db ${k2db} --use-daemon --threads \$numCores --classified-out \$krakenout/\${id}.fq --output "-" --report \${reportout}/\${id}.txt \$files
    done < "${runids_postqc_dir}/ids_single_short_postQC.txt"

    # also process singletons
    while IFS= read -r id; do
        files=\$(find ${fastp}/single/short/ -mindepth 1 -maxdepth 1 -type f -name "\$id*.clean.fastq.gz" | sort)
        echo "Running kraken2 on singleton \$id"
        krakenout="kraken2/single/short/\${id}"
        mkdir -p \$krakenout
        reportout="kraken2_reports/single/short"
        mkdir -p \$reportout
        k2 classify --db ${k2db} --use-daemon --threads \$numCores --classified-out \$krakenout/\${id}.fq --output "-" --report \${reportout}/\${id}.txt \$files
    done < "${runids_postqc_dir}/ids_singletons_postQC.txt"

    end2=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Processing single reads, long"
    start3=\$(date +"%Y-%m-%d %H:%M:%S")

    mkdir -p "kraken2/single/long/"
    while IFS= read -r id; do
        files=\$(find ${fastp}/single/long/ -mindepth 1 -maxdepth 1 -type f -name "\$id*.clean.fastq.gz" | sort)
        echo "Running kraken2 on \$id"
        krakenout="kraken2/single/long/\${id}"
        mkdir -p \$krakenout
        reportout="kraken2_reports/single/long"
        mkdir -p \$reportout
        k2 classify --db ${k2db} --use-daemon --threads \$numCores --classified-out \$krakenout/\${id}.fq --output "-" --report \${reportout}/\${id}.txt \$files
    done < "${runids_postqc_dir}/ids_single_long_postQC.txt"

    end3=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Start time: \$start"
    echo "End time: \$end"
    echo "Start2 time: \$start2"
    echo "End2 time: \$end2"
    echo "Start3 time: \$start3"
    echo "End3 time: \$end3"

    k2 clean --stop-daemon
    """
}
 
process makePatternFiles {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    val taxid

    output:
    path "taxid_patterns/", emit: patternsdir

    script:
    """
    projDir="${workflow.projectDir}"

    echo "Downloading children of taxid ${taxid}..."
    datasets download taxonomy taxon "${taxid}" --children --filename for_taxids.zip
    unzip -o -d for_taxids "for_taxids.zip"
    tail -n +2 "for_taxids/ncbi_dataset/data/taxonomy_summary.tsv" | cut -f 2 | sort -n > taxids.txt # the fixed part of the filename is ncbi_dataset/data/taxonomy_summary.tsv

    # to split the taxids
    echo "Splitting taxids..."
    n=1000
    mkdir -p taxids
    split -d -a 2 -l \$n "taxids.txt" "taxids/" --additional-suffix=.txt # to make it easier to extract the part numbers, I just made that the name of the file

    # make the JSONs (this requires yq to be installed)
    echo "Creating JSONs..."
    for f in taxids/*; do 
        bash "\$projDir/../scripts/taxids_to_patternsJson.sh" \$f 
    done # makes the jsons in the same dir as the taxid files

    # move the jsons to their own dir for convenience
    mkdir -p "taxid_patterns/"
    mv \$(find "taxids/" -name "*.json") "taxid_patterns/"

    # at this point you can delete the entire taxids dir.
    rm -rf "taxids/"

    echo "Done!"
    """
}

process filterByTaxid {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path kraken2
    path patternsdir

    output:
    path "grepq/", emit: grepq

    script:
    """
    projDir="${workflow.projectDir}"

    echo "Filtering paired-end reads..."
    start=\$(date +"%Y-%m-%d %H:%M:%S")
    bash "\$projDir/../scripts/run_grepq.sh" ${patternsdir} ${kraken2} "paired/"
    end=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Filtering single-end short reads..." # don't need to specify singletons; it processes all fastqs in this dir
    start2=\$(date +"%Y-%m-%d %H:%M:%S")
    bash "\$projDir/../scripts/run_grepq.sh" ${patternsdir} ${kraken2} "single/short/"
    end2=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Filtering single-end long reads..."
    start3=\$(date +"%Y-%m-%d %H:%M:%S")
    bash "\$projDir/../scripts/run_grepq.sh" ${patternsdir} ${kraken2} "single/long/"
    end3=\$(date +"%Y-%m-%d %H:%M:%S")

    echo "Start time: \$start"
    echo "End time: \$end"
    echo "Start2 time: \$start2"
    echo "End2 time: \$end2"
    echo "Start3 time: \$start3"
    echo "End3 time: \$end3"
    """
}