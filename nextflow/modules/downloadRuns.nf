/*
 * Download specified runs and sort them by read type. 
 */
process downloadRuns {
    conda "${workflow.projectDir}/envs/envs.yml"

    input:
    path run_list // user-curated newline-separated list of SRA accessions; no header
    path predownload_outputs // from part 1: the predownload pipeline
    val read_length_threshold // runs longer than this are considered long reads

    output:
    path "fasterq-dump/", emit: fqdump
    path "run_ids/", emit: runids_dir

    script:
    """
    pysradb="${predownload_outputs}/metadata_pysradb.tsv" # info on single vs paired reads
    esearch="${predownload_outputs}/metadata_esearch.csv" # avgLengths in col 2

    # download runs
    while IFS= read -r line; do
        echo "Prefetch: \$line"
        prefetch "\$line" -O prefetch
    done < "${run_list}"

    while IFS= read -r line; do
        echo "Downloading: \$line"
        fasterq-dump "prefetch/\${line}" -O "fasterq-dump" --split-3 --skip-technical
    done < "${run_list}"

    # the prefetch data can be deleted now that fasterq-dump has retrieved the full-size files
    rm -rf "prefetch/"

    # sort runs
    mkdir "run_ids/"

    # first separate by SE vs PE
    # the pysradb output columns aren't always in the same order so we need to find the column by name (library_layout)
    # first get the full set of IDs, no filtering whatsoever 
    > "run_ids/ids_single.txt"
    > "run_ids/ids_paired.txt"
    awk -F'\t' -v col1="library_layout" -v col2="run_accession" '
    NR==1 {
        # find column indices
        for (i=1; i<=NF; i++) {
            if (\$i == col1) c1 = i
            if (\$i == col2) c2 = i
        }
        next
    }
    \$c1 == "SINGLE" { print \$c2 >> "run_ids/ids_single.txt" }
    \$c1 == "PAIRED" { print \$c2 >> "run_ids/ids_paired.txt" }
    ' \$pysradb

    # further subdivide the single reads into long and short
    comm -12 <(cat "run_ids/ids_single.txt" | sort) <(awk -F ',' -v t="${read_length_threshold}" '\$2 <= t {print \$1}' "\$esearch" | sort) > "run_ids/ids_single_short.txt"
    comm -12 <(cat "run_ids/ids_single.txt" | sort) <(awk -F ',' -v t="${read_length_threshold}" '\$2 > t {print \$1}' "\$esearch" | sort) > "run_ids/ids_single_long.txt"

    # intersect with the sorted list of manually curated runs, since the lists above may contain runs that don't actually exist after filtering.
    cat "${run_list}" | sort > "run_ids/filtered_sorted.txt"
    comm -12 "run_ids/filtered_sorted.txt" <(cat "run_ids/ids_paired.txt" | sort) > "run_ids/ids_paired_filtered.txt"
    comm -12 "run_ids/filtered_sorted.txt" <(cat "run_ids/ids_single_long.txt" | sort) > "run_ids/ids_single_long_filtered.txt"
    comm -12 "run_ids/filtered_sorted.txt" <(cat "run_ids/ids_single_short.txt" | sort) > "run_ids/ids_single_short_filtered.txt"

    # move files accordingly
    IdLists=("run_ids/ids_paired_filtered.txt" "run_ids/ids_single_long_filtered.txt" "run_ids/ids_single_short_filtered.txt")
    Dirs=("fasterq-dump/paired" "fasterq-dump/single/long" "fasterq-dump/single/short")
    length=\${#Dirs[@]} # should be 3

    for ((i = 0; i < \$length; i++)); do
        echo "Index \$i: \${IdLists[i]} \${Dirs[i]}"
        mkdir -p "\${Dirs[i]}"

        while IFS= read -r id; do
            echo "Processing \$id"
            mv \$(find "fasterq-dump/" -type f -name "\$id*") "\${Dirs[i]}"
        done < "\${IdLists[i]}"
    done

    # handle singletons
    echo "Checking for singletons..."
    singletons=\$(find "fasterq-dump/paired/" -type f ! -name *_1.fastq ! -name *_2.fastq)
    find "fasterq-dump/paired/" -type f ! -name *_1.fastq ! -name *_2.fastq |  xargs -n 1 basename | sed 's/\\.[^.]*\$//' >> "run_ids/ids_singletons.txt" 
    # normally it would be sed 's/\.[^.]*$//', but in the Nextflow context, must escape backslash and dollar sign

    if [[ -z \$singletons ]]; then 
        echo "No singletons to move."
    else 
        echo "Singletons found: \$singletons"
        mv \$singletons "fasterq-dump/single/short/"        
    fi
    """
}