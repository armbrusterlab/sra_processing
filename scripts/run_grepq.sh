#!/bin/bash
# Intended to be run by Nextflow, so the output (grepqdir) is created in the current directory
set -euo pipefail

export patternsdir=$1
kraken=$2
batchname=$3 # paired/, single/short/, or single/long/

krakendir="${kraken}/${batchname}"
export grepqdir="grepq/${batchname}"
mkdir -p $grepqdir

numCores=$(( $(nproc) / 4 ))

# -L to make find follow symlinks, which is what Nextflow will pass to it. Needed for all find calls.
find -L "$krakendir" -type f -name "*.fq" | parallel --env grepqdir --env patternsdir -j "$numCores" '
    f={} 
    grepq_outputs="${grepqdir}/temp_$(basename "$f" .fq).txt"
    echo "Processing $f; output file list: ${grepq_outputs}"

    for patterns in $(find -L "$patternsdir" -type f | sort); do
        num=$(basename ${patterns%.*})
        base=$(basename "$f" .fq);
        id=${base%_*};
        mkdir -p "${grepqdir}/${id}/"
        fragment="${grepqdir}/${id}/$(basename "$f" .fq)_${num}.fq.gz"
        echo "Processing $f with IDs part $num"
        grepq -R --write-gzip --best "$patterns" "$f" > $fragment
        echo "$fragment" >> "$grepq_outputs"
    done

    cat $(cat "$grepq_outputs") > "${grepqdir}/${id}/$(basename "$f" .fq).fq.gz"
    rm -f $(cat "$grepq_outputs")
    rm -f "$grepq_outputs"
'