#!/bin/bash
# Intended to be run by Nextflow, so files and dirs are created in the current dir
# $1: NCBI query

echo "Fetching metadata from NCBI esearch..."
esearch -db sra -query "$1" | efetch -format runinfo | cut -d ',' -f 1,7,8 | awk 'NR==1{print;next} {print| "sort -k1,1n -t\",\""}' > metadata_esearch.csv

# if NCBI server access is interrupted, it'll write error messages to the file, so need to remove these lines to avoid downstream problems.
mv metadata_esearch.csv raw_metadata_esearch.csv
head -n 1 raw_metadata_esearch.csv > metadata_esearch.csv
grep "^[DES]R[RXSP]*" raw_metadata_esearch.csv >> metadata_esearch.csv

# must break file into smaller pieces to avoid error from pysradb
splitSize=100
split -d -a 2 -l ${splitSize} <(tail -n +2 metadata_esearch.csv) metadata_esearch_ --additional-suffix=.csv # tail gets rid of the header

echo "Fetching metadata from pysradb..."
for f in $(find -type f -name "metadata_esearch_*" | sort); do
    echo Processing $f
    temp=${f##*_}; num=${temp%.*}
    pysradb metadata $(cat $f | cut -f 1 -d ",") --detailed > metadata_pysradb_${num}.tsv
    sleep 1
done

# rejoin pieces: each pysradb output may have different columns, so you can't just use the column header from the first file.
echo "Merging pysradb outputs..."

python3 << 'EOF'
import pandas as pd
import glob

files = sorted(glob.glob("metadata_pysradb_*.tsv"))

dfs = []
all_columns = [] # to maintain expected order of columns, use list rather than set

# First pass: discover all columns
for f in files:
    df = pd.read_csv(f, sep="\t", dtype=str)
    dfs.append(df)
    all_columns += [col for col in df.columns if col not in all_columns]

# all_columns = sorted(all_columns)

# Second pass: align each df to the full column set
aligned = [df.reindex(columns=all_columns) for df in dfs]

# Merge
merged = pd.concat(aligned, ignore_index=True)

# move run_accession column to front
col = merged.pop('run_accession')
merged.insert(0, 'run_accession', col)

# Write final output
merged.to_csv("metadata_pysradb.tsv", sep="\t", index=False)
EOF

# delete temp files
rm -f $(find -name "metadata_esearch_*")
rm -f $(find -name "metadata_pysradb_*")

echo "Fetching taxonomy analysis metadata..."
mkdir -p "taxonomy_analysis/"
for run in $(tail -n +2 "metadata_esearch.csv" | cut -f 1 -d ","); do
    echo "Processing $run"
    wget -O "taxonomy_analysis/${run}.json" "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy?&acc=${run}&cluster_name=public"
    sleep 0.75
done

sleep 1

# reattempt downloads for runs ($9) of size ($4) == 0
for run in $(ls -lh "taxonomy_analysis/" | awk -F' ' '$5 == 0 {print $9}'); do
    run=${run%.*} # get rid of .json extension
    echo "Reattempting download for $run"
    wget -O "taxonomy_analysis/${run}.json" "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy?&acc=${run}&cluster_name=public"
    sleep 1
done

# if the second attempt fails, record the name in a file
for run in $(ls -lh "taxonomy_analysis/" | awk -F' ' '$5 == 0 {print $9}'); do
    echo ${run%.*} >> "taxonomy_analysis_download_failed.txt"
    # rm -f $(find taxonomy_analysis/ -name $run)
done