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
split -d -a 4 -l ${splitSize} <(tail -n +2 metadata_esearch.csv) metadata_esearch_ --additional-suffix=.csv # tail gets rid of the header

echo "Fetching metadata from pysradb..."
for f in $(find -type f -name "metadata_esearch_*" | sort); do
    echo Processing $f
    temp=${f##*_}; num=${temp%.*}
    pysradb metadata $(cat $f | cut -f 1 -d ",") --detailed > metadata_pysradb_${num}.tsv
    sleep 1
done


# # Validate pysradb outputs and re-fetch any that are broken.
# # A file is considered VALID if:
# #   - it has more than 1 tab-separated column, AND
# #   - its header line contains "run_accession", AND
# #   - it does not start with "ValueError" / XML error text.
# is_valid_pysradb_file() {
#     local f="$1"
#     [ -s "$f" ] || return 1                       # empty file
#     local header
#     header=$(head -n 1 "$f")
#     # Reject obvious error dumps
#     case "$header" in
#         ValueError*|*"Unable to parse xml"*|*"<?xml"*) return 1 ;;
#     esac
#     # Header must contain run_accession and have >1 tab-separated field
#     echo "$header" | grep -q "run_accession" || return 1
#     [ "$(echo "$header" | awk -F'\t' '{print NF}')" -gt 1 ] || return 1
#     return 0
# }

# echo "Validating pysradb outputs..."
# invalid_files=()
# for f in $(find -type f -name "metadata_pysradb_*.tsv" | sort); do
#     if is_valid_pysradb_file "$f"; then
#         echo "  OK      : $f"
#     else
#         echo "  INVALID : $f"
#         invalid_files+=("$f")
#     fi
# done

echo "Validating pysradb outputs..."
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo "Script location: $SCRIPT_DIR"
# Python validator prints OK/BAD lines and exits 1 if any are bad.
# Capture the list of bad files for retry.
python3 "$SCRIPT_DIR/validate_pysradb.py" > validation.log 2>&1 || true
cat validation.log

mapfile -t invalid_files < <(grep '^BAD ' validation.log | awk '{print $2}')

if [ ${#invalid_files[@]} -eq 0 ]; then
    echo "All pysradb files valid."
else
    echo "${#invalid_files[@]} invalid file(s) need retry."
fi

# # Retry invalid files once (or more, if you like)
# max_retries=2
# for attempt in $(seq 1 $max_retries); do
#     [ ${#invalid_files[@]} -eq 0 ] && break
#     echo "Retry attempt $attempt for ${#invalid_files[@]} invalid file(s)..."
#     still_invalid=()
#     for f in "${invalid_files[@]}"; do
#         # Recover the chunk number from the filename: metadata_pysradb_007.tsv -> 007
#         temp=${f##*_}; num=${temp%.*}
#         chunk="metadata_esearch_${num}.csv"
#         if [ ! -s "$chunk" ]; then
#             echo "  No matching chunk for $f (expected $chunk); skipping retry."
#             still_invalid+=("$f")
#             continue
#         fi
#         echo "  Re-fetching $f from $chunk (attempt $attempt)"
#         pysradb metadata $(cat "$chunk" | cut -f 1 -d ",") --detailed > "$f"
#         sleep 2
#         if is_valid_pysradb_file "$f"; then
#             echo "    -> now valid: $f"
#         else
#             echo "    -> STILL invalid: $f"
#             still_invalid+=("$f")
#         fi
#     done
#     invalid_files=("${still_invalid[@]}")
# done

# Path to the Python validator (adjust if you already define SCRIPT_DIR elsewhere)
VALIDATOR="$SCRIPT_DIR/validate_pysradb.py"

# Retry invalid files once (or more, if you like)
max_retries=2
for attempt in $(seq 1 $max_retries); do
    [ ${#invalid_files[@]} -eq 0 ] && break
    echo "Retry attempt $attempt for ${#invalid_files[@]} invalid file(s)..."
    still_invalid=()
    for f in "${invalid_files[@]}"; do
        # Recover the chunk number from the filename: metadata_pysradb_007.tsv -> 007
        temp=${f##*_}; num=${temp%.*}
        chunk="metadata_esearch_${num}.csv"
        if [ ! -s "$chunk" ]; then
            echo "  No matching chunk for $f (expected $chunk); skipping retry."
            still_invalid+=("$f")
            continue
        fi
        echo "  Re-fetching $f from $chunk (attempt $attempt)"
        pysradb metadata $(cat "$chunk" | cut -f 1 -d ",") --detailed > "$f"
        sleep 2
        if python3 "$VALIDATOR" "$f" > /dev/null 2>&1; then
            echo "    -> now valid: $f"
        else
            echo "    -> STILL invalid: $f"
            # Optional: show the reason
            python3 "$VALIDATOR" "$f" 2>&1 | sed 's/^/       /'
            still_invalid+=("$f")
        fi
    done
    invalid_files=("${still_invalid[@]}")
done

# Report anything that never succeeded
if [ ${#invalid_files[@]} -gt 0 ]; then
    echo "WARNING: the following pysradb files could not be fetched after $max_retries retries:"
    printf '  %s\n' "${invalid_files[@]}"
    echo "They will be excluded from the merged output."
    # Rename them so the Python merge won't pick them up
    for f in "${invalid_files[@]}"; do
        mv "$f" "${f}.bad"
    done
fi

# rejoin pieces: each pysradb output may have different columns, so you can't just use the column header from the first file.
echo "Merging pysradb outputs..."

python3 << 'EOF'
import pandas as pd
import glob
import csv, io

# this version of repair_tsv does not work
# def repair_tsv(src, dst):
#     """Rewrite a TSV so each record is one physical line, quoting minimal."""
#     with open(src, "r", newline="", encoding="utf-8", errors="replace") as fin, \
#          open(dst, "w", newline="", encoding="utf-8") as fout:
#         reader = csv.reader(fin, delimiter="\t", quotechar='"')
#         writer = csv.writer(fout, delimiter="\t", quoting=csv.QUOTE_MINIMAL,
#                             lineterminator="\n")
#         try:
#             for row in reader:
#                 writer.writerow(row)
#         except csv.Error as e:
#             # csv.Error here means a truly unrecoverable quoting issue
#             print(f"csv.Error in {src}: {e}")
#             return False
#     return True

files = sorted(glob.glob("metadata_pysradb_*.tsv"))

dfs = []
all_columns = [] # to maintain expected order of columns, use list rather than set

# First pass: discover all columns
# This "for" block works but if a file is corrupted it might skip lines downstream of the corrupted line
for f in files:
    # df = pd.read_csv(
    #     f, sep="\t", dtype=str,
    #     engine="python",
    #     on_bad_lines="skip",
    # )
    # df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False,  on_bad_lines="warn")
    # df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False, quoting=csv.QUOTE_NONE) # the quoting option solves issues where there are unmatched quotes in the metadata, but might introduce issues if the string "\t" appears in the metadata
    try:
        # print(f)
        df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False, quoting=csv.QUOTE_NONE)
        # df = df.drop_duplicates(subset="run_accession", keep="first")
        dfs.append(df)
        all_columns += [col for col in df.columns if col not in all_columns]
    except:
        try:
            print(f"There was initially an error with parsing {f}; it may contain '\\t' interpreted as tab. Reattempting...")
            df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
            # df = df.drop_duplicates(subset="run_accession", keep="first")
            dfs.append(df)
            all_columns += [col for col in df.columns if col not in all_columns]
        except:
            print(f"Failed to parse {f}")

# for f in files:
#     try:
#         df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False, on_bad_lines="warn")
#     except pd.errors.ParserError as e:
#         print(f"ParserError in {f}: {e}")
#         repaired = f + ".repaired"
#         if not repair_tsv(f, repaired):
#             print(f"Could not repair {f}, skipping.")
#             continue
#         print(f"Using repaired file {repaired}")
#         df = pd.read_csv(repaired, sep="\t", dtype=str, low_memory=False, on_bad_lines="warn")

# all_columns = sorted(all_columns)

# Second pass: align each df to the full column set
aligned = [df.reindex(columns=all_columns) for df in dfs]

# Merge
merged = pd.concat(aligned, ignore_index=True)

# move run_accession column to front
col = merged.pop('run_accession')
merged.insert(0, 'run_accession', col)

# earlier even when I dropped duplicates I still ended up with duplicate IDs, so I can only assume they were appearing in multiple files
merged = merged.drop_duplicates(subset="run_accession", keep="first")

# Write final output
merged.to_csv("metadata_pysradb.tsv", sep="\t", index=False)
EOF

# delete temp files
rm -f $(find -name "metadata_esearch_*") # it's fine to delete this because in debugging you can easily re-split the esearch file
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