from pathlib import Path
import pandas as pd
import csv
import json
import argparse

# example run:
# python "/home/kcw2/metadata-magnet/scripts/sra_processing/screen_sra_coverage.py" --taxdir "/home/kcw2/data/testing/testdir1/taxonomy_analysis/" --genome_length 6.26e6 --taxid 286

# # the function inputs
# taxdir = "/home/kcw2/data/testing/testdir1/taxonomy_analysis/"
# genome_length = 6.26e6 # scientific notation is saved as a float I think
# coverage_threshold = 30 # when parsing this, set the type as int or float
# taxid = 286 # int or string?

def join_data(taxdir, genome_length, taxid, coverage_threshold=30, outdir='.', showAllMetdata=False):
    outdir_path = Path(outdir)
    outdir_path.mkdir(exist_ok=True, parents=True)

    esearch_file = f"{taxdir}/../metadata_esearch.csv"
    esearch = pd.read_csv(esearch_file, sep=",", on_bad_lines="warn")

    pysradb_file = f"{taxdir}/../metadata_pysradb.tsv"
    #pysradb = pd.read_csv(pysradb_file, sep="\t", low_memory=False) # get a warning about mixed data types without low_memory=False
    pysradb = pd.read_csv(pysradb_file, sep="\t", low_memory=False,  on_bad_lines="warn", quoting=csv.QUOTE_NONE) # some lines contain unmatched quotes; this ignores them

    coverage_by_run = {}
    pass_by_run = {}
    spots_by_run = {}

    pathlist = sorted(Path(taxdir).rglob('*.json'), key=lambda p: p.stem)

    for path in pathlist:
        run_id = path.stem # filename without extension- this is the run ID

        matches = pysradb.loc[pysradb.run_accession == run_id, "library_layout"]
        if matches.empty: # e.g. if NCBI glitched on this accession and pysradb failed to get metadata for it
            print(f"WARNING: {run_id} missing from pysradb metadata; skipping")
            
            # fill with bogus values; even if in theory the coverage is high enough, the pysradb metadata is not re-acquired later, so to avoid problems skip this
            coverage_by_run[run_id] = 0
            pass_by_run[run_id] = False
            spots_by_run[run_id] = 0
            
            continue
            # spot_factor = 2 # for the calculation, assume that the run is paired because that requirement is more stringent
            # but it doesn't really matter either way because later an inner join is performed, removing rows without pysradb metadata
        else:
            read_type = matches.iloc[0]
            spot_factor = 1 + (read_type == "PAIRED") # in paired end runs, a single spot has two reads

        avg_length = esearch[esearch.Run == run_id].reset_index(drop=True).avgLength[0]

        # parse the json
        with open(path, 'r') as file:
            raw = file.read().strip()

            # if empty file, skip
            if not raw:
                print(f"WARNING: empty JSON for {run_id}")
                coverage_by_run[run_id] = 0
                pass_by_run[run_id] = False
                spots_by_run[run_id] = 0
                continue

            # if invalid json, skip
            try:
                data = json.loads(raw)
            except json.JSONDecodeError:
                print(f"WARNING: invalid JSON for {run_id}")
                coverage_by_run[run_id] = 0
                pass_by_run[run_id] = False
                spots_by_run[run_id] = 0
                continue

            match = next((item for item in data[0]["tax_table"] if item["tax_id"] == taxid), None)
            num_spots = 0 if match is None else match["total_count"]

        coverage = (spot_factor * num_spots * avg_length) / genome_length
        coverage_by_run[run_id] = coverage
        pass_by_run[run_id] = coverage >= coverage_threshold
        spots_by_run[run_id] = num_spots

    esearch[f"coverage_taxid_{taxid}"] = esearch.Run.astype(str).map(coverage_by_run)
    # esearch["coverage_pass"] = esearch.Run.astype(str).map(pass_by_run)
    # esearch["spots"] = esearch.Run.astype(str).map(spots_by_run)

    df = pd.DataFrame({
        'run_id': esearch.Run.astype(str),
        'coverage': esearch.Run.astype(str).map(coverage_by_run),
        'coverage_pass': esearch.Run.astype(str).map(pass_by_run),
        f"spots_for_taxid": esearch.Run.astype(str).map(spots_by_run)
    })

    passed = df[df.coverage_pass == True].run_id

    df.to_csv(f"{outdir}/coverage_taxid_{taxid}.tsv", sep='\t', index=False) # may just name this coverage.tsv, though it's possible to use wildcards to find this tsv
    
    print("Number of runs with sufficient coverage:", len(passed))

    # print(taxdir, genome_length, taxid, coverage_threshold)
    passed.to_csv(f"{outdir}/coverage_pass.txt", sep='\t', index=False, header=False)
    
    # join the passed IDs to the metadata
    if showAllMetdata:
        passed_df = df
    else:
        passed_df = df[df.coverage_pass == True]

    df_full = pd.merge(passed_df, esearch, left_on = "run_id", right_on = "Run").drop(columns="Run")
    df_full = pd.merge(df_full, pysradb, left_on = "run_id", right_on = "run_accession").drop(columns="run_id") # default how='inner'
    
    # which runs are most space-efficient for coverage?
    df_full["coverage_megabyte_ratio"] = df_full[f"coverage"] / df_full["size_MB"]
    
    # reorganize df: move run accession to front, and rename ambiguously-named columns
    for c in ["total_spots", "spots_for_taxid", "coverage_megabyte_ratio", "coverage", "run_accession"]:
      col = df_full.pop(c)
      df_full.insert(0, col.name, col)
    
    df_full.rename(columns={"avgLength":"avg_read_length", "size_MB": "size_megabytes"}, inplace=True)

    if showAllMetdata:
        df_full.to_csv(f"{outdir}/coverage_taxid_{taxid}_allWithMetadata.tsv", sep='\t', index=False)
    else:
        df_full.drop('coverage_pass', axis=1).to_csv(f"{outdir}/coverage_taxid_{taxid}_passedWithMetadata.tsv", sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to filter SRA runs by coverage for a taxid of choice.")

    # positional arguments (required)
    parser.add_argument("-t", "--taxdir", type=str, help='Directory containing SRA taxonomy analyses.')
    parser.add_argument("-L", "--genome_length", type=float, help="Length of reference genome in bp.")
    parser.add_argument("-i", "--taxid", type=int, help="Taxonomic id to reference when calculating coverage.") # must be int or else it won't match in the JSON
    parser.add_argument("-c", "--coverage_threshold", type=float, default=30, help="Coverage required for the taxonomic ID.")
    parser.add_argument("-o", "--outdir", type=str, default='.', help='Output directory (default current directory).')
    parser.add_argument("-a", "--showAllMetadata", type=str, default="False", help='If true, joins metadata to full dataset rather than just the ones that passed coverage check')

    args = parser.parse_args()
    showAll_bool = args.showAllMetadata == "True"
    print(args.taxdir, args.genome_length, args.taxid, args.coverage_threshold, args.outdir, showAll_bool)
    join_data(args.taxdir, args.genome_length, args.taxid, args.coverage_threshold, args.outdir, showAll_bool)