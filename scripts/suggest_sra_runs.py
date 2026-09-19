#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import argparse

# This suggestion assumes that you want to perform statistical comparison between environemts.
# It doesn't need to be followed, but it's a good start unless the user has clear ideas about what data they want to download 
def suggest_runs(f, outdir, terms_colname="terms_logistic_regression"):
    # f: the passedWithMetadata file from the coverage check
    # params to add when adding feature to filter out terms: category_colname = "category_logistic_regression", subcategory_colname = "subcategory_logistic_regression"
    outdir_path = Path(outdir)
    outdir_path.mkdir(exist_ok=True, parents=True)

    df = pd.read_csv(f, sep="\t")

    # first: assume the user wants to download only one run per study, as runs from the same study aren't independent
    # from each study, pick the run with the best ratio of coverage (of taxid chosen earlier) to file size so it's less likely to be filtered later
    bestRatio_df = (
        df
        .sort_values("coverage_megabyte_ratio", ascending=False)
        .drop_duplicates(subset="study_accession", keep="first")
    )

    # filter out the rows with no terms, since they can't be used in downstream statistical analysis anyway
    bestRatio_df = bestRatio_df[bestRatio_df[terms_colname].notnull()]
    bestRatio_df.to_csv(f"{outdir}/suggested_runs_metadata.tsv", sep='\t', index=False)
    bestRatio_df['run_accession'].to_csv(f"{outdir}/suggested_runs.txt", sep='\t', index=False, header=False) # this can be directly plugged into the download pipeline if desired

    # also calculate a range for disk space needed to run analysis on this set of suggested runs
    with open(f"{outdir}/estimated_space_required.txt", "w") as file:
        file.write("Note: runs are suggested based on the assumption that the user is performing statistical comparisons of mutations between environmental sources.\n")
        file.write("As such, only one run is selected from each study. Runs are selected based on best ratio of coverage (of the taxid specified earlier) to file size.\n")
        prefetch_sum = sum(bestRatio_df['size_megabytes']) / 1000 # convert MB to GB
        file.write(f"The total size of prefetch files for this selection of runs is {prefetch_sum:.2f} GB.\n")
        file.write(f"The actual run files are ~7 times larger than the prefetch files, so these will require about {7*prefetch_sum:.2f} GB.\n")
        file.write(f"Due to inevitable redundancies in file size, the total disk space requirement will likely be in the range of {20*prefetch_sum:.2f} GB to {30*prefetch_sum:.2f} GB.\n")
        file.write("Please use nextflow clean after you are done running the pipeline.\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to filter the passedWithMetadata table from the coverage check to a set of suggested runs, and to estimate the disk space required to analyze this set of runs.")

    parser.add_argument("-f", "--file", type=str, help='passedWithMetadata file from coverage check.')
    parser.add_argument("-o", "--outdir", type=str, default=".", help="Output dir.")
    parser.add_argument("-t", "--terms_colname", type=str, default="terms_logistic_regression", help="Name of column for predicted terms.")

    args = parser.parse_args()
    print(args.file, args.outdir, args.terms_colname)
    suggest_runs(args.file, args.outdir, args.terms_colname)