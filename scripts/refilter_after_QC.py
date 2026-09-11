from pathlib import Path
import pandas as pd
import json
import argparse

# After QC (fastp/fastplong) there may not be sufficient reads left in a dataset to run breseq downstream.
# This script estimates the coverage for the taxid of choice (set by user earlier in the pipeline)
# and returns lists of run ID files after filtering 

def refilter_after_qc(qcdir, runidsdir, metadata_file, genome_length, coverage_threshold):
    # process the PE, SE short, and SE long data dirs separately, and save to the current dir an updated ID file for each
    idfiles = [f"{runidsdir}/ids_paired_filtered.txt", f"{runidsdir}/ids_single_long_filtered.txt", f"{runidsdir}/ids_single_short_filtered.txt"]
    reportdirs = [f"{qcdir}/paired/reports", f"{qcdir}/single/long/reports", f"{qcdir}/single/short/reports"]
    isPaired = [True, False, False]
    outnames = ["run_ids_postqc/ids_paired_postQC.txt", "run_ids_postqc/ids_single_long_postQC.txt", "run_ids_postqc/ids_single_short_postQC.txt"]

    for input in zip(idfiles, reportdirs, isPaired, outnames):
        check_coverage(input[0], input[1], input[2], input[3], metadata_file, genome_length, coverage_threshold)

def check_coverage(idlist, reportdir, paired, outname, metadata_file, genome_length, coverage_threshold):
    metadata = pd.read_csv(metadata_file, sep="\t") # metadata_file is the passedWithMetadata output of screen_sra_coverage.py

    with open(outname, 'w') as outfile:
        with open(idlist, 'r') as ids:
            print(f"Writing to {outname}")
            # iterate thru IDs. Find the corresponding QC report in reportdir (ifelse based on paired). Load it.
            for line in ids:
                run_id = line.strip()
                reportfile = f"{reportdir}/{run_id}pe.fastq.json" if paired else f"{reportdir}/{run_id}.fastq.json"
                # print(f"Run id {run_id}, reportfile {reportfile}")

                # load the json
                with open(reportfile, "r") as f:
                    raw = f.read().strip()

                # if empty file, skip
                if not raw:
                    print(f"WARNING: empty JSON for {run_id}")
                    bp_remaining = 0

                # if invalid json, skip
                try:
                    data = json.loads(raw)
                    bp_remaining = data["summary"]["after_filtering"]["total_bases"]
                except json.JSONDecodeError:
                    print(f"WARNING: invalid JSON for {run_id}")
                    bp_remaining = 0

                # calculate the coverage remaining for this run ID after filtering
                # to estimate the proportion of the total bases that count toward coverage for the taxid of interest, get info from metadata file
                row = metadata[metadata["run_accession"] == run_id].iloc[0] # there should be only one row that matches the run accession anyway
                taxid_proportion = row["spots_for_taxid"] / row["total_spots"]
                coverage = (taxid_proportion * bp_remaining) / genome_length
                passed = coverage >= coverage_threshold


                # finally, calculate coverage 
                coverage = (taxid_proportion * bp_remaining) / genome_length
                passed = coverage >= coverage_threshold
                print(f"run_id {run_id}\tbp {bp_remaining}\tproportion of bp {taxid_proportion}\tcoverage {coverage}\tpassed? {passed}") # diagnostics

                if passed:
                    outfile.write(f"{run_id}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to re-filter SRA runs after QC, removing runs with insufficient estimated coverage.")

    parser.add_argument("-q", "--qcdir", type=str, help='Directory containing fastp/fastplong outputs organized as expected by the pipeline.')
    parser.add_argument("-i", "--runidsdir", type=str, help='Directory containing files with runids.')
    parser.add_argument("-m", "--metadata", type=str, help='Metadata file; the passedWithMetadata output of screen_sra_coverage.py.')
    parser.add_argument("-L", "--genome_length", type=float, help="Length of reference genome in bp.")
    parser.add_argument("-c", "--coverage_threshold", type=float, default=30, help="Coverage required for the taxonomic ID.") # optional

    args = parser.parse_args()
    print(args.qcdir, args.runidsdir, args.metadata, args.genome_length, args.coverage_threshold)
    refilter_after_qc(args.qcdir, args.runidsdir, args.metadata, args.genome_length, args.coverage_threshold)