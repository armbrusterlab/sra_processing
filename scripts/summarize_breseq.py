#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd
import re

def write_summary(breseq_dirs, outdir, n=1, filter_intergenic = False, filter_synonymous = False):
    outdir_path = Path(outdir)
    outdir_path.mkdir(exist_ok=True, parents=True)

    # it's generally bad practice to grow a df, but this will give flexibility wrt columns of the output df
    df_mutations = pd.DataFrame() 
    df_missingCoverage = pd.DataFrame() 
    df_newJunction = pd.DataFrame() 
    with open(breseq_dirs) as f:
        for line in f:
            dir = line.strip()
            print(f"Processing {dir}")

            # dir_path = Path(dir)
            # for item in dir_path.iterdir():
            #     print(item)

            summary = f"{dir}/index.html"
            tab = pd.read_html(summary, extract_links="body") # need link locations from the evidence column, but this turns all elements into tuples
            for i in range(1, len(tab)):
                df=tab[i]
                try:
                  table_type = list(df.iloc[-2])[0][0]

                  df.columns = [tup[0] for tup in df.iloc[-1]] # rename df columns so they're referenceable; colnames are in last row
  
                  if table_type == "Predicted mutations":
                      df["evidence_file"] = df["evidence"].apply(lambda x: x[1] if isinstance(x, tuple) else None) # get link locations
  
                  # clean up columns now- no more tuples
                  for col in df.columns:
                      if col != "evidence_file":
                          df[col] = df[col].apply(lambda x: x[0] if isinstance(x, tuple) else x)
  
                  df = df.iloc[:-2] # last two rows don't actually have data, so crop them.

                  if n == -1: # this is used in the Nextflow context. It gets the basename, which is the lowest dir.
                      name = Path(dir).parts[-1]
                  else:
                      name = "/".join(Path(dir).parts[-1 * (n+1):-1])
                  df["source"] = name
                except: # the above will fail if no mutations predicted
                  table_type = None
                    
                if table_type == 'Predicted mutations':
                    df_mutations = pd.concat([df_mutations, df], axis=0, join='outer', ignore_index=True)
                elif table_type == 'Unassigned missing coverage evidence':
                    df_missingCoverage = pd.concat([df_missingCoverage, df], axis=0, join='outer', ignore_index=True)
                elif table_type == 'Unassigned new junction evidence':
                    df_newJunction = pd.concat([df_newJunction, df], axis=0, join='outer', ignore_index=True)

    # now that all dfs have been joined, filter if requested.
    if filter_intergenic:
        print("Filtering out intergenic variants...")
        df_mutations = df_mutations[df_mutations['annotation'].str.contains('intergenic', na=False) == False]

    if filter_synonymous:
        print("Filtering out synonymous mutants...")
        df_mutations = df_mutations[df_mutations['annotation'].apply(find_synonymous) == False]

    # fix non-breaking spaces in mutations file (they appear in other files, but these aren't the priority)

    # for the following: the former string has character U+00a0 (non-breaking space)
    df_mutations = df_mutations.rename(columns={'seq id': 'seq_id'}) 
    df_mutations['annotation'] = df_mutations['annotation'].str.replace(' ', ' ')
    df_mutations['gene'] = df_mutations['gene'].str.replace(' ', ' ')
    df_mutations['mutation'] = df_mutations['mutation'].str.replace(' ', ' ')

    # # save dfs to files with utf-16 encoding to avoid displaying garbled text in Excel
    # df_mutations.to_csv(outdir_path.joinpath("mutations.tsv"), sep="\t", index=False, header=True, encoding="utf-16")
    # df_missingCoverage.to_csv(outdir_path.joinpath("missingCoverage.tsv"), sep="\t", index=False, header=True, encoding="utf-16")
    # df_newJunction.to_csv(outdir_path.joinpath("newJunction.tsv"), sep="\t", index=False, header=True, encoding="utf-16")

    # utf-16 encoding caused irreconcilable issues downstream, so the files will be saved without it
    df_mutations.to_csv(outdir_path.joinpath("mutations.tsv"), sep="\t", index=False, header=True)
    df_missingCoverage.to_csv(outdir_path.joinpath("missingCoverage.tsv"), sep="\t", index=False, header=True)
    df_newJunction.to_csv(outdir_path.joinpath("newJunction.tsv"), sep="\t", index=False, header=True)

    # the mutations df will be saved in a more comprehensive form later, so for now, save the others as xlsx
    df_missingCoverage.to_csv(outdir_path.joinpath("missingCoverage_utf16.tsv"), sep="\t", index=False, header=True, encoding="utf-16")
    df_newJunction.to_csv(outdir_path.joinpath("newJunction_utf16.tsv"), sep="\t", index=False, header=True, encoding="utf-16")

    print("Done!")

def find_synonymous(s):
    if s != s: # this is true only if s is NaN
        return False
    
    pattern = r"^([A-Z])[0-9]+([A-Z])?" # capturing groups: the two nucleotides flanking some number
    found = re.match(pattern, s)
    if found:
        return found.group(1) == found.group(2)
    else: # e.g. if it's an intergenic variant or blank
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to summarize breseq outputs.")

    # positional arguments (required)
    parser.add_argument("breseq_dirs", type=str, help='A file containing one breseq output dir per line. Each output dir must contain a breseq output file "index.html".')
    parser.add_argument("outdir", type=str, help="The name of the output directory to write the results to.")

    # named (optional) arguments
    parser.add_argument("-n", "--num", type=int, default=1, help="Number of levels of dir names to preserve in the output (source column).")
    parser.add_argument("-i", "--filter_intergenic", action="store_true", help="Filter out intergenic variants.")
    parser.add_argument("-s", "--filter_synonymous", action="store_true", help="Filter out synonymous mutations.")

    args = parser.parse_args()
    write_summary(args.breseq_dirs, args.outdir, args.num, args.filter_intergenic, args.filter_synonymous)