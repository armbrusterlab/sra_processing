import pandas as pd
import argparse

def join_data(f1, f2, f3, outname):
    df1=pd.read_csv(f1, sep='\t') # breseq output
    df2=pd.read_csv(f2, sep=',') # NCBI esearch metadata
    df3=pd.read_csv(f3, sep='\t') # pysradb metadata

    df_full = pd.merge(df1, df2, left_on = "source", right_on = "Run").drop(columns="source")
    df_full = pd.merge(df_full, df3, left_on = "Run", right_on = "run_accession").drop(columns="Run")

    # reorganize df for legibility
    # move important columns to the front
    for c in ["evidence_file", "run_accession"]:
        col = df_full.pop(c)
        df_full.insert(0, col.name, col)

    # replace ambiguous column names
    df_full.rename(columns={"avgLength":"avg_read_length", "size_MB": "size_megabytes"}, inplace=True)

    df_full.to_csv(outname, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to summarize breseq outputs.")

    # positional arguments (required)
    parser.add_argument("breseq_summary", type=str, help='Output of summarize_breseq.py.')
    parser.add_argument("esearch_metadata", type=str, help="File containing SRA metadata from NCBI esearch.")
    parser.add_argument("pysradb_metadata", type=str, help="File containing SRA metadata from pysradb.")
    parser.add_argument("outname", type=str, help="The name of the file to write the results to.")

    args = parser.parse_args()
    join_data(args.breseq_summary, args.esearch_metadata, args.pysradb_metadata, args.outname)