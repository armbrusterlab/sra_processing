#!/usr/bin/env python

import pandas as pd
import json
import argparse

def get_sources_per_mutant(f, outname = "mutation_frequencies.tsv", category_colname = "category_logistic_regression", subcategory_colname = "subcategory_logistic_regression"):
    # f is breseq_summary_withMetadata.tsv
    # category_colname and subcategory_colname have default values according to the models in the v19 and v20 dirs,
    # but if you use your own model and the model has a different name, please change these values accordingly.
    df = pd.read_csv(f, sep="\t")

    # mutation labels: mutate a column that joins gene and description
    # id_cols = ["gene", "annotation"]
    id_cols = ["seq_id", "annotation"] # I debated whether to also join with the mutation column but decided against it because different mutations may be synonymous

    df['mutation_id'] = (
        df[id_cols]
        .fillna('') # NA's become empty strings
        .agg(': '.join, axis=1)  # join the strings
    )

    df['mutation_id'] = [s.replace(u'\xa0', ' ') for s in df['mutation_id']] # \xa0 is a non-breaking space, but for this purpose a regular space will do

    # parse lists from string
    df[category_colname] = [json.loads(s.replace("'", '"')) for s in df[category_colname]]
    df[subcategory_colname] = [json.loads(s.replace("'", '"')) for s in df[subcategory_colname]]

    freqs = {id : {} for id in set(df['mutation_id'])}

    for _, row in df.iterrows():
        id = row['mutation_id']
        for i in range(len(row[category_colname])): 
            k = (row[category_colname][i], row[subcategory_colname][i])
            freqs[id][k] = freqs[id].get(k, 0) + 1

    # convert the nested dictionary to df
    data = {'Mutation' : [], 'Category' : [], 'Subcategory' : [], 'Count' : []}
    for id in freqs.keys():
        for label in freqs[id].keys():
            data['Mutation'].append(id)
            data['Category'].append(label[0])
            data['Subcategory'].append(label[1])
            data['Count'].append(freqs[id][label])

    pd.DataFrame(data).to_csv(outname, sep="\t", index=False) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to count environmental sources per mutation.")

    parser.add_argument("-f", "--mutants_file", type=str, help='breseq_summary_withMetadata.tsv containing environmental sources predicted from metadata.')
    parser.add_argument("-o", "--outname", type=str, default="mutation_frequencies.tsv", help="Output filename.")
    parser.add_argument("-c", "--category_colname", type=str, default="category_logistic_regression", help="Name of column for predicted category.")
    parser.add_argument("-s", "--subcategory_colname", type=str, default="subcategory_logistic_regression", help="Name of column for predicted subcategory.")

    args = parser.parse_args()
    print(args.mutants_file, args.outname, args.category_colname, args.subcategory_colname)
    get_sources_per_mutant(args.mutants_file, args.outname, args.category_colname, args.subcategory_colname)