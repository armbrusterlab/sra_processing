#!/usr/bin/env python

import pandas as pd
import json
import re
import argparse

def override(predictions_file, outname,
             keywords_add_file = "", keywords_subtract_file = "", 
             terms_colname = "terms_logistic_regression", category_colname = "category_logistic_regression", subcategory_colname = "subcategory_logistic_regression"):
    # the default values for keywords_add_file and keywords_subtract_file are "" because it's easier to set it up that way with Nextflow
    df = pd.read_csv(predictions_file, sep="\t")

    if keywords_add_file != "":
        keywords_add = pd.read_csv(keywords_add_file, sep="\t")
    if keywords_subtract_file != "":
        keywords_subtract = pd.read_csv(keywords_subtract_file, sep="\t")

    # parse lists from string
    df[category_colname] = [json.loads(s.replace("'", '"')) for s in df[category_colname]]
    df[subcategory_colname] = [json.loads(s.replace("'", '"')) for s in df[subcategory_colname]]

    # first subtract based on keywords
    if keywords_subtract_file != "":
        print("Subtracting terms...")

        for i,row in df.iterrows():
            t = row['text_for_prediction']
            i_to_remove = []
            # print(t)

            for n,kw in enumerate(keywords_subtract['keyword']):
                kw_cat = keywords_subtract['category'][n]
                kw_subcat = keywords_subtract['subcategory'][n]
                # print(kw, kw_cat, kw_subcat)

                if re.search(kw, t.lower()): # if a keyword match is detected, need to remove the corresponding items from category and subcategory
                    # print("Keyword match:", kw)
                    for j in range(len(row[category_colname])):
                        c = row[category_colname][j]
                        s = row[subcategory_colname][j]
                        # print(i, j, c, ":", s)

                        if (kw_subcat == "*" or s == kw_subcat) and c == kw_cat:
                            i_to_remove.append(j)
                            # print("Removing ", kw_cat, kw_subcat)
            df.at[i, category_colname] = [row[category_colname][i] for i in range(len(row[category_colname])) if i not in i_to_remove]
            df.at[i, subcategory_colname] = [row[subcategory_colname][i] for i in range(len(row[category_colname])) if i not in i_to_remove]
            # df.at instead of df.loc, otherwise it will complain that the length of the iterable being assigned doesn't match

    # then add based on keywords (if the subcategory isn't a duplicate)
    if keywords_add_file != "":
        print("Adding terms...")

        for i,row in df.iterrows():
            t = row['text_for_prediction']

            for n,kw in enumerate(keywords_add['keyword']):
                if re.search(kw, t.lower()):
                    kw_cat = keywords_add['category'][n]
                    kw_subcat = keywords_add['subcategory'][n]

                    if kw_subcat not in row[subcategory_colname]:
                        df.at[i, category_colname].append(kw_cat)
                        df.at[i, subcategory_colname].append(kw_subcat)

    # finally, produce the terms_colname column anew based on updated values of category and subcategory columns
    df[terms_colname] = [", ".join(
            [df[category_colname][i][j] + ": " + df[subcategory_colname][i][j] for j in range(len(df[category_colname][i]))]
        ) for i in range(len(df))]

    pd.DataFrame(df).to_csv(outname, sep="\t", index=False) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to edit predicted environmental sources using regex string matching.")

    parser.add_argument("-p", "--predictions_file", type=str, help='Output of predict_environmental_source.py.')
    parser.add_argument("-o", "--outname", type=str, help="Output filename.")

    # not strictly required for the function to run
    parser.add_argument("-a", "--keywords_add_file", type=str, default="", help="Table of keywords for terms to add.")
    parser.add_argument("-m", "--keywords_subtract_file", type=str, default="", help="Table of keywords for terms to subtract.")

    # optional; most likely will not have to modify
    parser.add_argument("-t", "--terms_colname", type=str, default="terms_logistic_regression", help="Terms column name.")
    parser.add_argument("-c", "--category_colname", type=str, default="category_logistic_regression", help="Category column name.")
    parser.add_argument("-s", "--subcategory_colname", type=str, default="subcategory_logistic_regression", help="Subategory column name.")

    args = parser.parse_args()
    override(args.predictions_file, args.outname, 
             keywords_add_file=args.keywords_add_file, keywords_subtract_file=args.keywords_subtract_file, 
             terms_colname=args.terms_colname, category_colname=args.category_colname, subcategory_colname=args.subcategory_colname)

