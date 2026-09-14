import pandas as pd

# myfile = r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\isolation_sources_bacdive_2026-06-25.csv" # raw string so I don't have to change the backslashes to forward slashes
myfile = r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\isolation_sources_bacdive_2026-08-12.csv" # raw string so I don't have to change the backslashes to forward slashes
print(f"file is {myfile}")

# with open(myfile, "r") as f:
#     lines = f.readlines()
#     for line in lines[1:]: # skip header
#         t = line.split(",") # this doesn't parse correctly when a cell in the csv contains a comma
#         print(t[0], t[3], t[6], t[7], t[8])

df = pd.read_csv(myfile)

# for i in range(len(df)):
#     line = df.iloc[i]
#     print(line["ID"], line["Isolation source"], line["Category 1"], line["Category 2"], line["Category 3"])
#     print(df.iloc[i]["Isolation source"].isna())

# Combine rows with same 'ID' and concatenate other values
# combined_df = (
#     df.groupby('ID', as_index=False)
#       .agg({
#           "Isolation source": lambda x: list(x)[0], # the rest of the items in the list are NaN's
#           "Category 1": lambda x: '###'.join([s.lstrip('#') for s in x if isinstance(s, str)]),
#           "Category 2": lambda x: '###'.join([s.lstrip('#') for s in x if isinstance(s, str)]),
#           "Category 3": lambda x: '###'.join([s.lstrip('#') for s in x if isinstance(s, str)])
#       })
# )

# actually, I want to join categories 1 and 2, THEN collapse them.
combined_df = (
    df.groupby('ID', as_index=False)
      .agg({
          "Species": lambda x: list(x)[0],
          "Culture collection number": lambda x: list(x)[0],
          "Isolation source": lambda x: list(x)[0], # the rest of the items in the list are NaN's
          "Country": lambda x: list(x)[0],
          "Continent": lambda x: list(x)[0],
          "Category 1": lambda x: [s.lstrip('#') for s in x if isinstance(s, str)],
          "Category 2": lambda x: [s.lstrip('#') for s in x if isinstance(s, str)],
          "Category 3": lambda x: [s.lstrip('#') for s in x if isinstance(s, str)]
      })
)

combined_df["joined_1_2"] = [
    # [a + "@@@" + b for a, b in zip(row["Category 1"], row["Category 2"])]
    # there are redundant items in the joined_1_2 column so let's make it non-redundant
    list(set([a + "@@@" + b for a, b in zip(row["Category 1"], row["Category 2"])])) 
    for _, row in combined_df.iterrows()
]

# convert lists to reliably-parseable string format
# note that these lists don't necessarily have the same number of elements as each other, nor as the isolation source column
combined_df["Category 1"] = [
    '###'.join(row["Category 1"])
    for _, row in combined_df.iterrows()
]

combined_df["Category 2"] = [
    '###'.join(row["Category 2"])
    for _, row in combined_df.iterrows()
]

combined_df["Category 3"] = [
    '###'.join(row["Category 3"])
    for _, row in combined_df.iterrows()
]

combined_df["joined_1_2"] = [
    '###'.join(row["joined_1_2"])
    for _, row in combined_df.iterrows()
]

len(df) # 150860
len(set(df["ID"])) # 63296
len(combined_df) # 63296

# outfile = r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\bacdiveReformat_2026-06-25.tsv"
outfile = r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\bacdiveReformat_2026-08-12.tsv"
combined_df.to_csv(outfile, sep="\t", header = True, index = False)
# the idea is that the joined_1_2 column will serve as input to the ML model (after being re-separated into two columns, each with lists as elements)