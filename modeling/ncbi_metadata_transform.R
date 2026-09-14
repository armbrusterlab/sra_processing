library(tidyverse)

transform_metadata <- function(metadata_file, outname) {
  # first flatten it to one row per record
  df=read.csv(metadata_file, header=TRUE, sep="\t", na.strings=c("","NA")) # note that this enters periods in place of spaces in the colnames, and turns the blank strings to NA values
  df2 <- df |> 
    pivot_wider(
      names_from="Assembly.BioSample.Attribute.Name",
      values_from="Assembly.BioSample.Attribute.Value",
      values_fn = \(x) paste(unique(x), collapse = "; ")
    )

  # for reference, this is how you'd filter to rows where a column, containing a space in its name, doesn't have an NA: you wrap the column name in backticks
  # df2 |> filter(!is.na(`Isolation Site`))

  # extract the info from the titles
  # this should handle NA's gracefully, and if I later decide I don't want to use any one of these columns it's easy to remove it
  title_cols = c("Assembly.BioProject.Lineage.Title", "Assembly.BioSample.BioProject.Title") # for now remove Assembly.BioSample.Description.Title
  titles <- df2 |> 
    select(all_of(title_cols)) |> 
    mutate(across(title_cols, ~ str_extract(.x, "(?<=from ).*") |> 
    str_replace_na(replacement = "")))

  # for convenience, grab the other useful columns from df2
  # use any_of() to handle it gracefully if there ends up being no "Isolation Site" column (as it is pivoted out of other metadata and may not exist in the dataset the user obtains
  # additionally, turn NA's to blank spaces here as well so they don't show up in the string
  df_useful <- df2 |> 
    select(any_of(c("Assembly.BioSample.Isolation.source", "Assembly.BioSample.Host.disease", "Isolation Site"))) |> 
    mutate(across(everything(), ~ str_replace_na(.x, replacement = "")))

  # join the strings and simultaneously add this new column to df
  # can't use "df2$foo <-" notation because the output of unite already has a name that you have to provide, and things get weird when you do that
  df2 <- cbind(df2, unite(cbind(titles, df_useful), "joined_string", sep=" "))

  # if for whatever reason an accession is duplicated, assume the rows are all duplicates, and keep only the first instance
  df2 <- df2 |>
    group_by(.data[["Assembly.Accession"]]) |> 
    slice(1) |>
    ungroup()

  # save to file
  df2 |> 
    write.table(file=outname, quote=FALSE, sep='\t', row.names=FALSE) # rownames false so that it doesn't add a 1, 2, 3... to the left of the table
}

### Process CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 2) {
  stop('Please provide two arguments: <metadata_file> <outname>')
}

metadata_file <- args[1]
outname <- args[2]

cat("Metadata file provided:", metadata_file, "\n")
cat("Output file:", outname, "\n")

transform_metadata(metadata_file, outname)