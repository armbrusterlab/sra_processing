#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())

run_tests <- function(breseq_mutants_file, predictions_file, outdir, terms_colname = "terms_logistic_regression", adjust = "fdr", report_all = FALSE) {
    breseq_mutants=read.csv(breseq_mutants_file, header=TRUE, sep='\t')
    predictions=read.csv(predictions_file, header=TRUE, sep='\t') |>
        select(run_accession, !!sym(terms_colname))

    # pivots longer from the items in the terms_colname column (which are comma-separated lists)
    # and renames the column to "term"
    predictions_long <- predictions |>
        rename(term = !!sym(terms_colname)) |>
        mutate(
            term = na_if(trimws(term), "")
        ) |>
        separate_longer_delim(
            term,
            delim = ", "
        )
    
    df <- inner_join(breseq_mutants, predictions_long, by=c("source" = "run_accession"), relationship = "many-to-many")

    # before running any tests, write a basic summary table
    df |> 
        group_by(seq.id, gene, annotation, term) |> summarize(n=n()) |>
        write.table(file.path(outdir, "term_counts_by_mutant.tsv"), sep='\t', row.names = FALSE, quote = FALSE, na = "")

    # another possible summary:  df |> group_by(seq.id, gene, annotation) |> summarize(n=n(), n_terms=(length(unique(term))))

    # TODO use tryCatch on both test functions, where the catch is to write an error message
}

# ### Process CLIs
# args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

# if (length(args) < 3) {
#   stop('Please provide three arguments: <mutation_frequencies_file> <mutant> <outname>')
# }

# mutation_frequencies_file <- args[1]
# mutant <- args[2]
# outname <- args[3]

# envsource_boxplots(mutation_frequencies_file, mutant, outname)