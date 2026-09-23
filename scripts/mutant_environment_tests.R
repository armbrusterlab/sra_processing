#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

run_tests <- function(breseq_mutants_file, predictions_file, outdir, terms_colname = "terms_logistic_regression", adjust = "fdr", report_all = TRUE) {
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

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
    
    df <- inner_join(breseq_mutants, predictions_long, by=c("source" = "run_accession"), relationship = "many-to-many") # |>
        # mutate(term = tidyr::replace_na(term, "no term")) 

    # before running any tests, write a basic summary table
    # for now I am opting not to do per-mutation tests because that is a lot of comparisons; the user can use this table to run their own tests if they want
    df |> 
        group_by(seq_id, gene, annotation, term) |> summarize(n=n()) |>
        write.table(file.path(outdir, "term_counts_by_mutation.tsv"), sep='\t', row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-16")

    # the below is an even more general summary, on the level of genes
    df |>
        group_by(seq_id, term) |> summarize(n=n()) |> 
        write.table(file.path(outdir, "term_counts_by_gene.tsv"), sep='\t', row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-16")

    # another possible summary:  df |> group_by(seq_id, gene, annotation) |> summarize(n=n(), n_terms=(length(unique(term))))

    ### Run statistical tests
    df <- df |>
        filter(!is.na(term)) # NA's would cause issues with statistical testing, and it doesn't make sense to compare these anyway

    tryCatch(
        {
        test_by_gene(df, outdir, adjust = adjust, report_all = report_all)
        }, 
        error = function(e) {
            writeLines(paste("Error:", e$message), file.path(outdir, "test_by_gene_error.txt"))
    })

    tryCatch(
        {
        test_by_mutation(df, outdir, adjust = adjust, report_all = report_all)
        }, 
        error = function(e) {
            writeLines(paste("Error:", e$message), file.path(outdir, "test_by_mutation_error.txt"))
    })
    # test_by_gene(df, outdir, adjust = adjust, report_all = report_all)
    # test_by_mutation(df, outdir, adjust = adjust, report_all = report_all)
}

test_by_gene <- function(df, outdir, adjust="fdr", alpha=0.05, report_all = TRUE) {
    run_fisher_2x2_tests(df, "seq_id", outdir, "stats_by_gene.tsv", adjust, alpha, report_all)
}

test_by_mutation <- function(df, outdir, adjust="fdr", alpha=0.05, report_all = TRUE) {
    df <- df |> 
        mutate(mutation = paste(seq_id, gene, annotation, sep="; "))
    
    run_fisher_2x2_tests(df, "mutation", outdir, "stats_by_mutation.tsv", adjust, alpha, report_all)
}

# this function was written with the assistance of Copilot
run_fisher_2x2_tests <- function(df, data_col, outdir, basename, adjust="fdr", alpha=0.05, report_all = TRUE) {
    data_values <- sort(unique(df[[data_col]]))
    terms <- sort(unique(df$term))

    results <- expand.grid(
        tempname = data_values,
        term = terms,
        stringsAsFactors = FALSE
    ) |>
        rowwise() |>
        mutate(
            p.value = {
                current_value <- tempname
                current_term <- term

                is_value <- df[[data_col]] == current_value
                is_term <- df$term == current_term

                fisher.test(
                    matrix(
                        c(
                            sum(is_value & is_term),
                            sum(is_value & !is_term),
                            sum(!is_value & is_term),
                            sum(!is_value & !is_term)
                        ),
                        nrow = 2
                    )
                )$p.value
            }
        ) |>
        ungroup() |>
        mutate(
            p.adj = p.adjust(p.value, method = adjust),
            significant = (p.adj < alpha)
        )

     names(results)[names(results) == "tempname"] <- data_col # rename tempname to the actual name of the column analyzed

    if (!report_all) {
        results <- results |>
            filter(significant)
    }
    
    results |> 
        write.table(file.path(outdir, basename), sep='\t', row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-16")
}

### Process CLIs (from Nextflow)
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 6) {
  stop('Not enough arguments provided')
}

breseq_mutants_file <- args[1]
predictions_file <- args[2]
outdir <- args[3]
terms_colname <- args[4]
adjust <- args[5]
report_all <- args[6] == "TRUE" # it's read as a string from the command line

run_tests(breseq_mutants_file, predictions_file, outdir, terms_colname, adjust, report_all)