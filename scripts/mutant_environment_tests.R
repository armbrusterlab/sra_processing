#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)

run_tests <- function(breseq_mutants_file, predictions_file, outdir, terms_colname = "terms_logistic_regression", adjust = "fdr", report_all = TRUE, adjust_separately = FALSE) {
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
    
    df <- inner_join(breseq_mutants, predictions_long, by=c("source" = "run_accession"), relationship = "many-to-many") |>
        mutate(
            annotation = if_else(
                str_detect(annotation, "^([A-Za-z])\\d+([A-Za-z])\\s") & !str_detect(annotation, "^([A-Za-z])\\d+\\1\\s"),
                str_remove(annotation, "\\s.*$"),
                annotation
            )
        ) |> # for non-synonymous mutations, drop the part specifying the exact nucleotide mutated; more concerned with the amino acid change
        mutate(annotation = ifelse(grepl("coding", annotation), paste(mutation, annotation, sep=": "), annotation)) # annotations that include "coding" require more info from the mutation column which would be too specific for other annotations
        # |>
        # mutate(term = tidyr::replace_na(term, "no term")) 

    # before running any tests, write a basic summary table
    summary_by_mutant <- df |> 
        group_by(seq_id, gene, annotation, term) |> 
        summarize(n=n(), .groups = "drop") |>
        complete(nesting(seq_id, gene, annotation), term, fill = list(n = 0)) # for combinations of mutations and terms not observed, set value of 0 
        
    summary_by_mutant |>
        write.table(file.path(outdir, "term_counts_by_mutation.tsv"), sep='\t', row.names = FALSE, quote = FALSE, na = "")

    con <- file(file.path(outdir, "UTF-16LE_term_counts_by_mutation.tsv"), open = "wb")
    writeBin(as.raw(c(0xFF, 0xFE)), con)  # UTF-16LE BOM
    close(con)

    write.table(
        summary_by_mutant,
        file.path(outdir, "UTF-16LE_term_counts_by_mutation.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = "",
        fileEncoding = "UTF-16LE",
        append = TRUE
    )

    # the below is an even more general summary, on the level of genes
    summary_by_gene <- df |>
        group_by(seq_id, term) |> 
        summarize(n=n(), .groups = "drop") |>
        complete(seq_id, term, fill = list(n = 0)) # for combinations of genes and terms not observed, set value of 0 

    summary_by_gene |> 
        write.table(file.path(outdir, "term_counts_by_gene.tsv"), sep='\t', row.names = FALSE, quote = FALSE, na = "")

    con <- file(file.path(outdir, "UTF-16LE_term_counts_by_gene.tsv"), open = "wb")
    writeBin(as.raw(c(0xFF, 0xFE)), con)  # UTF-16LE BOM
    close(con)

    write.table(
        summary_by_gene,
        file.path(outdir, "UTF-16LE_term_counts_by_gene.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = "",
        fileEncoding = "UTF-16LE",
        append = TRUE
    )

    # another possible summary:  df |> group_by(seq_id, gene, annotation) |> summarize(n=n(), n_terms=(length(unique(term))))

    ### Run statistical tests
    df <- df |>
        filter(!is.na(term)) # NA's would cause issues with statistical testing, and it doesn't make sense to compare these anyway

    tryCatch(
        {
            results <- test_by_gene(df, outdir, adjust = adjust, report_all = report_all) |>
                filter(significant == TRUE)

            # new: downstream, to reduce the number of comparisons made, only run stat tests for mutants from genes which were
            # discovered to be enriched in any terms
            df <- df |> filter(seq_id %in% results$seq_id)
        }, 
        error = function(e) {
            writeLines(paste("Error:", e$message), file.path(outdir, "test_by_gene_error.txt"))
    })

    tryCatch(
        {   # now testing only on genes which passed the first test, as explained above
            results <- test_by_mutation(df, outdir, adjust = adjust, report_all = report_all, adjust_separately = adjust_separately)
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

test_by_mutation <- function(df, outdir, adjust="fdr", alpha=0.05, report_all = TRUE, adjust_separately = FALSE) {
    df <- df |> 
        mutate(mutation = paste(seq_id, gene, annotation, sep="; "))
    if (adjust_separately) {
        subdir = file.path(outdir, "stats_by_mutation")
        if (!dir.exists(subdir)) {
            dir.create(subdir, recursive = TRUE)
        }

        df |> 
            group_by(seq_id) |> 
            group_map(function(group_data, group_keys) {
                curr_seq = group_keys$seq_id
                cat(curr_seq)
                cat(group_data |> nrow())
                run_fisher_2x2_tests(group_data, "mutation", subdir, paste0("mutations_", curr_seq, ".tsv"), adjust, alpha, report_all)
        })
    } else {
        run_fisher_2x2_tests(df, "mutation", outdir, "stats_by_mutation.tsv", adjust, alpha, report_all)
    }
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
                    ),
                    alternative = "greater"
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
        write.table(file.path(outdir, basename), sep='\t', row.names = FALSE, quote = FALSE, na = "")

    # results |> 
    #     write.table(file.path(outdir, paste0("UTF-16LE_", basename)), sep='\t', row.names = FALSE, quote = FALSE, na = "", fileEncoding = "UTF-8")

    con <- file(file.path(outdir, paste0("UTF-16LE_", basename)), open = "wb")

    writeBin(as.raw(c(0xFF, 0xFE)), con)  # UTF-16LE BOM

    close(con)

    write.table(
        results,
        file.path(outdir, paste0("UTF-16LE_", basename)),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = "",
        fileEncoding = "UTF-16LE",
        append = TRUE
    )

    return(results)
}

### Process CLIs (from Nextflow)
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 7) {
  stop('Not enough arguments provided')
}

breseq_mutants_file <- args[1]
predictions_file <- args[2]
outdir <- args[3]
terms_colname <- args[4]
adjust <- args[5]
report_all <- args[6] == "TRUE" # it's read as a string from the command line
adjust_separately <- args[7] == "TRUE" # it's read as a string from the command line

run_tests(breseq_mutants_file, predictions_file, outdir, terms_colname, adjust, report_all, adjust_separately)