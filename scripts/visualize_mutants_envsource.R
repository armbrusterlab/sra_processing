#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
theme_set(theme_classic())

envsource_boxplots <- function(f, mutant, outname) {
    df <- read.csv(f, header=TRUE, sep="\t")

    df |>
        filter(Mutation == mutant) |>
        ggplot(aes(x=Subcategory, y=Count, fill=Category)) +
        geom_bar(stat = "identity") +
        facet_wrap(vars(Category), scales = "free_x") +
        labs(title=paste0("Box plots of environmental source labels for ", {mutant})) +
        theme(legend.position = "none")
    
    ggsave(outname)
}

### Process CLIs
args <- commandArgs(trailingOnly = TRUE) # only get the CLIs that come after the name of the script

if (length(args) < 3) {
  stop('Please provide three arguments: <mutation_frequencies_file> <mutant> <outname>')
}

mutation_frequencies_file <- args[1]
mutant <- args[2]
outname <- args[3]

envsource_boxplots(mutation_frequencies_file, mutant, outname)