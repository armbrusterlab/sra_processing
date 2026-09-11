#!/usr/bin/env Rscript
# use this script to trim all short paired reads to uniform length for downstream use with GangSTR
# example run:
# Rscript "/home/kcw2/metadata-magnet/scripts/sra_processing/trim_reads.R" -d "/home/kcw2/data/testing/sra_wbpI_fha1/cf/grepq/paired/" -o "/home/kcw2/data/testing/sra_wbpI_fha1/cf/trimmedReads"

library(ShortRead)

trim_reads <- function(r1, r2, outdir, minWidth=NA) {
    fq1 = readFastq(r1)
    fq2 = readFastq(r2)

    print("Summary of input data:")
    print(summary(width(fq1)))
    cat(length(fq1), "reads in R1\n")
    print(summary(width(fq2)))
    cat(length(fq2), "reads in R2\n")

    if (is.na(minWidth)) {
        minWidth <- trunc(min(summary(width(fq1))["Mean"], summary(width(fq2))["Mean"]))
        cat("Inferring minWidth; using", minWidth)
    }

    # remove reads shorter than minWidth
    my_filter <- srFilter(function(x) width(x) >= minWidth, name = "MinLength")
    fq1_filtered <- fq1[my_filter(fq1)]
    fq2_filtered <- fq2[my_filter(fq2)]
    cat("After filtering reads shorter than minWidth=", minWidth, ":\n")
    print(length(fq1_filtered))
    print(length(fq2_filtered))

    # crop the remaining reads
    fq1_cropped <- narrow(fq1_filtered, start=1, end=minWidth)
    fq2_cropped <- narrow(fq2_filtered, start=1, end=minWidth)

    # remove reads that don't match between the two datasets
    # headers have been modified by kraken2, so if the two reads in the pair don't have the exact same kraken mapping they won't be identical.
    # Must remove everything past the name of the read.
    new_ids1 <- BStringSet(gsub(" .*", "", id(fq1_cropped)))
    fq1_renamed <- ShortReadQ(sread(fq1_cropped), quality(fq1_cropped), new_ids1)
    new_ids2 <- BStringSet(gsub(" .*", "", id(fq2_cropped)))
    fq2_renamed <- ShortReadQ(sread(fq2_cropped), quality(fq2_cropped), new_ids2)

    intersected <- intersect(id(fq1_renamed), id(fq2_renamed))
    fq1_id <- which(id(fq1_renamed) %in% intersected)
    fq2_id <- which(id(fq2_renamed) %in% intersected)

    fq1_final <- fq1_renamed[fq1_id]
    fq2_final <- fq2_renamed[fq2_id]

    print("After filtering to reads that appear in both mates:")
    print(length(fq1_final))
    print(length(fq2_final))

    # save outputs
    print("Saving outputs...")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    writeFastq(fq1_final, file.path(outdir, "r1.fq.gz"), compress = TRUE)
    writeFastq(fq2_final, file.path(outdir, "r2.fq.gz"), compress = TRUE)
}

library(optparse)
# option_list <- list(
#     make_option(c("-a", "--read1"), type = "character", help = "name of read 1 fastq file (required)"),
#     make_option(c("-b", "--read2"), type = "character", help = "name of read 1 fastq file (required)"),
#     make_option(c("-o", "--outdir"), type = "character", help = "output directory name (required)"),
#     make_option(c("-w", "--minWidth"), type = "numeric", default = 150, help = "threshold below which reads should be discarded")
# )
# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)
# if (is.null(opt$read1) | is.null(opt$read2) | is.null(opt$outdir)) {
#  print_help(opt_parser)
# }
# trim_reads(opt$read1, opt$read2, opt$outdir, minWidth = opt$minWidth)

# might want to instead take as input an entire dir of reads (or list of them) and run trim_reads on each so we only have to load ShortRead once 
option_list <- list(
    make_option(c("-d", "--fqdir"), type = "character", help = "name of dir in which every subdir contains *_1.fq.gz and *_2.fq.gz paired reads (required)"),
    make_option(c("-o", "--outdir"), type = "character", help = "output directory name (required)"),
    make_option(c("-w", "--minWidth"), type = "numeric", default = NA, help = "threshold below which reads should be discarded")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$fqdir) | is.null(opt$outdir)) {
 print_help(opt_parser)
}

dirs <- Sys.glob(paste(opt$fqdir, "/*/", sep=""))
for (dir in dirs) {
    run <- basename(dir)
    cat("Processing", run)
    out <- file.path(opt$outdir, run)
    r1 <- Sys.glob(file.path(dir, "*_1.fq.gz"))
    r2 <- Sys.glob(file.path(dir, "*_2.fq.gz"))

    trim_reads(r1, r2, out, minWidth = opt$minWidth)
}