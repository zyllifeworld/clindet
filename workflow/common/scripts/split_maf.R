library(readr)
library(dplyr)

input_maf <- snakemake@inputs[['maf']]
output_snp <- snakemake@inputs[['snp']]
output_indel <- snakemake@inputs[['indel']]

lines <- readLines(if (grepl("\\.gz$", input_maf)) gzfile(input_maf, "rt") else input_maf)

comment_idx <- grepl("^#", lines)
first_data_idx <- which(!comment_idx)[1]

comment_lines <- if (first_data_idx > 1) lines[1:(first_data_idx - 1)] else character()
data_text <- paste(lines[first_data_idx:length(lines)], collapse = "\n")

maf <- read_tsv(I(data_text), show_col_types = FALSE)

if (!"Variant_Type" %in% colnames(maf)) {
    stop("Variant_Type column not found in MAF header.")
}

maf_snp <- maf %>% filter(Variant_Type == "SNP")
maf_indel <- maf %>% filter(Variant_Type %in% c("INS", "DEL"))

write_with_comments <- function(df, outfile, comments) {
    con <- file(outfile, "wt")
    on.exit(close(con), add = TRUE)
    if (length(comments) > 0) writeLines(comments, con)
    write_tsv(df, con)
}

write_with_comments(maf_snp, output_snp, comment_lines)
write_with_comments(maf_indel, output_indel, comment_lines)