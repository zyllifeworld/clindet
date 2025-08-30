library(tidyverse)
library(glue)

input_cnv <- snakemake@input[['cnv']]
output_ascat <- snakemake@output[['ascat']]
brass_cnv <- snakemake@params[['brass_cnv']]
if (brass_cnv == 'ASCAT'){
    load(input_cnv)
    sink(output_ascat,append=FALSE,split = FALSE)
    glue('rho ',res$summary$allSols.refitted$purity)
    glue('Ploidy ',res$summary$allSols.refitted$ploidy)
    glue('GenderChr X')
    glue('GenderChrFound N')
    sink()
} else if (brass_cnv == 'purple') {
    df <- read_tsv(input_cnv)
    gender <- if_else(df$gender == 'MALE','Y','X')
    sink(output_ascat,append=FALSE,split = FALSE)
    glue('rho ',df$purity[1])
    glue('Ploidy ',df$ploidy[1])
    glue('GenderChr {gender}')
    glue('GenderChrFound N')
    sink()
}
