library(tidyverse)

trust4_files <- list.files('/public/ClinicalExam/lj_sih/projects/project_pipeline/RNA/results/IG/TRUST4','*report.tsv',full.name = T)

BCR_result <- list()
for (i in seq_len(length(trust4_files))){
    file <- trust4_files[i]
    BCR <- read_tsv(file)
    sample_name <- basename(file) %>% str_remove('_report.tsv')
    BCR$sample_name <- sample_name
    BCR_result[[i]] <- BCR
}

ig_data <- reduce(BCR_result,rbind)

# ig_data <- readRDS('results/Ruijin/RNA/IG.rds')

ig_data_top <- ig_data %>% group_by(sample_name) %>%
  top_n(5,wt = `#count`) %>%  mutate(clonetype = paste0('Clone_',row_number())) %>%
  mutate(cid = fct_lump(cid, n = 2,w = `#count`))

colnames(ig_data_top)[1] <- 'count'

ig_data_top <- ig_data_top %>% mutate(frequency = count/sum(count))
write_csv(ig_data_top,'RNA_BCR_frac.csv')