#parse args
args <- commandArgs(TRUE)

isocan <- args[1]
input_files <- args[2] #input
original_gtf <- args[3] #input
original_corr <- args[4]
original_fa <- args[5]
filter_res <- args[6] #sqanti result
gtf_out <- args[7] #update out
corr_out <- args[8]
fa_out <- args[9]
count_unique_out <- args[10] #used
existence_out <- args[11]

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(purrr)
  pkgload::load_all(isocan, export_all = FALSE, quiet = TRUE)
})

input_files <- read_csv(input_files, col_types = cols(.default = col_character()))

#filter
message("Filtering for potential SQANTI isoforms...")
selected_isoform <- filter_res %>% read_tsv(show_col_types = FALSE) %>% filter(filter_result == "Isoform") %>% `[[`("isoform")
  
#update gtf
message("Updating transcriptome annotation...")
updated_gtf <- 
  original_gtf %>% 
  isocan::read_gtf() %>% 
  filter(transcript_id %in% selected_isoform)
multi_exon_isoform <- 
  updated_gtf %>% 
  count(transcript_id) %>% 
  filter(n >= 2L) %>% 
  `[[`("transcript_id")
isocan::write_gtf(updated_gtf, gtf_out)

#update corr
corr <- original_corr %>% read_tsv(col_types = cols(.default = col_character())) %>% filter(transcript_id %in% selected_isoform)
write_tsv(corr, corr_out)

unique_corr <- anti_join(corr, corr %>% select(seq_id, sample) %>% mutate(dupli = duplicated(.)) %>% filter(dupli), by = c("seq_id", "sample"))
#write_tsv(get_count(unique_corr), count_unique_out)

message("judge FL existence")
bind_rows(
  corr %>% filter(type == "intronic_match_include" & transcript_id %in% multi_exon_isoform),
  unique_corr %>% filter(!transcript_id %in% multi_exon_isoform)
) %>% 
  mutate(existence = TRUE) %>% 
  group_by(gene_id, transcript_id, sample) %>% 
  summarise(existence = any(existence)) %>% 
  ungroup() %>% 
  spread(key = "sample", value = "existence", fill = FALSE) %>% 
  write_tsv(existence_out)

#update fa
message("Updating transcript FASTA...")
bs_fa <- Biostrings::readDNAStringSet(original_fa)
bs_fa[names(bs_fa) %in% selected_isoform] %>% Biostrings::writeXStringSet(fa_out)
message("Complete.")
