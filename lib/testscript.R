#!/usr/bin/Rscript --vanilla

# Example DTU analysis script for MuSTA results
# Modify the R script and execute in RStudio as necessary

# Reference: Tekath & Dugas, Bioinformatics 37(21): 3781-3787, 2021
# Example largely follows sample code on https://tobitekath.github.io/DTUrtle/

library(DTUrtle)
library(BiocParallel)
biocpar <- BiocParallel::MulticoreParam(4) # 4 threads multicore processing

f_input_gtf <- "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/musta-anno.gtf"

# Import annotation from GTF & map 1-to-1 for transcript and gene symbols
tx2gene <- import_gtf(gtf_file = f_input_gtf) 
tx2gene <- move_columns_to_front(
	df = tx2gene, columns = c("transcript_name", "gene_name"))
tx2gene$gene_name <- one_to_one_mapping(
	name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(
	name = tx2gene$transcript_name, id = tx2gene$transcript_id)

# Import Salmon transcript-level quantification data and scale
files <- Sys.glob("/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/samples/*/musta.salmon/quant.sf")
names(files) <- gsub(".*/", "", gsub("/musta.salmon/quant.sf", "", files))
cts <- import_counts(files = files, type = "salmon")

# Create a sample data to specify which sample/cell belongs to which group
# Expand this section as necessary
samps_key <- list("bc1005" = "NC_SHC", "bc1006" = "NC_SH3", "bc1008" = "IFNg_SHC", "bc1012" = "IFNg_SH3")
samps <- data.frame(sample_id = colnames(cts))
samps$group <- as.factor(as.character(samps_key[match(gsub("[abc]$", "", samps$sample_id), names(samps_key))]))
samps$replicates <- sprintf("%s_%s", samps$group, gsub("^.*([abc])$", "\\1", samps$sample_id))

pd <- data.frame(id = colnames(cts), group = samps$group, stringsAsFactors = FALSE)

#pd <- data.frame(
#	id = colnames(cts),       # Sample IDs from MuSTA input CSV
#	group = ...,              # e.g. group = c("Normal", "Tumor", "Normal", ...)
#	stringsAsFactors = FALSE)

dturtle  <- run_drimseq(counts = cts, tx2gene = tx2gene, pd = pd, id_col = "id", 
	cond_col = "group", filtering_strategy = "bulk",
	BPPARAM = biocpar, cond_levels = c("NC_SHC", "NC_SH3"))
#	BPPARAM = biocpar, cond_levels = pd$group)
	
dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)
dturtle <- create_dtu_table(dturtle = dturtle) # DTU output table
setwd("/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result") # move to output folder to organize output data

# Create ouputs, save to ~savepath~ and link to `dtu_table`.
dturtle <- plot_proportion_barplot(
	dturtle = dturtle, savepath = "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/images", 
	add_to_table = "barplot", BPPARAM = biocpar)

dturtle <- plot_proportion_pheatmap(
	dturtle = dturtle, savepath = "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/images",
	include_expression = TRUE,
	add_to_table = "pheatmap", BPPARAM = biocpar)

dturtle <- plot_transcripts_view(
	dturtle = dturtle, gtf = gtf, genome = genome, one_to_one = TRUE,
	savepath = "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/images", 
	add_to_table = "transcript_view", BPPARAM = biocpar)
	
dturtle$dtu_table$musta_gid <- 
	tx2gene$gene_id[match(dturtle$dtu_table$gene_ID, tx2gene$gene_name)]
	
# Create interactive HTML-table from results data frame
# Optional: specify colorful column formatters
column_formatter_list <- list(
	"gene_qvalue" = table_pval_tile("white", "orange", digits = 3),
	"minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 3),
	"number_tx" = formattable::color_tile('white', "lightblue"),
	"number_significant_tx" = formattable::color_tile('white', "lightblue"),
	"max(Condition1-Condition2)" = table_percentage_bar('lightgreen', "#FF9999", digits = 2))

plot_dtu_table(dturtle = dturtle, savepath = "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/musta-anno_dtuResult.html", column_formatters = column_formatter_list)

# Writing DTU analysis output table	
write.table(apply(dturtle$dtu_table, 2, as.character), 
	file = "/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq_reorder-20220518/result/musta-anno_dtuTable.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

## Optional: DGE analysis via DESeq2
# Import counts and summarize to gene-level
cts_dge <- import_dge_counts(files, type = "salmon", tx2gene = tx2gene)

# Perform the DGE calling with DESeq2
# Appropriate parameters are chosen based on bulk dge_calling_strategy
dturtle$dge_analysis <- run_deseq2(counts = cts_dge, pd = pd, id_col = "id",
	cond_col = "group", lfc_threshold = 0.5, sig_threshold = 0.01, 
	dge_calling_strategy = "bulk", BPPARAM = biocpar)
