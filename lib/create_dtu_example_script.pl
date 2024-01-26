#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename 'dirname';

unless (@ARGV) {
	die "Usage: dtu_example.pl <MuSTA GTF> <MuSTA prefix> [DTU script output]\n";
}

my ($gtf_anno, $prefix, $f_out) = @ARGV;
my @scripts;

(my $d_musta = $gtf_anno) =~ s/\/result.*gtf$//;
(my $table_out = $gtf_anno) =~ s/\.gtf$/_dtuTable.txt/;
(my $html_out = $gtf_anno) =~ s/\.gtf$/_dtuResult.html/;

my %subst = (
	'___ANNOGTF___' => $gtf_anno, # -anno.gtf instead of .gtf in file$sqanti_gtf
	'___MUSTAOUTPUTDIR___' => $d_musta, # output_dir from launcher
	'___MUSTAPREFIX___' => $prefix, # args$prefix (from musta)
	'___DTUOUTPUTDIR___' => dirname($gtf_anno), # same directory as musta prefix
	'___DTUTABLEOUTPUT___' => $table_out,
	'___DTUHTMLOUTPUT___' => $html_out
);

while (<DATA>) {
	chomp;
	my $s = $_;
	foreach my $i (keys %subst) {
		$s =~ s/$i/$subst{$i}/g;
	}
	push @scripts, $s;
}

if (defined $f_out) {
	open my $F_OUTPUT, '>', $f_out or die "Can't create script!\n";
	print $F_OUTPUT "$_\n" for @scripts;
	close $F_OUTPUT;
} else {
	print "$_\n" for @scripts;
}

__DATA__
#!/usr/bin/Rscript --vanilla

# Example DTU analysis script for MuSTA results
# Modify the R script and execute in RStudio as necessary

# Reference: Tekath & Dugas, Bioinformatics 37(21): 3781-3787, 2021
# Example largely follows sample code on https://tobitekath.github.io/DTUrtle/

library(DTUrtle)
library(BiocParallel)
biocpar <- BiocParallel::MulticoreParam(4) # 4 threads multicore processing

ref_genome <- "hg38" # modify as necessary
f_input_gtf <- "___ANNOGTF___"

# Import annotation from GTF & map 1-to-1 for transcript and gene symbols
tx2gene <- import_gtf(gtf_file = f_input_gtf) 
tx2gene <- move_columns_to_front(
	df = tx2gene, columns = c("transcript_name", "gene_name"))
tx2gene$gene_name <- one_to_one_mapping(
	name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(
	name = tx2gene$transcript_name, id = tx2gene$transcript_id)

# Import Salmon transcript-level quantification data and scale
files <- Sys.glob("___MUSTAOUTPUTDIR___/samples/*/___MUSTAPREFIX___.salmon/quant.sf")
names(files) <- gsub(".*/", "", gsub("/___MUSTAPREFIX___.salmon/quant.sf", "", files))
cts <- import_counts(files = files, type = "salmon")

# Create a sample data to specify which sample/cell belongs to which group
# Expand this section as necessary
pd <- data.frame(
	id = colnames(cts),       # Sample IDs from MuSTA input CSV
	group = ...,              # e.g. group = c("Normal", "Tumor", "Normal", ...)
	stringsAsFactors = FALSE)

dturtle  <- run_drimseq(counts = cts, tx2gene = tx2gene, pd = pd, id_col = "id", 
	cond_col = "group", filtering_strategy = "bulk",
	BPPARAM = biocpar, cond_levels = pd$group)

dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)
dturtle <- create_dtu_table(dturtle = dturtle) # DTU output table
setwd("___DTUOUTPUTDIR___") # move to output folder to organize output data

# Create ouputs, save to ~savepath~ and link to `dtu_table`.
dturtle <- plot_proportion_barplot(
	dturtle = dturtle, savepath = "___DTUOUTPUTDIR___/images", 
	add_to_table = "barplot", BPPARAM = biocpar)

dturtle <- plot_proportion_pheatmap(
	dturtle = dturtle, savepath = "___DTUOUTPUTDIR___/images",
	include_expression = TRUE,
	add_to_table = "pheatmap", BPPARAM = biocpar)

dturtle <- plot_transcripts_view(
	dturtle = dturtle, gtf = f_input_gtf, genome = ref_genome, one_to_one = TRUE,
	savepath = "___DTUOUTPUTDIR___/images", 
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

plot_dtu_table(dturtle = dturtle, savepath = "___DTUHTMLOUTPUT___", column_formatters = column_formatter_list)

# Writing DTU analysis output table	
write.table(apply(dturtle$dtu_table, 2, as.character), 
	file = "___DTUTABLEOUTPUT___", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Optional: DGE analysis via DESeq2
# Import counts and summarize to gene-level
cts_dge <- import_dge_counts(files, type = "salmon", tx2gene = tx2gene)

# Perform the DGE calling with DESeq2
# Appropriate parameters are chosen based on bulk dge_calling_strategy
dturtle$dge_analysis <- run_deseq2(counts = cts_dge, pd = pd, id_col = "id",
	cond_col = "group", lfc_threshold = 0.5, sig_threshold = 0.01, 
	dge_calling_strategy = "bulk", BPPARAM = biocpar)
