#!/usr/bin/env -S Rscript --no-save --no-restore
args <- commandArgs(trailing = TRUE)

gtf_input  <- args[1] # make_qscripts.R: MuSTA output GTF: file$sqanti_gtf
gtf_anno   <- args[2] # parse_args.R: Gencode GTF: args$gtf (inherited from launcher)
mc_cores   <- floor(min(as.numeric(args[3]), as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))/2, na.rm = TRUE))

if (is.na(mc_cores)) {
	message("# Number of threads not specified, reverting to single-threaded mode.")
	mc_cores <- 1
} else {
	message(sprintf("# Multithreading mode: %d threads", mc_cores))
}

if ((!file.exists(gtf_input)) & (!file.exists(gtf_anno))) {
	stop("Both query and database GTFs need to be specified.", call. = FALSE)
	q(55)
}

gtf_output <- gsub("\\.gtf$", "-anno.gtf", gtf_input)

Sys.setenv("OPENBLAS_NUM_THREADS" = mc_cores)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
data.table::setDTthreads(mc_cores)

message("# Generating annotations...")
anno.tmpfile <- tempfile()
cmd <- sprintf("%s intersect -a %s -b %s -wao > %s", Sys.which("bedtools"), gtf_input, gtf_anno, anno.tmpfile)
s.cmd <- system(cmd)

if (s.cmd > 0) {
	stop("Annotation error.", call. = FALSE)
	q(s.cmd)
}

message("# Processing MuSTA output...")
musta.gtf <- data.table::fread(anno.tmpfile, header = FALSE, sep = "\t")
musta.gtf$gene_id <- 
	as.character(mclapply(musta.gtf$V9, function(x) trimws(gsub("\"", "", gsub("gene_id ", "", grep(pattern = "gene_id", x = unlist(strsplit(x, ";")), value = TRUE)))), mc.cores = 16))
musta.gtf$transcript_id <- 
	as.character(mclapply(musta.gtf$V9, function(x) trimws(gsub("\"", "", gsub("transcript_id ", "", grep(pattern = "transcript_id", x = unlist(strsplit(x, ";")), value = TRUE)))), mc.cores = 16))
musta.gtf$gene_name <- 
	as.character(mclapply(musta.gtf$V18, function(x) trimws(gsub("\"", "", gsub("gene_name ", "", grep(pattern = "gene_name", x = unlist(strsplit(x, ";")), value = TRUE)))), mc.cores = 16))
musta.gtf$transcript_name <- 
	as.character(mclapply(musta.gtf$V18, function(x) trimws(gsub("\"", "", gsub("transcript_name ", "", grep(pattern = "transcript_name", x = unlist(strsplit(x, ";")), value = TRUE)))), mc.cores = 16))
musta.gtf$Start <- musta.gtf$V4
musta.gtf$End <- musta.gtf$V5

musta.gtf.transcript <-
	musta.gtf[,list(V2 = ".", V6 = ".", V8 = ".", V3 = "transcript", Start = min(V4), End = max(V5), gene_symbol = paste(unique(gene_name), collapse = "/")), by = c("V1", "V7", "transcript_id", "gene_id")]
musta.gtf.genes <-
	musta.gtf[,list(V2 = ".", V6 = ".", V8 = ".", V3 = "gene", Start = min(V4), End = max(V5), gene_symbol = paste(unique(gene_name), collapse = "/")), by = c("V1", "V7", "gene_id")]

musta.gtf.transcript$V9 <- 
	gsub("character\\(0\\)", "", sprintf("gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\";", 
		as.character(musta.gtf.transcript$gene_id), 
		as.character(musta.gtf.transcript$transcript_id), 
		as.character(musta.gtf.transcript$gene_symbol), 
		as.character(musta.gtf.transcript$transcript_id)))
musta.gtf.genes$V9 <- 
	gsub("character\\(0\\)", "", sprintf("gene_id \"%s\"; gene_name \"%s\";", 
		as.character(musta.gtf.genes$gene_id),
		as.character(musta.gtf.genes$gene_symbol)))
musta.gtf$V9 <- 
	gsub("character\\(0\\)", "", sprintf("gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\";",
		as.character(musta.gtf$gene_id), 
		as.character(musta.gtf$transcript_id), 
		as.character(musta.gtf$gene_name), 
		as.character(musta.gtf$transcript_id)))
	
musta.gtf <- musta.gtf[,c("V1", "V2", "V3", "Start", "End", "V6", "V7", "V8", "V9")]
musta.gtf.genes <- musta.gtf.genes[,c("V1", "V2", "V3", "Start", "End", "V6", "V7", "V8", "V9")]
musta.gtf.transcript <- musta.gtf.transcript[,c("V1", "V2", "V3", "Start", "End", "V6", "V7", "V8", "V9")]
musta.anno <- rbind(musta.gtf.genes, musta.gtf.transcript, musta.gtf)
musta.anno <- unique(musta.anno)
message(sprintf("# Writing annotated MuSTA output to %s...", gtf_output))
write.table(musta.anno, file = gtf_output, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
if (file.exists(anno.tmpfile)) {
	invisible(file.remove(anno.tmpfile))
}
