#!/usr/bin/env -S Rscript --no-save --no-restore
# modified 5/18 to add workaround for autorenaming of header row

# OPENBLAS_NUM_THREADS=12  /home/jason/miniconda3/envs/musta2/bin/Rscript --slave --vanilla /rgdata/home/jason/musta2/lib/recalc_pbcount.R /rgdata/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq-20220518.csv /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.correlation.txt /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.existence.txt 12 0 /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.salmon.tpm.txt /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.salmon.count.txt

#args <- "/rgdata/home/jason/musta_kk1402-20220406/run_kk1402_rnaseq-20220518.csv /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.correlation.txt /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.existence.txt 12 0 /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.salmon.tpm.txt /rgdata/home/jason/musta_kk1402-20220406/test_musta2-20220615_retest/result/sepFalse.salmon.count.txt"
#arg0 <- unlist(strsplit(args, " "))
arg0 <- commandArgs(trailing = TRUE)
f.input <- list(config = arg0[1], corr = arg0[2], exist = arg0[3], mcores = as.numeric(arg0[4]), dupl.rm = as.numeric(arg0[5]), salmonTPM = arg0[6], salmonCount = arg0[7])
f.output <- list(
	ccsList     = gsub("correlation.txt$", "ccsList.txt",        f.input$corr),
	pbCount     = gsub("correlation.txt$", "pbcounts.txt",       f.input$corr),
	pbReads     = gsub("correlation.txt$", "ccsReads.txt",       f.input$corr),
	IS3Comp     = gsub("correlation.txt$", "IS3Comp.txt",        f.input$corr),
	salmonCount = gsub("count.txt$",       "count-filtered.txt", f.input$salmonCount),
	salmonTPM   = gsub("tpm.txt$",         "tpm-filtered.txt",   f.input$salmonTPM))

if (is.na(f.input$mcores)) {
	#f.input$mcores <- max(1, as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS")))
	f.input$mcores <- 12
	message("# Invalid thread count option specified; setting to 1")
}
if (is.na(f.input$dupl.rm)) {
	f.input$dupl.rm <- NULL
} else if (tolower(f.input$dupl.rm) == "auto") {
	f.input$dupl.rm <- NA
	# not implemented as of Apr 2022
} else {
	f.input$dupl.rm <- min(f.input$dupl.rm, 1)
}

for (prereq in unlist(f.input)[1:3]) {
	if (!file.exists(prereq)) {
		stop(sprintf("File not found: %s", prereq), call. = FALSE)
	}
}

if (!file.exists(f.output$pbCount)) {
	suppressPackageStartupMessages(library(data.table))
	suppressPackageStartupMessages(library(parallel))
	data.table::setDTthreads(f.input$mcores)

	musta.corr  <- fread(f.input$corr, header = TRUE, sep = "\t")
	musta.exist <- fread(f.input$exist, header = TRUE, sep = "\t")
	csv.input   <- fread(f.input$config)
	csv.input$cluster_report <- gsub("^\\s", "", csv.input$cluster_report)
	bs <- unique(names(table(musta.corr$sample)))
	bs.is_id <- vector(mode = "list", length = length(bs))
	bs.cr    <- vector(mode = "list", length = length(bs))
	bs.pbt   <- vector(mode = "list", length = length(bs))
	names(bs.pbt) <- bs

	for (i in 1:length(bs)) {
		message(sprintf("Importing IsoSeq CCS read information for sample %s...", bs[i]))
		bs.pbt[[i]]   <- musta.exist$transcript_id[unlist(musta.exist[,eval(bs[i]), with = FALSE])]
		bs.is_id[[i]] <- musta.corr[sample == bs[i],]
		bs.is_id[[i]] <- unique(bs.is_id[[i]][transcript_id %in% bs.pbt[[i]],list(transcript_id, seq_id, g_id_sample, tx_id_sample)])
		bs.is_id[[i]]$key <- sprintf("%s:%s", bs[i], gsub("^.*_transcript\\/", "transcript/", bs.is_id[[i]]$seq_id))
		bs.cr[[i]] <- fread(csv.input$cluster_report[match(bs[i], csv.input$sample)])
		bs.cr[[i]]$key <- sprintf("%s:%s", bs[i], bs.cr[[i]]$cluster_id)
	}
	cr <- unique(data.table::rbindlist(bs.cr)) # remove duplicate entries as a possible consequence of technical replicates?
	ccs.samples <- matrix(0, nrow = dim(cr)[1], ncol = length(bs))
	for (i in 1:length(bs)) {
		ccs.samples[which(cr$read_id %in% bs.cr[[i]]$read_id), i] <- 1
	}
	rownames(ccs.samples) <- cr$read_id

	message("Correlating IsoSeq CCS read information...")
	musta.corr$key <- sprintf("%s:%s", musta.corr$sample, musta.corr$seq_id)
	m_match <- musta.corr[,list(Gene = paste(unique(gene_id), collapse = "/"), Transcript = paste(unique(transcript_id), collapse = "/")), by = list(key)]
	m_match <- m_match[match(cr$key, m_match$key),]
	crf <- unique(data.frame(cr$read_id, cr$cluster_id, m_match$Gene, m_match$Transcript, ccs.samples))
	colnames(crf) <- c("CCS", "IsoSeq3", "MuSTA_gene", "MuSTA_transcript", bs)

	if (!is.null(f.input$dupl.rm)) {
		# duplicate CCS reads will not be counted
		crf_dupl <- crf$CCS[grep("/", crf$MuSTA_transcript)]
	} else {
		crf_dupl <- NULL
	}

	# For auto score assignment
	if (is.na(f.input$dupl.rm)) {
		mts <- musta.corr[,list(NQT = length(unique(transcript_id))), c("key")]
		mts <- mts[match(cr$key, mts$key), list(NQT)]
		mts <- ifelse(is.na(mts$NQT), yes = 0, no = mts$NQT)
		names(mts) <- cr$read_id
		mts <- mts[mts > 1]
	} else {
		mts <- NULL
	}

	musta.tc <- matrix(0, nrow = dim(musta.exist)[1], ncol = length(bs))
	colnames(musta.tc) <- bs
	rownames(musta.tc) <- sprintf("%s:%s", musta.exist$gene_id, musta.exist$transcript_id)

	musta.tid <- matrix("", nrow = dim(musta.exist)[1], ncol = length(bs))
	colnames(musta.tid) <- bs
	rownames(musta.tid) <- sprintf("%s:%s", musta.exist$gene_id, musta.exist$transcript_id)

	get.ccs <- function(j, id, cluster.report, lookup, hash, dupl = NULL, score = 0, dupl.occurrences = NULL) {
		cur.reads <- cluster.report[key %in% unique(lookup[transcript_id %in% hash$transcript_id[j],]$key),]$read_id
		excl.reads <- NULL
		if (!is.null(dupl)) {
			excl.reads <- cur.reads[cur.reads %in% dupl]
			cur.reads  <- cur.reads[!(cur.reads %in% dupl)]
		}
		if (length(excl.reads) < 1) {
			excl.reads <- NULL
		}
	
		x <- list(Reads = cur.reads, Duplicates = excl.reads)
		if (is.null(x$Duplicates)) {
			x$Score <- length(x$Reads)
		} else if (is.na(score)) {
			if (is.null(dupl.occurrences)) {
				x$Score <- length(x$Reads)
			} else {
				# auto mode, not implemented in Apr 2022
				x$Score <- (length(x$Reads) * 1) + sum(1/(dupl.occurrences[names(dupl.occurrences) %in% x$Duplicates]))
			}		
		} else {
			x$Score <- (length(x$Reads) * 1) + (length(x$Duplicates) * score)
		}
		x$uniqueReads <- paste(unique(x$Reads), collapse = ";")
		return(x)
	}
	
	for (i in 1:length(bs)) {
		message(sprintf("Processing BioSample %s...", bs[i]))
		j <- which(unlist(musta.exist[,eval(bs[i]), with = FALSE]))
		cur.ccs <- parallel::mclapply(j, get.ccs, 
			cluster.report = bs.cr[[i]], 
			lookup = bs.is_id[[i]], 
			hash = musta.exist, 
			dupl = crf_dupl,
			score = f.input$dupl.rm,
			dupl.occurrences = mts,
			mc.cores = f.input$mcores)
			musta.tc[j,i]  <- as.vector(sapply(cur.ccs, function(x) x$Score))
			musta.tid[j,i] <- as.vector(sapply(cur.ccs, function(x) x$uniqueReads))
			invisible(rm(cur.ccs))
			invisible(gc())
		}
		invisible(gc())

	musta.tcf <- data.frame(musta.exist$gene_id, musta.exist$transcript_id, musta.tc)
	musta.idf <- data.frame(musta.exist$gene_id, musta.exist$transcript_id, musta.tid)
	musta.tcf <- musta.tcf[apply(musta.tc, 1, sum) > 0,] 
	musta.idf <- musta.idf[apply(musta.tc, 1, sum) > 0,]
	colnames(musta.tcf) <- c("gene_id", "transcript_id", bs)
	colnames(musta.idf) <- c("gene_id", "transcript_id", bs)

	crf$MuSTA_gene       <- ifelse(is.na(crf$MuSTA_gene),       yes = "", no = crf$MuSTA_gene)
	crf$MuSTA_transcript <- ifelse(is.na(crf$MuSTA_transcript), yes = "", no = crf$MuSTA_transcript)

	message("Filtering Salmon outputs...")
	bs_tempHeader <- tempfile()
	cat(sprintf("ID\t%s\nID\t%s\n", paste(bs, collapse = "\t"), paste(bs, collapse = "\t")), sep = "\n", file = bs_tempHeader)
	bs_header <- read.table(file = bs_tempHeader, sep = "\t", header = TRUE)
	
	colnames_count_rewrite <- TRUE
	colnames_tpm_rewrite <- TRUE
	if (!identical(colnames(bs_header), as.character(bs_header[1,]))) {
		message("Warning: IDs may be modified during data export due to incompatible headers.")
		colnames_count_rewrite <- as.character(bs_header[1,])
		colnames_tpm_rewrite <- as.character(bs_header[1,])
	}
	if (dim(bs_header)[2] == length(bs) + 1) {
		invisible(file.remove(bs_tempHeader))
	}
	if (file.exists(f.input$salmonCount) & file.exists(f.input$salmonTPM)) {
		salmon_count <- read.table(f.input$salmonCount, header = TRUE, sep = "\t")
		if (colnames(salmon_count)[1] == "ID") {
			salmon_count <- salmon_count[match(musta.tcf$transcript_id, salmon_count$ID),]
			salmon_count <- salmon_count[,match(colnames(bs_header), colnames(salmon_count))]
			colnames_count_rewrite <- as.character(bs_header[1,])
		} else {
			# remove from final version - keep just in case for now
			salmon_count <- salmon_count[match(musta.tcf$transcript_id, rownames(salmon_count)),]
			salmon_count <- salmon_count[,match(colnames(bs_header)[-1], colnames(salmon_count))]
			colnames_count_rewrite <- as.character(bs_header[1,-1])
		}
		#salmon_count <- as.data.frame(na.omit(salmon_count)) # added 6/7

		salmon_tpm <- read.table(f.input$salmonTPM, header = TRUE, sep = "\t")
		if (colnames(salmon_tpm)[1] == "ID") {
			salmon_tpm <- salmon_tpm[match(musta.tcf$transcript_id, salmon_tpm$ID),]
			salmon_tpm <- salmon_tpm[,match(colnames(bs_header), colnames(salmon_tpm))]
			colnames_tpm_rewrite <- as.character(bs_header[1,])
		} else {
			salmon_tpm <- salmon_tpm[match(musta.tcf$transcript_id, rownames(salmon_tpm)),]
			salmon_tpm <- salmon_tpm[,match(colnames(bs_header)[-1], colnames(salmon_tpm))]
			colnames_tpm_rewrite <- as.character(bs_header[1,-1])
		}
		#salmon_tpm <- as.data.frame(na.omit(salmon_tpm)) # added 6/7
		write.table(salmon_count, file = f.output$salmonCount, quote = FALSE, sep = "\t", row.names = TRUE, col.names = colnames_count_rewrite)
		write.table(salmon_tpm, file = f.output$salmonTPM, quote = FALSE, sep = "\t", row.names = TRUE, col.names = colnames_tpm_rewrite)
	}

	message("Writing MuSTA PBcounts output...")
	write.table(crf, file = f.output$ccsList, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	write.table(musta.tcf, file = f.output$pbCount, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	write.table(musta.idf, file = f.output$pbReads, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	message("PBcounts calculation complete.")
	q(status = 0)
} 
message(sprintf("File with the same name as '%s' already exists...\nExiting.", f.output$pbCount)) 
