Sys.setenv("http_proxy"="http://proxy.chiba-cc.jp:8080")

args <- commandArgs(TRUE)
pipeline_name <- args[1] # pipe
srcdir        <- args[2] # abs_path to fafile_dir
dcf_path      <- args[3] # path to write DCF (seed) file
pkgdir        <- args[4] # abs_path to write package files

date <- Sys.Date()
user <- "user"
fa_name <- fs::path_file(srcdir) #example: hg38
suffix <- ".fa"
chroms <- fs::path_ext_remove(fs::path_file(fs::dir_ls(srcdir, type = "file")))
chroms <- gtools::mixedsort(chroms)
chroms_Rexp <- sprintf("c('%s')", paste(chroms, collapse = "', '"))

pkgname <- sprintf("BSgenome.%s.%s.%s", fa_name, user, pipeline_name)
pkgpath <- fs::path(pkgdir, pkgname)

# Rsubread adjustments: last updated 6/1/2023
# pkgpath -> dir$genome
genome_fa <- args[5]
memLimit <- 1024 * min(as.numeric(args[6]), 1024)
index_base <- paste(rev(rev(unlist(strsplit(genome_fa, "\\.")))[-1]), collapse = ".")
index_subread <- sprintf("%s/rsubread_%s", pkgdir, basename(index_base))
index_subread_fileList <- list.files(path = dirname(index_subread), pattern = basename(index_subread), full.names = TRUE)
index_exists <- (length(file.exists(index_subread_fileList)) > 0) & as.logical(prod(file.exists(index_subread_fileList)))
index_sizes  <- length(which(file.info(index_subread_fileList)$size >= 1048576))
if (as.logical(args[7])) {
	if (!(index_exists & (index_sizes > 1))) {
		message("## RNA-seq genome index not found, creating...")
		Rsubread::buildindex(basename = index_subread, reference = genome_fa, memory = memLimit)
	}
} else {
	message("## Skipping RNA-seq index construction")
}

# release_name has apparently been deprecated
# release_name: {pipeline_name}
# provider_version: {pipeline_name}
# circ_seqs: needed for non UCSC or NCBI; for now assume character(0)

dcf <- c(
	sprintf("Package: %s", pkgname),
	sprintf("Title: Full genome sequences for %s", fa_name),
	sprintf("Description: Full genome sequences for %s as provided by %s and stored in Biostrings objects.", fa_name, user),
	"Version: 0.0.1",
	sprintf("Author: %s", user),
	sprintf("Maintainer: %s <your_email_address>", user),
	"License: Depend on your sequence data file",
	sprintf("organism: %s", fa_name),
	sprintf("organism_biocview: %s", fa_name),
	"common_name: unknown",
	sprintf("genome: %s", pipeline_name),
	sprintf("provider: %s", user),
	sprintf("release_date: %s", date),
	sprintf("BSgenomeObjname: %s", fa_name),
	sprintf("seqs_srcdir: %s", srcdir),
	sprintf("seqnames: %s", chroms_Rexp),
	"circ_seqs: character(0)",
	sprintf("seqfiles_suffix: %s", suffix))
	cat(dcf, file = dcf_path, sep = "\n")

if (fs::dir_exists(pkgpath)) fs::dir_delete(pkgpath)
BSgenome::forgeBSgenomeDataPkg(dcf_path, destdir = pkgdir)


