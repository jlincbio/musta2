POSTINSTALLCHECK <- TRUE # to be modified by install.pl
if (!POSTINSTALLCHECK) {
	required_packages_cran <- c('argparse', 'BiocManager', 'caret', 'clisymbols', 'config', 'devtools', 'dplyr', 'drake', 'DT', 'e1071', 'forcats', 'fs', 'furrr', 'future', 'ggplot2', 'ggplotify', 'gridBase', 'gridExtra', 'htmltools', 'jsonlite', 'optparse', 'pkgload', 'plotly', 'plyr', 'pROC', 'purrr', 'randomForest', 'readr', 'remotes', 'reshape', 'rlang', 'rmarkdown', 'scales', 'stringi', 'stringr', 'tibble', 'tidyr')
	required_packages_bioc <- c("BSgenome", "Biostrings", "rtracklayer", "GenomicRanges", "GenomicFeatures", "NOISeq")
	
	message_bioc <- ""
	message_cran <- ""
	
	####CRAN####
	not_installed_packages_cran <- required_packages_cran[!(required_packages_cran %in% rownames(installed.packages()))]
	if (length(not_installed_packages_cran) > 0L) {
		not_installed_flag <- TRUE
		message_cran <- sprintf("install.packages(\'%s\', dependencies = TRUE)", not_installed_packages_cran)
	}
	
	####Bioconductor####
	not_installed_packages_bioc <- required_packages_bioc[!(required_packages_bioc %in% rownames(installed.packages()))]
	if (length(not_installed_packages_bioc) > 0L) {
		message_bioc <- sprintf("BiocManager::install(\'%s\')", not_installed_packages_bioc)
	}
	
	####Summary####
	if ((length(not_installed_packages_cran) > 0) & (length(not_installed_packages_bioc) > 0)) {
		stop(sprintf("\n### MISSING R PACKAGES DETECTED: INSTALL WITH THE FOLLOWING COMMANDS IN R:\n%s\n###\n",
			paste(c(message_cran, message_bioc), collapse = "\n")), call. = FALSE)
	}
}


