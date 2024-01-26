#!/usr/bin/env -S Rscript --no-save --no-restore

# MuSTA2 launcher, modified 2022/06/08
musta2_invoked_cmds <- paste("musta2", commandArgs(trailing = TRUE), collapse = " ")
Sys.setenv(SGE_ROOT = "/opt/sge", SMRT_ROOT = "/opt/pacbio/smrtlink", CUR_MUSTA2_INVOKED = musta2_invoked_cmds)

.name    <- "MuSTA"
.version <- "2.0.0a"

this_file <- function() {
	cmdArgs <- commandArgs(trailingOnly = FALSE)
	needle <- "--file="
	match <- grep(needle, cmdArgs)
	if (length(match) > 0) {
		# Rscript
		res <- sub(needle, "", cmdArgs[match])
	} else if (identical(cmdArgs, c("RStudio", "--interactive"))) {
		# Rstudio interactive mode
		res <- rstudioapi::getActiveDocumentContext()$path
	} else {
		# 'source'd via R console
		res <- sys.frames()[[1]]$ofile
	}
	if (res == "") {
		return(".")
	}
	return(res)
}
is_fun <- function(obj_name) is.function(eval(parse(text = obj_name)))
dir_log_path <- function(softname, o_dir = output_dir) path(o_dir, "log", softname)

## Dependency check and library initialization
.this_dir <- gsub("\\/lib$", "", normalizePath(dirname(this_file()))) # not detected by ls()
source(sprintf("%s/lib/check_dependencies.R", .this_dir))
suppressWarnings(suppressPackageStartupMessages({library(argparse, quietly = TRUE)}))
source(sprintf("%s/lib/parse_arg.R", .this_dir))
if (Sys.getenv("CONDA_DEFAULT_ENV") == "") {
	stop("MuSTA2 Miniconda environment is required to continue.", call. = FALSE)
}

if (file.exists(args$input)) {
	args$input <- normalizePath(args$input)
}

if (file.exists(args$gtf)) {
	args$gtf <- normalizePath(args$gtf)
}

suppressPackageStartupMessages(library(rlang, quietly = TRUE))
if (!identical(base::grep("/", args$prefix), integer(0))) abort("The value of '--prefix' must NOT contain '/'", "requirement error")

suppressPackageStartupMessages({
	library(stringr, quietly = TRUE)
	library(readr, quietly = TRUE)
	library(dplyr, quietly = TRUE)
	library(tidyr, quietly = TRUE)
	library(purrr, quietly = TRUE)
	library(fs, quietly = TRUE)
	library(ggplot2, quietly = TRUE)
	library(drake, quietly = TRUE)
})

obj <- ls()
obj <- setdiff(obj[!map_lgl(obj, is_fun)], "args")
rm(list = c(obj, "obj"))

## Importing other components
source(paste0(.this_dir, "/lib/file_location.R"))  # written as relative paths from THIS file
pipeline_dir <- path_abs(.this_dir) # get the directory path of THIS file
lib <- map(lib, ~ path(pipeline_dir, .)) # as absolute paths
src <- map(src, ~ path(pipeline_dir, .)) # as absolute paths

source(src$utils.R)
load_jobwatcher_chr <- sprintf("suppressPackageStartupMessages(pkgload::load_all('%s', export_all = FALSE, quiet = TRUE))", lib$jobwatcher)
eval(parse(text = load_jobwatcher_chr))
jobwatcher_mode <- get_jobwatcher_mode()

## Display parameters and start MuSTA
cat(sprintf("Parameters for %s, version %s:\n ", .name, .version))
cat(sprintf(" %13s: %s\n", names(args), args))
cat("\n====== MUSTA PROCESS BEGINS ======\n")
if (args$dry_run) {
	print(if (requireNamespace("sessioninfo", quietly = TRUE)) sessioninfo::session_info() else sessionInfo())
	cat("\n====== MUSTA DRY-RUN BEGINS ======\n")
}

source(src$parse_config.R) # use absolute paths instead; first implemented 20211029
for (i in 1:length(soft)) {
	soft[[i]][1] <- Sys.which(soft[[i]][1])
}

output_subdir <- c("log", "script", "report", "genome", "samples", "merge", "result", "fusion", "merge/salmon", "merge/sqanti_filter_supplement") %>% path(output_dir, .)
dir <- map(output_subdir, force) %>% setNames(str_remove(output_subdir, "^.*/"))
sample_dir_root <- sprintf("%s/samples", output_dir)
sample_dirs <- path(dir$samples, samples)
log_categories <- c("file_management", "interleave", "lordec", "minimap2", "samtools", "merge", "salmon", "sqanti", "fusion", "others")
dir_logs <- map(log_categories, dir_log_path) %>% setNames(log_categories)

if (!args$dry_run) {
	# originally permissions set to 775, changed to 777 to avoid unexplained write errors with Salmon 1.5.2
	# solution 1 (testing): force newer version of Salmon in environment as of 6/13/2022
	walk(dir, dir_create, mode = "0777")
	walk(sample_dirs, dir_create, mode = "0777")
	walk(dir_logs, dir_create, mode = "0777")
	inform(sprintf("Output folder created: %s", output_dir))
}

qsub_proxy_param <- ''
cur_proxy_setting <- Sys.getenv('HTTP_PROXY')
if (args$proxy != '') {
	if (cur_proxy_setting == '') {
		Sys.setenv("HTTP_PROXY"=args$proxy)
		qsub_proxy_param <- "#$ -v HTTP_PROXY"
		message("Temporarily setting system proxy settings...")
	}
}

####intermediate file location####
source(src$intermediate_files.R)

####tidy IO files####
source(src$update_input_files.R)
IO_summary_file <- path(dir$report, "IO.summary.csv")
if (!args$dry_run) readr::write_csv(input_files, IO_summary_file)
####Qscripts####
source(src$make_qscript.R)

####plan####
source(src$plan.R)
eval(parse(text = drake_txt))
if (args$force && !args$dry_run) clean(garbage_collection = TRUE, jobs = 1, purge = TRUE)
wkflow_plot <- drake_ggraph(drake_config(df_plan), targets_only = TRUE, mode = "out", label_nodes = TRUE) + theme_void() + coord_flip() + ggtitle(sprintf("%s (v %s) workflow", .name, .version))
suppressMessages({wkflow_plot <- wkflow_plot + scale_x_continuous(trans = "reverse")})
wkflow_plot$layers <- map(wkflow_plot$layers, ~ {if (is(.x[["geom"]], "GeomText")) .x[["aes_params"]][["size"]] <- 3L;.x})
ggsave(if (!args$dry_run) path(dir$report, "plan_pre_run.pdf") else "plan_pre_run.pdf", wkflow_plot, width = 32, height = 18, units = "cm")

####check file existence####
## below added 01/12/22 to allow possible short-read skipping
qrecall_files <- c(input_files$long_read_hq, ref_gtf, genome_fa)
if (!args$no_salmon) {
	qrecall_files <- c(qrecall_files, input_files$short_read_1, input_files$short_read_2)
}
## above added 01/12/22 to allow possible short-read skipping
if (args$use_lq) qrecall_files <- c(qrecall_files, input_files$long_read_lq)
not_exist_files <- qrecall_files[!fs::file_exists(qrecall_files)]
if (length(not_exist_files) > 0L) abort(paste0("These files are not found: ", str_c(not_exist_files, collapse = ", ")), "requirement error")

## Dry run
if (args$dry_run) {
	cat(print_list(file))
	cat("\n======  MUSTA DRY-RUN ENDS  ======\n")
	cat(drake_txt)
	cat("\n")
	quit(save = "no", status = 1, runLast = FALSE)
}

## Actual run & drake report generation
future::plan(future::multicore, workers = 1)
make(df_plan, jobs = 1, parallelism = "future", prework = load_jobwatcher_chr)
wkflow_plot <- drake_ggraph(drake_config(df_plan), targets_only = TRUE, mode = "out", label_nodes = TRUE) + theme_void() + coord_flip() + ggtitle(sprintf("%s (v %s) workflow", .name, .version))
suppressMessages({wkflow_plot <- wkflow_plot + scale_x_continuous(trans = "reverse")})
wkflow_plot$layers <- map(wkflow_plot$layers, ~ {if (is(.x[["geom"]], "GeomText")) .x[["aes_params"]][["size"]] <- 3L;.x})
ggsave(path(dir$report, "plan_post_run.pdf"), wkflow_plot, width = 32, height = 18, units = "cm")
#Sys.setenv("HTTP_PROXY"=cur_proxy_setting)
cat("\n======  MUSTA RUN FINISHED  ======\n")
quit(save = "no", status = 0)
