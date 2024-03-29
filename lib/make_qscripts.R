# Last modified 2022/06/17: SQANTI3 compatibility & MuSTA2 call invocation
os_req <- sprintf("## MuSTA call: %s", Sys.getenv("CUR_MUSTA2_INVOKED"))
n_sample <- nrow(input_files)
ary <- arrayjob_option(n_sample)
SGE_TASK_ID <- if (jobwatcher_mode == "hgc") "$SGE_TASK_ID" else seq_along(samples)

if (Sys.getenv("R_LIBS") != "") {
	r_lib_str <- "R_LIBS=''"
} else {
	r_lib_str <- ""
}
Rcmd <- sprintf("%s %s --slave --vanilla", r_lib_str, Sys.which("Rscript"))

qsub_header_param <- paste(
	c(qsub_proxy_param, "#$ -l h=mk01", "#$ -v PYENV_ROOT", "#$ -v PATH", "#$ -v PYTHONPATH",
	"#$ -sync y", sprintf("# TBI: %s", os_req), "set -eu\n", as_bash_array(input_files)), collapse = "\n")

qsub_template <- function(..., script_path, script_dir = dir$script, 
	name = NA_character_, first_line = binbash(), parallel = parallel_option(),
	arrayjob = arrayjob_option(), directory = directory_option(),
	use_bash_profile = FALSE, other_req = qsub_header_param, recursive = TRUE,
	add_time = TRUE, qsub_args = "", additional_args = list(
		watch = FALSE, verbose = args$verbose, modify_req = FALSE)) {
	
	force(list(name, first_line, parallel, arrayjob, directory, use_bash_profile, other_req, script_path, script_dir, recursive, add_time, qsub_args))
	function(moot) {
		qsubfile <- make_qsubfile(..., name = name, first_line = first_line, parallel = parallel, arrayjob = arrayjob, directory = directory, use_bash_profile = use_bash_profile, other_req = other_req)
		path <- write_qsubfile(x = qsubfile, path = fs::path(script_dir, script_path), recursive = recursive, add_time = add_time)
		res <- do.call(qsub, c(list(path = path, args = qsub_args), additional_args))
		if (all(res$exit_code != 0L)) abort(sprintf("The job %s at %s ended with potential errors.", res$ID, res$path), "Subprocess error")
		return(res)
	}
}

# added 4/14/2023 to extract memory and thread settings
get_req <- function(x, r, que_req_obj = que_req) {
	return(que_req_obj[[x]][[r]])
}

parallel_option_req <- function(x, que_req_obj = que_req) {
        parallel_option(slot = que_req_obj[[x]][["slot"]], memory = que_req_obj[[x]][["memory"]], ljob = que_req_obj[[x]][["ljob"]])
}

parallel_option_req <- function(x, que_req_obj = que_req) {
	parallel_option(slot = que_req_obj[[x]][["slot"]], memory = que_req_obj[[x]][["memory"]], ljob = que_req_obj[[x]][["ljob"]])
}

qslot <- function(x, que_req_obj = que_req) que_req_obj[[x]][["slot"]]

#debug
pl_debug <- qsub_template("echo Hello World!", script_path = "debug")
pl_err   <- qsub_template("CommandNotExist", script_path = "err")

####genome.fa2BSgenome####
# genome_fa added for forgeBSgenome on 6/1/2023
pl_BSgenome <- qsub_template(
    script_path = "fasta2BSgenome",
    parallel = parallel_option_req("fasta2BSgenome"),
    directory = directory_option(out = dir_logs[["others"]]),
    str_c("bash", lib$multifa2singlefa.sh, genome_fa, file$genome_dir, sep = " "),
    str_c(Rcmd, lib$forgeBSgenome.R, .name, file$genome_dir, file$genome_seed, dir$genome, genome_fa, get_req(x = "fasta2BSgenome", r = "memory"), as.logical(!args$no_salmon), sep = " "),
    str_c("rm -f", file$genome_seed, sep = " "),
    str_c("rm -rf", file$genome_dir, sep = " ")
)
# Subread index: 6/1/2023
subread_idx <- sprintf("%s/rsubread_%s", dir$genome, paste(rev(rev(unlist(strsplit(genome_fa, "\\.")))[-1]), collapse = "."))

####Copy####
pl_cp_1 <-  #copy
  qsub_template(
    script_path = "cp_1",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("cp", input_files$long_read_hq, input_files$long_copy_hq, sep = " ", collapse = "\n"),
    str_c("cp", input_files$short_read_1, input_files$short_copy_1, sep = " ", collapse = "\n")
  )

pl_cp_2 <-  #copy
  qsub_template(
    script_path = "cp_2",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("cp", input_files$short_read_2, input_files$short_copy_2, sep = " ", collapse = "\n")
  )

pl_cp_lq <-  #copy
  qsub_template(
    script_path = "cp_2",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("cp", input_files$long_read_lq, input_files$long_copy_lq, sep = " ", collapse = "\n"),
  )

pl_fq2fa_hq <- 
  qsub_template(
    arrayjob = ary,
    script_path = "fq2fa_hq",
    directory = directory_option(out = dir_logs[["file_management"]]),
    paste0("#", soft$seqkit, " fq2fa ${long_copy_hq[", SGE_TASK_ID, "]} > ${lordec_hq[", SGE_TASK_ID, "]}")
  )

pl_fq2fa_lq <- 
  qsub_template(
    arrayjob = ary,
    script_path = "fq2fa_lq",
    directory = directory_option(out = dir_logs[["file_management"]]),
    paste0("#", soft$seqkit, " fq2fa ${long_copy_lq[", SGE_TASK_ID, "]} > ${lordec_lq[", SGE_TASK_ID, "]}")
  )

####Minimap2####
pl_minimap2_make_mmi <- 
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    script_path = "minimap2_make_mmi",
    directory = directory_option(out = dir_logs[["minimap2"]]),
    str_c(soft$minimap2, "-x splice:hq -uf ",
          "-d", file$genome_mmi, genome_fa,
          sep = " ")
  )

pl_minimap2_hq_sam <- 
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    arrayjob = ary,
    script_path = "minimap2_hq_sam",
    directory = directory_option(out = dir_logs[["minimap2"]]),
    str_c(soft$minimap2, " -ax splice:hq -uf ",
          file$genome_mmi, 
          " ${lordec_hq[", SGE_TASK_ID, "]} > ${minimap2_hq_sam[", SGE_TASK_ID, "]}"
          )
  )

pl_minimap2_hq_paf <- 
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    arrayjob = ary,
    script_path = "minimap2_hq_paf",
    directory = directory_option(out = dir_logs[["minimap2"]]),
    str_c(soft$minimap2, " -cx splice:hq -uf --cs ",
          file$genome_mmi, 
          " ${lordec_hq[", SGE_TASK_ID, "]} > ${minimap2_hq_paf[", SGE_TASK_ID, "]}"
          )
  )

pl_minimap2_lq_sam <- 
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    arrayjob = ary,
    script_path = "minimap2_lq_sam",
    directory = directory_option(out = dir_logs[["minimap2"]]),
    str_c(soft$minimap2, " -ax splice:hq -uf ",
          file$genome_mmi, 
          " ${lordec_lq[", SGE_TASK_ID, "]} > ${minimap2_lq_sam[", SGE_TASK_ID, "]}"
          )
  )

pl_minimap2_lq_paf <- 
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    arrayjob = ary,
    script_path = "minimap2_lq_paf",
    directory = directory_option(out = dir_logs[["minimap2"]]),
    str_c(soft$minimap2, " -cx splice:hq -uf --cs ",
          file$genome_mmi, 
          " ${lordec_lq[", SGE_TASK_ID, "]} > ${minimap2_lq_paf[", SGE_TASK_ID, "]}"
          )
  )

pl_qual_filter <-
  qsub_template(
    parallel = parallel_option_req("qual_filter"),
    script_path = "qual_filter",
    directory = directory_option(out = dir_logs[["others"]]),
    str_c(Sys.which("perl"), lib$quality_filter_paf_sam,
          IO_summary_file,
          args$map_qual,
          sep = " ")
  )

#Samtools
pl_samtools_hq <- 
  qsub_template(
    parallel = parallel_option_req("samtools"),
    arrayjob = ary,
    script_path = "samtools_hq",
    directory = directory_option(out = dir_logs[["samtools"]]),
    str_c(soft$samtools,
          " sort -@ ", qslot("samtools"), " -O bam -o ${samtools_hq_bam[", SGE_TASK_ID, "]}",
          " ${qfilt_hq_sam[", SGE_TASK_ID, "]}"
          ),
    str_c(soft$samtools, " index ${samtools_hq_bam[", SGE_TASK_ID, "]}")
  )

pl_samtools_lq <- 
  qsub_template(
    parallel = parallel_option_req("samtools"),
    arrayjob = ary,
    script_path = "samtools_lq",
    directory = directory_option(out = dir_logs[["samtools"]]),
    str_c(soft$samtools,
          " sort -@ ", qslot("samtools"), " -O bam -o ${samtools_lq_bam[", SGE_TASK_ID, "]}",
          " ${qfilt_lq_sam[", SGE_TASK_ID, "]}"
          ),
    str_c(soft$samtools, " index ${samtools_lq_bam[", SGE_TASK_ID, "]}")
  )

pl_samtools_merge <- 
  qsub_template(
    parallel = parallel_option_req("samtools"),
    arrayjob = ary,
    script_path = "samtools_merge",
    directory = directory_option(out = dir_logs[["samtools"]]),
    paste0(soft$samtools, " sort -@ ", qslot("samtools"), " -O bam -o ${samtools_bam[", SGE_TASK_ID, "]} ${merged_sam[", SGE_TASK_ID, "]}"),
    paste0(soft$samtools, " index ${samtools_bam[", SGE_TASK_ID, "]}")
  )

####merge hq and lq
pl_merge_fa <- 
  qsub_template(
    arrayjob = ary,
    script_path = "merge_fa",
    directory = directory_option(out = dir_logs[["lordec"]]),
    paste0("cat ${lordec_hq[", SGE_TASK_ID, "]} ${lordec_lq[", SGE_TASK_ID, "]} > ${merged_fa[", SGE_TASK_ID, "]}")
  )

pl_merge_sam <- 
  qsub_template(
    parallel = parallel_option_req("samtools"),
    arrayjob = ary,
    script_path = "merge_sam",
    directory = directory_option(out = dir_logs[["samtools"]]),
    str_c(soft$samtools, " merge -@ ", qslot("samtools"), " -f -O SAM",
          " ${merged_sam[", SGE_TASK_ID, "]}",
          " ${samtools_hq_bam[", SGE_TASK_ID, "]}",
          " ${samtools_lq_bam[", SGE_TASK_ID, "]}"
          )
  )

####Intra/Inter Merge####
pl_intra_merge <- 
  qsub_template(
    parallel = parallel_option_req("intra_merge"),
    arrayjob = ary,
    script_path = "intra_merge",
    directory = directory_option(out = dir_logs[["merge"]]),
    str_c(
      Rcmd, " ", lib$intra_merge.R, " ",  
      lib$isocan, " ",
      lib$merge_iso_commonpart.R,
      " ${qfilt_hq_paf[", SGE_TASK_ID, "]}",#inputs
      (if (is.null(input_files[["qfilt_lq_paf"]]) || any(is.na(input_files[["qfilt_lq_paf"]]))) " NULL" else str_c(" ${qfilt_lq_paf[", SGE_TASK_ID, "]}")),
      " ${intra_merge_corr[", SGE_TASK_ID, "]}",#outputs
      " ${intra_merge_gtf[", SGE_TASK_ID, "]}",
      " ${variant_local[", SGE_TASK_ID, "]}",
      " ${chimera_paf[", SGE_TASK_ID, "]}",
      " ", qslot("intra_merge"),
      " TRUE", #mono multi merge
      " FALSE" #collapse mono range
    )
  )

pl_inter_merge <- 
  qsub_template(
    parallel = parallel_option_req("inter_merge"),
    script_path = "inter_merge",
    directory = directory_option(out = dir_logs[["merge"]]),
    str_c(
      Rcmd, lib$inter_merge.R,
      lib$isocan,
      lib$merge_iso_commonpart.R,
      IO_summary_file,
      "intra_merge_corr",
      "intra_merge_gtf",
      file$inter_merge_corr,
      file$inter_merge_gtf,
      "fl_count",
      qslot("inter_merge"),
      "TRUE", #mono multi merge
      "FALSE", #collapse mono range
      sep = " "
    )
  )

# problem with conda container; BSgenome only 1.60.0 and rlang dynamic libraries cannot be loaded
# possible memory problem; setting threshold to 64GB temporarily alleviates issue
pl_gtf2fasta <- 
  qsub_template(
    parallel = parallel_option_req("gtf2fasta"),
    script_path = "gtf2fasta",
    directory = directory_option(out = dir_logs[["merge"]]),
    str_c(Rcmd, lib$gtf2fa.R, file$inter_merge_gtf, file$genome_pkg, file$inter_merge_fa, sep = " ")
  )

####salmon####
pl_salmon_idx_init <- 
  qsub_template(
    parallel = parallel_option_req("salmon_index"),
    script_path = "salmon_idx_init",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c("## DEPRECATION NOTICE: RNA-SEQ INDEX CREATION HAS BEEN ALLOCATED TO A DIFFERENT PART OF THE WORKFLOW", 
	  sprintf("## INDEX LOCATION: %s", file$salmon_index_pre_sqanti),
          "## THIS SCRIPT IS ONLY LEFT HERE AS A REFERENCE",  sep = "\n"))
#    str_c(soft$salmon, "index",
#          "-t", file$inter_merge_fa,
#          "-i", file$salmon_index_pre_sqanti, 
#          "--keepDuplicates", #added. 22samples contained only 2 duplicated isoform-pairs which were very similar with but 1 or ~10 bp diff against each other
#          sep = " ")

if (args$no_salmon) {
	message('Option "--no-salmon" specified; skipping salmon quant.')
	pl_salmon_execute_check <- "#"
} else {
        pl_salmon_execute_check <- ""
}

#pl_salmon_quant_init <- 
#  qsub_template(
#    parallel = parallel_option_req("salmon_quant"),
#    arrayjob = ary,
#    script_path = "salmon_quant_init",
#    directory = directory_option(out = dir_logs[["salmon"]]),
#    str_c(pl_salmon_execute_check, "mkdir -p ${sample_dir[", SGE_TASK_ID, "]}/", salmon_dirname$pre_sqanti, " && ", soft$salmon, " quant",
#          " -i ", file$salmon_index_pre_sqanti,
#          " -l A",
#          " -1 ${short_copy_1[", SGE_TASK_ID, "]}",
#          " -2 ${short_copy_2[", SGE_TASK_ID, "]}",
#          " -p ", qslot("salmon_quant"), 
#          " --gcBias --seqBias --validateMappings",
#          " -o ${sample_dir[", SGE_TASK_ID, "]}/", salmon_dirname$pre_sqanti)
#  )

# 6/1/23
pl_salmon_quant_init <-
  qsub_template(
    parallel = parallel_option_req("minimap2"),
    script_path = "salmon_quant_init",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(pl_salmon_execute_check,
          str_c(Rcmd, lib$subread_quant.R, 
                      file$salmon_index_pre_sqanti, IO_summary_file, file$inter_merge_gtf, salmon_dirname$pre_sqanti, 
                      "intermerge", get_req(x = "minimap2", r = "memory"), get_req(x = "minimap2", r = "slot"), genome_fa, sep = " ")
  )

pl_merge_salmon_init <- 
  qsub_template(
    parallel = parallel_option_req("merge_salmon"),
    script_path = "merge_salmon_init",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(pl_salmon_execute_check, Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_pre_sqanti TPM", file$salmon_tpm_pre_sqanti,
          sep = " "),
    str_c(pl_salmon_execute_check, Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_pre_sqanti NumReads", file$salmon_count_pre_sqanti,
          sep = " ")
  )

####SQANTI####
sq_header_cupcake <- sprintf("PYTHONPATH=%s", paste(unique(c(args$cupcake, unlist(strsplit(Sys.getenv("PYTHONPATH"), ":")))), collapse = ":"))
sq_qsub_param <- str_c(qsub_header_param, "\n#conda activate musta2\n", sprintf("export PYTHONPATH=%s/cDNA_Cupcake/sequence\n", .this_dir))
fl_count_all <- sprintf("%s/fl_count_all.txt", sample_dir_root) # current SQANTI3 version apparently only accepts a single FL_COUNT file

if (args$proxy != '') {
	sq_proxy_param <- sprintf("HTTP_PROXY=%s", args$proxy)
} else {
	sq_proxy_param <- ''
}

# str_c(input_files$fl_count, collapse = ","), #comma_separated list of FL file
## modified 6/8 to account for the lack of expression matrix from short read
if (args$no_salmon) {
        # no shortread processing
	message('Option "--no-salmon" specified; skipping salmon expression quantification for SQANTI3.')
        pl_sqanti_shortread_arg <- ""
} else {
	pl_sqanti_shortread_arg <- sprintf("-e %s", file$salmon_tpm_pre_sqanti)
}
if (file.exists(file$sqanti_class)) {
	sq_init_qc_run <- '#'
	message('Skipping initial SQANTI3 classification as output file appears to be present...')
} else {
	sq_init_qc_run <- ''
}
if (file.exists(file$sqanti_filter_result)) {
	sq_init_filter_run <- '#'
	message('Skipping SQANTI3 filtering step as output file appears to be present...')
} else {
	sq_init_filter_run <- ''
}
# to-do: modify lib source to account for addition of fl_count_all.R, 6/14/2022
if (args$sqanti2_mode) {
	sqanti3_filter_addl_param <- "-j 0.5 -i 80"
} else {
	sqanti3_filter_addl_param <- "-f TRUE"
}
pl_sqanti <- 
  qsub_template(
    parallel = parallel_option_req("sqanti"),
    script_path = "sqanti",
    directory = directory_option(out = dir_logs[["sqanti"]]),
    other_req = sq_qsub_param,
    "echo '======   MERGING FL COUNTS  ======'",
    "echo '======   MERGING FL COUNTS  ======' 1>&2",
    paste(Rcmd, lib$fl_count_all.R, fl_count_all, paste(input_files$fl_count, collapse = ","), sep = " "), 
    "echo '====== STARTING SQANTI3 QC  ======'",
    "echo '====== STARTING SQANTI3 QC  ======' 1>&2",
    str_c(sq_init_qc_run, sq_header_cupcake, args$python3, soft$sqanti_qc.py, "--force_id_ignore",
          pl_sqanti_shortread_arg, # modified 1/13/22 (see above)
          "-d", dir$merge, "-o", sqanti_prefix,
          ifelse(!is.null(input_files$SJ) && !any(is.na(input_files$SJ)), paste0("-c ", str_c(input_files$SJ, collapse = ",")), ""), #comma_separated list of SJ.out.tab
          "-fl", fl_count_all,
          file$inter_merge_gtf,
          ref_gtf,
          genome_fa,
          sep = " "), #isoforms.gtf, ref.gtf, ref.fa
    "echo '====== APPLY SQANTI3 FILTER ======'",
    "echo '====== APPLY SQANTI3 FILTER ======' 1>&2",
    str_c(sq_init_filter_run, sq_proxy_param, sq_header_cupcake, args$python3, soft$sqanti_filter.py,
          "ML", sqanti3_filter_addl_param, "-d", file$sqanti_filter_supplement_folder, #${out_dir}/${prefix}_.
          file$sqanti_class,
          sep = " "),
    str_c(sq_init_filter_run, "mv -f",
          paste0(file$sqanti_filter_supplement_folder, "/", sqanti_prefix, "_MLresult_classification*"),
          dir$merge,
          sep = " "), #move main results to ../ dir.
    "echo '====== UPDATING SQ ISOFORMS ======'",
    "echo '====== UPDATING SQ ISOFORMS ======' 1>&2",
    str_c(Rcmd, lib$update_post_sqanti.R,
          lib$isocan,
          IO_summary_file, 
          file$inter_merge_gtf,
          file$inter_merge_corr,
          file$inter_merge_fa,
          file$sqanti_filter_result,
          file$sqanti_gtf,
          file$sqanti_corr,
          file$sqanti_fa,
          file$tx_count_unique,
          file$tx_existence,
          sep = " ")
  )

####re salmon####
# see above for the execute checkbit
pl_salmon_index <- 
  qsub_template(
    parallel = parallel_option_req("salmon_index"),
    script_path = "salmon_index",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(soft$salmon, "index",
          "-t", file$sqanti_fa,
          "-i", file$salmon_index_post_sqanti,
          sep = " ")
  )

pl_salmon_quant <- 
  qsub_template(
    parallel = parallel_option_req("salmon_quant"),
    arrayjob = ary,
    script_path = "salmon_quant",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(pl_salmon_execute_check, 
          "mkdir -p ${sample_dir[", SGE_TASK_ID, "]}/", salmon_dirname$post_sqanti, " && ",
          soft$salmon, " quant",
          " -i ", file$salmon_index_post_sqanti,
          " -l A",
          " -1 ${short_copy_1[", SGE_TASK_ID, "]}",
          " -2 ${short_copy_2[", SGE_TASK_ID, "]}",
          " -p ", qslot("salmon_quant"), 
          " --gcBias --seqBias --validateMappings",
          " -o ${sample_dir[", SGE_TASK_ID, "]}/", salmon_dirname$post_sqanti
    )
  )

pl_merge_salmon <- 
  qsub_template(
    parallel = parallel_option_req("merge_salmon"),
    script_path = "merge_salmon",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(pl_salmon_execute_check, Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_post_sqanti TPM", file$salmon_tpm_post_sqanti,
          sep = " "),
    str_c(pl_salmon_execute_check, Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_post_sqanti NumReads", file$salmon_count_post_sqanti,
          sep = " ")
  )


####reSQANTI####
## modified 1/13 to account for the lack of expression matrix from short read
if (args$no_salmon) {
        # no shortread processing
	message('Option "--no-salmon" specified; skipping salmon expression quantification for SQANTI3.')
        pl_sqanti_shortread_redo_arg <- ""
} else {
        pl_sqanti_shortread_redo_arg <- sprintf("-e %s", file$salmon_tpm_post_sqanti)
}
pl_re_sqanti_qc <- 
  qsub_template(
    parallel = parallel_option_req("sqanti"),
    script_path = "re_sqanti_qc",
    directory = directory_option(out = dir_logs[["sqanti"]]),
    other_req = sq_qsub_param,
    str_c(sq_header_cupcake, args$python3, soft$sqanti_qc.py, "--force_id_ignore",
	  pl_sqanti_shortread_redo_arg, # modified 1/13
          "-d", dir$merge, "-o", re_sqanti_prefix, #${out_dir}/${prefix}_classification.txt will be generated.
          ifelse(!is.null(input_files$SJ) && !any(is.na(input_files$SJ)), paste0("-c ", str_c(input_files$SJ, collapse = ",")), ""), #comma_separated list of SJ.out.tab
          "-fl", fl_count_all, #str_c(input_files$fl_count, collapse = ","), #comma_separated list of FL file
          file$sqanti_gtf,
          ref_gtf,
          genome_fa,
          sep = " "), #isoforms.gtf, ref.gtf, ref.fa
    str_c("mv",
          path(dir$merge, paste0(re_sqanti_prefix, "_*")),
          dir$result,
          sep = " ")
  )

####salmon ref gtf####

pl_salmon_index_ref <- 
  qsub_template(
    parallel = parallel_option_req("salmon_index"),
    script_path = "salmon_index_ref",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(Rcmd, lib$gtf2fa.R, ref_gtf, file$genome_pkg, file$ref_fa, sep = " "),
    str_c(soft$salmon, "index",
          "-t", file$ref_fa,
          "-i", file$salmon_index_ref,
          sep = " ")
  )

pl_salmon_quant_ref <- 
  qsub_template(
    parallel = parallel_option_req("salmon_quant"),
    arrayjob = ary,
    script_path = "salmon_quant_ref",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(soft$salmon, " quant",
          " -i ", file$salmon_index_ref,
          " -l A",
          " -1 ${short_copy_1[", SGE_TASK_ID, "]}",
          " -2 ${short_copy_2[", SGE_TASK_ID, "]}",
          " -p ", qslot("salmon_quant"), 
          " --gcBias --seqBias --validateMappings",
          " -o ${sample_dir[", SGE_TASK_ID, "]}/", salmon_dirname$ref
    )
  )

pl_merge_salmon_ref <- 
  qsub_template(
    parallel = parallel_option_req("merge_salmon"),
    script_path = "merge_salmon_ref",
    directory = directory_option(out = dir_logs[["salmon"]]),
    str_c(Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_ref TPM", file$salmon_tpm_ref,
          sep = " "),
    str_c(Rcmd, lib$merge_salmon_iso.R,
          IO_summary_file, "salmon_exp_ref NumReads", file$salmon_count_ref,
          sep = " ")
  )

####Link Original####
pl_filter_variant <- #filter
  qsub_template(
    parallel = parallel_option_req("filter_variant"),
    script_path = "filter_variant",
    directory = directory_option(out = dir_logs[["others"]]),
    str_c(Rcmd, lib$filter_variants.R,
          lib$isocan,
          file$genome_pkg,
          IO_summary_file,
          "variant_local",
          "mutation",
          file$sqanti_corr,
          file$variant_filtered,
          qslot("filter_variant"),
          sep = " ")
  )

pl_link_original_range <- 
  qsub_template(
    parallel = parallel_option_req("link_original_range"),
    script_path = "link_original_range",
    directory = directory_option(out = dir_logs[["others"]]),
    str_c(Rcmd, lib$link_original_range.R,
          IO_summary_file,
          file$sqanti_corr,
          "variant_local",
          file$original_range_multiex,
          sep = " ")
  )
####Fusion Screening####
pl_fusion_bind <- 
  qsub_template(
    parallel = parallel_option_req("fusion_bind"),
    script_path = "fusion_bind",
    directory = directory_option(out = dir_logs[["fusion"]]),
    str_c(Rcmd, lib$bind_chimera.R,
          IO_summary_file,
          "chimera_paf",
          file$chimera_binded,
          sep = " ")
  )

pl_fusion_parse <- 
  qsub_template(
    parallel = parallel_option_req("fusion_parse"),
    script_path = "fusion_parse",
    directory = directory_option(out = dir_logs[["fusion"]]),
    str_c(Rcmd, lib$parse_chimera.R,
          lib$isocan,
          file$genome_pkg,
          file$chimera_binded,
          IO_summary_file,
          "sv",
          ref_gtf,
          file$re_sqanti_junc,
          file$fusion_rds,
          file$fusion_gene_fragment_gtf,
          "100",  #TODO parameter config
          "100000",
          qslot("fusion_parse"),
          sep = " ")
  )

pl_fusion_qsqanti_qc <- 
  qsub_template(
    parallel = parallel_option_req("fusion_sqanti"),
    script_path = "fusion_sqanti",
    directory = directory_option(out = dir_logs[["fusion"]]),
    other_req = sq_qsub_param,
    paste0("if [ $( cat ", file$fusion_gene_fragment_gtf, " | wc -w ) -gt 0 ]; then"),
    # run only when fusion_gtf contains at least one fusion.
    str_c(sq_header_cupcake, args$python3, soft$sqanti_qc.py, "--force_id_ignore -g", #gtf
          "-d", dir$fusion, 
          "-o", fusion_sqanti_prefix,
          file$fusion_gene_fragment_gtf,
          ref_gtf,
          genome_fa,
          sep = " "),
    "fi"
  )

pl_fusion_summary <- 
  qsub_template(
    parallel = parallel_option_req("fusion_summary"),
    script_path = "fusion_summary",
    directory = directory_option(out = dir_logs[["fusion"]]),
    paste0("if [ $( cat ", file$fusion_gene_fragment_gtf, " | wc -w ) -gt 0 ]; then"),
    # run only when fusion_gtf contains at least one fusion.
    str_c(Rcmd, lib$summarise_chimera.R,
          lib$isocan,
          IO_summary_file, 
          file$fusion_rds,
          ref_gtf,
          file$fusion_sqanti_class,
          file$fusion_summary,
          file$fusion_gtf,
          file$fusion_sv,
          sep = " "),
    "fi"
  )

####report####
pl_report_number <- 
  qsub_template(
    parallel = parallel_option_req("report_number"),
    script_path = "report_number",
    directory = directory_option(out = dir_logs[["others"]]),
    str_c(Rcmd, lib$report_number.R,
          IO_summary_file, 
          file$inter_merge_gtf,
          file$sqanti_gtf,
          file$number_report_samples,
          file$number_report_merged,
          qslot("report_number"),
          sep = " ")
  )

####clean up####
pl_cleanup_interleave <- 
  qsub_template(
    script_path = "cleanup_interleave",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -f", input_files$interleave, sep = " ", collapse = "\n"),
    str_c("rm -f", input_files$interleave_lordectmp, sep = " ", collapse = "\n")
  )

pl_cleanup_cp_long_hq <- 
  qsub_template(
    script_path = "cleanup_cp_long_hq",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -f", input_files$long_copy_hq, sep = " ", collapse = "\n")
  )

pl_cleanup_cp_long_lq <- 
  qsub_template(
    script_path = "cleanup_cp_long_lq",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -f", input_files$long_copy_lq, sep = " ", collapse = "\n")
  )

pl_cleanup_cp_short <- 
  qsub_template(
    script_path = "cleanup_cp_short",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -f", input_files$short_copy_1, sep = " ", collapse = "\n"),
    str_c("rm -f", input_files$short_copy_2, sep = " ", collapse = "\n")
  )

pl_cleanup_minimap2 <- 
  qsub_template(
    script_path = "cleanup_minimap2",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -f", file$genome_mmi, sep = " ", collapse = "\n"),
    str_c("rm -f", input_files$minimap2_hq_sam, sep = " ", collapse = "\n"),
    ifelse(!is.null(input_files[["minimap2_lq_sam"]]), str_c("rm -f", input_files$minimap2_lq_sam, sep = " ", collapse = "\n"), ""),
    str_c("rm -f", input_files$minimap2_hq_paf, sep = " ", collapse = "\n"),
    ifelse(!is.null(input_files[["minimap2_lq_paf"]]), str_c("rm -f", input_files$minimap2_lq_paf, sep = " ", collapse = "\n"), "")
  )

pl_cleanup_salmon_init <- 
  qsub_template(
    script_path = "cleanup_salmon_init",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -rf", file$salmon_index_pre_sqanti, sep = " ", collapse = "\n")
  )

pl_cleanup_sqanti <- 
  qsub_template(
    script_path = "cleanup_sqanti",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -rf", fs::path(input_files$sample_dir, salmon_dirname$pre_sqanti), sep = " ", collapse = "\n"),
    str_c("rm -f", file$inter_merge_corrected, sep = " ", collapse = "\n"),
    str_c("rm -f", file$salmon_tpm_pre_sqanti, sep = " ", collapse = "\n"),
    str_c("rm -f", file$salmon_count_pre_sqanti, sep = " ", collapse = "\n")
  )

pl_cleanup_salmon <- 
  qsub_template(
    script_path = "cleanup_salmon",
    directory = directory_option(out = dir_logs[["file_management"]]),
    str_c("rm -rf", path(dir$merge, "salmon"), sep = " ", collapse = "\n")
  )

## Revised "Clean-Up" Step
qsub_openblas_param <- sprintf("OPENBLAS_NUM_THREADS=%d", qslot("filter_variant"))
dupl.rm <- 0 # not yet implemented as argument
if ((file.exists(file$salmon_tpm_post_sqanti)) & (file.exists(file$salmon_count_post_sqanti))) {
	qsub_cleanup_salmon_param <- paste(file$salmon_tpm_post_sqanti, file$salmon_count_post_sqanti)
} else {
	qsub_cleanup_salmon_param <- ""
}
pl_cleanup_re_sqanti <- 
  qsub_template(
    parallel = parallel_option_req("filter_variant"),
    script_path = "cleanup_re_sqanti",
    directory = directory_option(out = dir_logs[["file_management"]]),
    "echo '===== COUNTING MUSTA ISOFORMS ===='",
    "echo '===== COUNTING MUSTA ISOFORMS ====' 1>&2",
    str_c(qsub_openblas_param, Rcmd, lib$recalc_pbcount.R,
          args$input,
          file$sqanti_corr,
          file$tx_existence,
          qslot("filter_variant"),
	  0, # dupl.rm
          qsub_cleanup_salmon_param, sep = " "),
    "echo '==== ANNOTATING TRANSCRIPTOME ===='",
    "echo '==== ANNOTATING TRANSCRIPTOME ====' 1>&2",    
    str_c(qsub_openblas_param, Rcmd, lib$annotate_output_gtf.R,
          file$sqanti_gtf,
          args$gtf, 
          qslot("filter_variant"), sep = " "),
    str_c("rm -f", file$re_sqanti_corrected, sep = " ", collapse = "\n")
  )


## DEPRECATED FUNCTIONS: TO BE REMOVED BEFORE PUBLIC RELEASE
pl_interleave <- qsub_template(
	parallel = parallel_option_req("lordec"),
	arrayjob = ary,
	script_path = "interleave",
	directory = directory_option(out = dir_logs[["interleave"]]),
	paste0("#", lib$interleave_fastq, " ${short_copy_1[", SGE_TASK_ID, "]} ${short_copy_2[", SGE_TASK_ID, "]} > ${interleave[", SGE_TASK_ID, "]}")
)
pl_lordec_build <- qsub_template(
	parallel = parallel_option_req("lordec_build"),
	script_path = "lordec_build",
	directory = directory_option(out = dir_logs[["lordec"]]),
	str_c("#",soft$lordec_build_sr_graph, "-T", qslot("lordec_build"), "-2", input_files$interleave, "-k 19 -s 3 -g", input_files$interleave_lordectmp, sep = " ")
)
pl_lordec_hq <- qsub_template(
	parallel = parallel_option_req("lordec"), arrayjob = ary, script_path = "lordec_hq", directory = directory_option(out = dir_logs[["lordec"]]),
	str_c("#",soft$lordec_correct, " -T ", qslot("lordec"), " -i ${long_copy_hq[", SGE_TASK_ID, "]}", " -k 19 -s 3 -2 ${interleave[", SGE_TASK_ID, "]}", " -o ${lordec_hq[", SGE_TASK_ID, "]}")
)
pl_lordec_lq <- qsub_template(
	parallel = parallel_option_req("lordec"), arrayjob = ary, script_path = "lordec_lq", directory = directory_option(out = dir_logs[["lordec"]]),
	str_c("#", soft$lordec_correct, " -T ", qslot("lordec"), " -i ${long_copy_lq[", SGE_TASK_ID, "]}",
		" -k 19 -s 3 -2 ${interleave[", SGE_TASK_ID, "]}", " -o ${lordec_lq[", SGE_TASK_ID, "]}")
)

