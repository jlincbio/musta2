drake_txt <- ''

plan_add(
  TRUE, 
  '#in "()", write direct dependencies. They will not be used actually.
  BSgenome = pl_BSgenome(),
')

plan_add(
  args$copy_input,
  '#copy
  cp_1 = pl_cp_1(),
  cp_2 = pl_cp_2(),
')

plan_add(
  args$copy_input && args$use_lq,
  '  cp_lq = pl_cp_lq(),
')

plan_add(
  args$use_lq,
  '
  merge_fa = pl_merge_fa(c(lordec_hq, lordec_lq)),
')

plan_add(
  args$copy_input,
  '
  cleanup_cp_long_hq = pl_cleanup_cp_long_hq(lordec_hq),
')

plan_add(
  args$copy_input && args$use_lq,
  '  cleanup_cp_long_lq = pl_cleanup_cp_long_lq(lordec_lq),
')

plan_add(
  TRUE,
  '#minimap2
  minimap2_make_mmi = pl_minimap2_make_mmi(),
  minimap2_hq_paf = pl_minimap2_hq_paf(c(lordec_hq, minimap2_make_mmi)),
  minimap2_hq_sam = pl_minimap2_hq_sam(c(lordec_hq, minimap2_make_mmi)),
')

plan_add(
  args$use_lq,
  '
  minimap2_lq_paf = pl_minimap2_lq_paf(c(lordec_lq, minimap2_make_mmi)),
  minimap2_lq_sam = pl_minimap2_lq_sam(c(lordec_lq, minimap2_make_mmi)),
')

plan_add(
  !args$keep,
  '
  cleanup_minimap2 = pl_cleanup_minimap2(report_number),
  '
)

plan_add(
  TRUE,
  str_glue(
    '
    qual_filter = pl_qual_filter({dep}),
    samtools_hq = pl_samtools_hq(qual_filter),
    ',
    dep = if (args$use_lq) 'c(minimap2_hq_paf, minimap2_lq_paf, minimap2_hq_sam, minimap2_lq_sam)' else 'c(minimap2_hq_paf, minimap2_hq_sam)'
  )
)

plan_add(
  args$use_lq,
  '
  samtools_lq = pl_samtools_lq(qual_filter),
  merge_sam = pl_merge_sam(c(samtools_hq, samtools_lq)),
  samtools_merge = pl_samtools_merge(merge_sam),
')

plan_add(
  TRUE,
  '#merge
  intra_merge = pl_intra_merge(qual_filter),
  fusion_bind = pl_fusion_bind(intra_merge),
  inter_merge = pl_inter_merge(intra_merge),
  gtf2fasta = pl_gtf2fasta(c(inter_merge, BSgenome)),
#salmon_init
  salmon_idx_init = pl_salmon_idx_init(gtf2fasta),
  salmon_quant_init = pl_salmon_quant_init(salmon_idx_init),
  merge_salmon_init = pl_merge_salmon_init(salmon_quant_init),
#sqanti
  sqanti = pl_sqanti(merge_salmon_init),
#salmon
  salmon_index = pl_salmon_index(sqanti),
  salmon_quant = pl_salmon_quant(salmon_index),
  merge_salmon = pl_merge_salmon(salmon_quant),
#re_sqanti
  re_sqanti_qc = pl_re_sqanti_qc(merge_salmon),
#others
  link_original_range = pl_link_original_range(sqanti),
  report_number = pl_report_number(sqanti)' # EOF, no comma.
)

plan_add(
  !args$keep,
  ',
  cleanup_salmon_init = pl_cleanup_salmon_init(salmon_quant_init),
  cleanup_sqanti = pl_cleanup_sqanti(report_number),
  cleanup_salmon = pl_cleanup_salmon(salmon_quant),
  cleanup_re_sqanti = pl_cleanup_re_sqanti(re_sqanti_qc)
  '
)

plan_add(
  !all(is.na(input_files$mutation)),
  ',
#variant
  filter_variant = pl_filter_variant(sqanti)' # implicitly relies on BSgenome
)

plan_add(
  !all(is.na(input_files$sv)),
  ',
#fusion
  fusion_parse = pl_fusion_parse(c(fusion_bind, sqanti)),
  fusion_sqanti_qc = pl_fusion_sqanti_qc(fusion_parse),
  fusion_summary = pl_fusion_summary(fusion_sqanti_qc)' # fusion_parse implicitly relies on BSgenome
)

plan_add(
  args$salmon_ref,
  ',
#ref
  salmon_index_ref = pl_salmon_index_ref(BSgenome),
  salmon_quant_ref = pl_salmon_quant_ref(salmon_index_ref),
  merge_salmon_ref = pl_merge_salmon_ref(salmon_quant_ref)'
)

plan_add(
  args$copy_input,
  str_glue(
  ',
  cleanup_cp_short = pl_cleanup_cp_short({dep})
',
  dep = if (args$salmon_ref) "c(salmon_quant, salmon_quant_ref)" else "salmon_quant")
)

drake_txt <- str_replace_all(drake_txt, ",#", ",\n#")
drake_txt <- str_replace_all(drake_txt, "[:space:]*\n[:space:]*", "\n  ")
drake_txt <- str_replace_all(drake_txt, "  #", "#")
drake_txt <- str_replace_all(drake_txt, "([:space:]|\n)+,", ",")

drake_txt <- paste0(
  'df_plan <- drake_plan(
  ', drake_txt, '\n)'
)
