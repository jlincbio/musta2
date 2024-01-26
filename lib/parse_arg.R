parser <- ArgumentParser(prog = "musta2", description = sprintf("%s, v%s: Multi-Sample Transcriptome Assembly", .name, .version))

parser$add_argument('-v', '--version', action = 'version', version = .version)

#### INPUT ####
group_input <- parser$add_argument_group('INPUT')

group_input$add_argument(
  '-i', '--input', required = TRUE,
  help = 'Comma-separated (CSV) file describing sample IDs and input paths; at minimum this should contain the following header fields: "sample", "long_read_hq", "short_read_1", "short_read_2"; optional fields may include: "cluster_report", "long_read_lq", "SJ", "mutation", "sv"'
)

group_input$add_argument(
  '-t', '--gtf', required = TRUE,
  help = 'Reference transcriptome GTF, e.g. Encode Release 28'
)

group_input$add_argument(
  '-g', '--genome', required = TRUE,
  help = 'Reference genome FASTA'
)

#### OUTPUT ####
group_output <- parser$add_argument_group('OUTPUT')

group_output$add_argument(
  '-o', '--output', default = 'output',
  help = 'Output location; defaults to "output"'
)

group_output$add_argument(
  '-p', '--prefix', default = 'musta2',
  help = 'Result file prefix; defaults to "musta2"'
)

#### CONFIGURATION ####
group_configuration <- parser$add_argument_group('CONFIGURATION')

group_configuration$add_argument(
  '--copy-input', action = "store_true",
  help = 'Copy input files to a temporary location instead of using them as is; useful when inputs are stored on network shares or external drives with slow I/O.'
)

group_configuration$add_argument(
  '--no-salmon', action = "store_true",
  help = 'Perform MuSTA analysis with long-read data only; skips Salmon quantification steps.'
)

group_configuration$add_argument(
  '--sqanti2-mode', action = "store_true",
  help = 'Sets the default SQANTI classification threshold probability to 0.7 and intrapriming length to 80 bp (i.e., SQANTI2 defaults).'
)

group_configuration$add_argument(
  '--use-lq', action = "store_true",
  help = 'Use both low and high quality IsoSeq reads; "long_read_lq" field is required.'
)

group_configuration$add_argument(
  '--map-qual', type = 'integer', default = 50L,
  help = 'Mapping quality threshold (default: 50)'
)

group_configuration$add_argument(
  '--salmon-ref', action = "store_true",
  help = 'Run salmon against a reference GTF as well as against a MuSTA-derived GTF.'
)

group_configuration$add_argument(
  '--keep', action = "store_true", 
  help = 'Keep intermediate files'
)


#### DEPENDENCIES ####
group_ext <- parser$add_argument_group('DEPENDENCIES', 'Note: most of these should have already been preconfigured during the installation process.')
group_ext$add_argument(
  '--sqanti-filter', default = sprintf("%s/sqanti3_filter.py", .this_dir),
  help = 'Path to "sqanti3_filter.py" in SQANTI3'
)

group_ext$add_argument(
  '--env', default = 'musta2', help = 'Name of Miniconda Environment; defaults to "musta2"'
)

group_ext$add_argument(
  '--base', default = Sys.getenv("MUSTA2_BASEDIR"), help = 'Location of MuSTA2 library'
)

group_ext$add_argument(
  '--sqanti-qc', default = sprintf("%s/sqanti3_qc.py", .this_dir),
  help = 'Path to "sqanti3_qc.py" in SQANTI3'
)

group_ext$add_argument(
  '--cupcake', default = sprintf("%s/cDNA_Cupcake/sequence", .this_dir),
  help = 'Path to PacBio "Cupcake" sequence scripts for SQANTI3'
)  

group_ext$add_argument(
  '--python3', default = Sys.which("python3"),
  help = 'Path to Python 3'
)

group_ext$add_argument(
  '--proxy', default = "", help = 'http proxy server settings.'
)

#group_ext$add_argument(
#  '--lordec-build-SR-graph', default = 'lordec-build-SR-graph',
#  help = 'Deprecated (no functionality)'
#)

#group_ext$add_argument(
#  '--lordec-correct', default = 'lordec-correct',
#  help = 'Deprecated (no functionality)'
#)

group_ext$add_argument(
  '--minimap2', default = Sys.which("minimap2"),
  help = 'Path to minimap2'
)

group_ext$add_argument(
  '--samtools', default = Sys.which("samtools"),
  help = 'Path to samtools'
)

group_ext$add_argument(
  '--salmon', default = Sys.which("salmon"),
  help = 'Path to salmon'
)
group_ext$add_argument(
  '--seqkit', default = Sys.which("seqkit"),
  help = 'Path to seqkit'
)

#### Miscellaneous ####
group_miscellaneous <- parser$add_argument_group('MISCELLANEOUS')

group_miscellaneous$add_argument(
  '--force', action = "store_true",
  help = 'Restart the entire pipeline rather than resuming from where the last run terminated.'
)

group_miscellaneous$add_argument(
  '--verbose', action = "store_true", 
  help = 'Verbose mode (for debugging)'
)

group_miscellaneous$add_argument(
  '--dry-run', action = "store_true",
  help = 'Run MuSTA2 in dry-run mode; a flight plan will be generated but MuSTA componants are not executed.'
)

args <- parser$parse_args()
