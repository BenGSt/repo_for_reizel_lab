
# -ncpus: number of cores to use (default: 6)
# -samplesheet: path to a samplesheet csv file. (default:./samplesheet.csv)
# -output_dir: path to save bed files and other outputs to (defualt: ./dmrseq_output)
# -min_qval: minimum qval for a DMR to be considered significant (default: 0.05)
# -save_full_dmr_data: flag, save the full DMR data including statistics as a TSV file
# -save_rds: flag, save the resulting objects as RDS files (can be loaded individually in R)
# -save_rdata: flag, save the R environment (can be loaded with load())
# [dmrseq_args]: additional arguments passed to dmrseq

Rscript find_dmrs_dmrseq.R \
  -ncpus 6 \
  -samplesheet ./example_samplesheet.csv \
  -output_dir ./example_output \
  -min_qval 0.05 \
  cutoff 0.1
