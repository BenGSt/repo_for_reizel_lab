#!/bin/bash

write_rrbs_jobs_args() # <1 for single end or 2 for piared end> <raw_data_dir>
{
  for sample_name in $(ls $2); do
    if [[ $1 == 1 ]]; then
      echo $sample_name, -output_dir $sample_name -single_end -input_fastq_file $(realpath $2/$sample_name/*)
    else
      echo $sample_name, -output_dir $sample_name -paired-end -paired_input_fastq_files $(realpath $2/$sample_name/*)
    fi
  done
}

write_condor_submition_files()
{
  PATH_TO_EXECUTABLES=/srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag
submission_file_names=(
                trim_illumina_adaptors_jobs.sub \
                trim_diversity_adaptors_jobs.sub \
                align_jobs.sub
                meth_calling_jobs.sub \
                make_tiles_jobs.sub
)
executables=(
            trim_illumina_adaptors.sh \
            trim_diversity_adaptors.sh \
            align_to_genome.sh \
            methylation_calling.sh \
            combine_methylation_coverage_to_tiles.sh \
)
cores=(4 1 10 1 1)
rams=(500MB 500MB 32GB 10GB 500MB)

for i in "${!submission_file_names[@]}"; do

  cat << EOF > ${submission_file_names[$i]}
executable = $PATH_TO_EXECUTABLES/${executables[$i]}
Initialdir = $(pwd)
request_cpus = ${cores[$i]}
RequestMemory = ${rams[$i]}
Arguments = \$(name) \$(args)
universe = vanilla
log = ./\$(name)/condor_logs/\$(name)_align.log
output = ./\$(name)/condor_logs/\$(name)_align.out
error = ./\$(name)/condor_logs/\$(name)_align.error
queue name, args from rrbs_jobs.args
EOF

done
}

