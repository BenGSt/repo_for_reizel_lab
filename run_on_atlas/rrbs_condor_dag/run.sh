#!/bin/bash

main()
{
  if [[ $# -lt 2 ]]; then
    echo USAGE: $0 \<1 for single end or 2 for piared end\> \<\raw_data_dir\>
    echo run from the directory you wish the output to be written to
    exit 1
  fi
  write_rrbs_jobs_args $1 $2
  write_condor_submition_files
  write_condor_dag
  make_dirs

  echo Submit the jobs by running: condor_submit_dag rrbs_jobs.dag
  echo Good Luck!
}


write_rrbs_jobs_args() # <1 for single end or 2 for piared end> <raw_data_dir>
{
  for sample_name in $(ls $2); do
    if [[ $1 == 1 ]]; then
      echo $sample_name, -output_dir $(pwd)/$sample_name -single-end -input_fastq_file $(realpath $2/$sample_name/*.fastq.gz) >> rrbs_jobs.args
    else
      echo $sample_name, -output_dir $(pwd)/$sample_name -paired-end -paired_input_fastq_files $(realpath $2/$sample_name/*.fastq.gz) >> rrbs_jobs.args
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
cores=(8 1 10 10 1)
rams=(300MB 10MB 16GB 3GB 3GB)

for i in "${!submission_file_names[@]}"; do

  cat << EOF > ${submission_file_names[$i]}
executable = $PATH_TO_EXECUTABLES/${executables[$i]}
Initialdir = $(pwd)
request_cpus = ${cores[$i]}
RequestMemory = ${rams[$i]}
Arguments = \$(args)
universe = vanilla
log = ./\$(name)/condor_logs/\$(name)_$(echo ${executables[$i]} | sed 's/.sh//').log
output = ./\$(name)/condor_logs/\$(name)_$(echo ${executables[$i]} | sed 's/.sh//').out
error = ./\$(name)/condor_logs/\$(name)_$(echo ${executables[$i]} | sed 's/.sh//').out
queue name, args from rrbs_jobs.args
EOF

done
}


write_condor_dag()
{
  cat << EOF > rrbs_jobs.dag
JOB trim_illumina trim_illumina_adaptors_jobs.sub
JOB trim_ovation trim_diversity_adaptors_jobs.sub
JOB align align_jobs.sub
JOB meth_call meth_calling_jobs.sub
JOB tiles make_tiles_jobs.sub

PARENT trim_illumina  CHILD trim_ovation
PARENT trim_ovation CHILD align
PARENT align CHILD meth_call
PARENT meth_call CHILD tiles
EOF
}

make_dirs()
{
  awk -F , '{print "mkdir -p "$1"/condor_logs"}' rrbs_jobs.args | bash
}

main "$@"