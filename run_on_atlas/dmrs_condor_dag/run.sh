#!/bin/bash

main()
{
  if [[ $# -lt 2 ]]; then
    echo USAGE: $0
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


write_dmr_jobs_args() # <1 for single end or 2 for piared end> <raw_data_dir>
{
 e.g 56n vs 57n, --meth_call_files_dir /storage/bfe_reizel/bengst/analyzed_data/Nadav_betaRepl/56n_vs_57n/cov_files --samp_ids \"56n_1 56n_2 57n_1 57n_2\" --treatments \"1 1 0 0\" --pipeline bismarkCoverage --output_dir /storage/bfe_reizel/bengst/analyzed_data/Nadav_betaRepl/56n_vs_57n --known_genes_file /storage/bfe_reizel/bengst/analyzed_data/mm10KnownGenes.bed
}


write_dmr_jobs_sub_file()
{
    cat << EOF > dmr_jobs.sub
executable = dmr_job.sh
log = ./\$(name)/\$(name).log
output = ./\$(name)/\$(name).out
error = ./\$(name)/\$(name).error
request_cpus = 1
Initialdir = $(pwd)
Arguments = \$(args)
RequestMemory = 4GB
universe = vanilla
queue name,args from dmr_jobs.args
EOF
}

write_condor_submition_files()
{
PATH_TO_EXECUTABLES=/srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag
submission_file_names=(
                dmr_jobs.sub
)
executables=(
            dmr_job.sh
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