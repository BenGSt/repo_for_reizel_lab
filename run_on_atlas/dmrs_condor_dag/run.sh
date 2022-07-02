#!/bin/bash

main()
{
  if [[ $# -lt 2 ]]; then
    echo USAGE: $0
    echo run from the directory you wish the output to be written to
    exit 1
  fi

  write_dmr_jobs_sub_file
  make_dirs

  echo Submit the jobs by running: condor_submit_dag rrbs_jobs.dag
  echo Good Luck!
}


write_dmr_jobs_args() # <1 for single end or 2 for piared end> <raw_data_dir>
{
 echo treated_vs_control, \
  --meth_call_files_dir /storage/bfe_reizel/bengst/analyzed_data/.../bismark_cov_files \
 --samp_ids id1-id2-id3-id4 \
 --treatments 1-1-0-0 \
 --pipeline bismarkCoverage \
 --output_dir /storage/bfe_reizel/bengst/analyzed_data/Nadav_betaRepl/56n_vs_57n \
 --known_genes_file /storage/bfe_reizel/bengst/genomic_reference_data/mm10KnownGenes.bed > dmr_jobs.args

 echo You\'re going to have to manually edit dmr_jobs.args
}

write_heatmap_jobs_args() # <1 for single end or 2 for piared end> <raw_data_dir>
{
 echo 97_vs_91, /storage/bfe_reizel/bengst/analyzed_data/KKTR-TargetingMafAMotifWithTet/dmrs_01.07.2022/all_samples_100bp_tiles_each_run_as_separate_samples.bed \
      97_vs_91 > heatmap_jobs.args

 echo edit heatmap_jobs.args if you want to select part of the samples \(see /scripts_2022/make_heatmap.R\)
}


write_dmr_jobs_sub_file()
{
    cat << EOF > dmr_jobs.sub
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/dmr_job.sh
log = ./\$(name)/condor_logs/dmrs_\$(name).log
output = ./\$(name)/condor_logs/dmrs_\$(name).out
error = ./\$(name)/condor_logs/dmrs_\$(name).out
request_cpus = 1
Initialdir = $(pwd)
Arguments = \$(args)
RequestMemory = 8GB
universe = vanilla
queue name,args from dmr_jobs.args
EOF
}


write_heatmap_jobs_sub_file()
{
  cat << EOF > heatmap_jobs.sub
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/heatmap_job.sh
log = ./\$(name)/condor_logs/heatmap_\$(name).log
output = ./\$(name)/condor_logs/heatmap_\$(name).out
error = ./\$(name)/condor_logs/heatmap_\$(name).out
request_cpus = 1
Initialdir = $(pwd)
Arguments = \$(args)
RequestMemory = 8GB
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
  awk -F , '{print "mkdir -p "$1"/condor_logs"}' dmr_jobs.args | bash
}

main "$@"