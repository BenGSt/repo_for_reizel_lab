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


write_dmr_jobs_args()
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


write_heatmap_jobs_args()
{
 echo 97_vs_91, /storage/bfe_reizel/bengst/analyzed_data/KKTR-TargetingMafAMotifWithTet/dmrs_01.07.2022/all_samples_100bp_tiles_each_run_as_separate_samples.bed \
      97_vs_91 \
      --sample_names 91E-97B-91D-91A-97G-97A-97C-97F-91B-97E-91C-91A_2-97G_2-97A_2-97E_2-91C_2 > heatmap_jobs.args

 echo edit heatmap_jobs.args if you want to select part of the samples \(see /scripts_2022/make_heatmap.R\)
}


write_homer_jobs_args()
{
 echo 97_vs_91 > homer_jobs.args
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
RequestMemory = 1GB
universe = vanilla
queue name,args from heatmap_jobs.args
EOF
}


write_homer_jobs_sub_file()
{
  cat << EOF > homer_jobs.sub
executable = /srv01/technion/bengst/scripts/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/homer_job.sh
log = ./\$(name)/condor_logs/homer_\$(name).log
output = ./\$(name)/condor_logs/homer_\$(name).out
error = ./\$(name)/condor_logs/homer_\$(name).out
request_cpus = 1
Initialdir = $(pwd)
Arguments = \$(name)
RequestMemory = 4GB
universe = vanilla
queue name from homer_jobs.args
EOF
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