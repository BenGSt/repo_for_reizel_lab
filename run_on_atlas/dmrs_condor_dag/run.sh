#!/bin/bash

main()
{
#  if [[ $# -lt 2 ]]; then
#    echo USAGE: $0
#    echo first manually edit dmr_jobs.args - then run this script from the directory you wish the output to be written to
#    exit 1
#  fi

  if [[ -f dmr_jobs.args ]]; then
    write_heatmap_jobs_args
    write_homer_jobs_args

    write_dmr_jobs_sub_file
    write_heatmap_jobs_sub_file
    write_homer_jobs_sub_file
    write_condor_dag

    make_dirs

    echo Submit the jobs by running: condor_submit_dag dmr_pipline_jobs.dag
    echo Good Luck!
  else
    echo writing an example dmr_jobs.args, edit it then rerun this script.
    write_dmr_jobs_args
  fi
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
  #format: <name_for_condor_logs>,  <path to all_samples_100bp_tiles.bed> <sample_dir> [args for make_heatmap.R]
  all_samp_tiles=/storage/bfe_reizel/bengst/analyzed_data/KKTR-TargetingMafAMotifWithTet/dmrs_01.07.2022/all_samples_100bp_tiles_each_run_as_separate_samples.bed
  cat dmr_jobs.args | awk -F , -v all_samp_tiles=$all_samp_tiles 'match($0, /--samp_ids ([^ ]*)/, array)  {print $1",",  all_samp_tiles, $1, "--sample_names " array[1]}' > heatmap_jobs.args

 echo edit heatmap_jobs.args if you want to select part of the samples \(see /scripts_2022/make_heatmap.R\)
}


write_homer_jobs_args()
{
#format: sample_dir
 cat dmr_jobs.args | awk -F , '{print $1}' > homer_jobs.args

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
  cat << EOF > dmr_pipline_jobs.dag
JOB find_dmrs dmr_jobs.sub
JOB heatmap heatmap_jobs.sub
JOB homer  homer_jobs.sub

PARENT find_dmrs  CHILD heatmap homer
EOF
}


make_dirs()
{
  awk -F , '{print "mkdir -p "$1"/condor_logs"}' dmr_jobs.args | bash
}


main "$@"