#!/bin/bash


#How to run:
# make sure you have the cov files you need in a dir
# run run.sh once
# edit dmr_jobs.args, heatmap_jobs.args - make sure heatmap_jobs.args sample_dir is output dir of dmr_job
# run run.sh 2nd time
# submit
#
###
main()
{
  if [[ $# -lt 1 ]]; then
    echo USAGE: $0 \<all_samples_100bp_tiles.bed\>
    exit 1
  fi

  all_samples_100bp_tiles=$1
  if [[ -f dmr_jobs.args ]]; then
    write_homer_jobs_args
    write_heatmap_jobs_args "$all_samples_100bp_tiles"

    write_dmr_jobs_sub_file
    write_heatmap_jobs_sub_file
    write_homer_jobs_sub_file
    write_condor_dag

    make_dirs

    echo Submit the jobs by running: condor_submit_dag dmr_pipline_jobs.dag
    echo Good Luck!
    echo
  else
    echo Writing example dmr_jobs.args file. Edit it, one line per job, then rerun this script.
    echo
    write_dmr_jobs_args
    echo
    echo
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
 --known_genes_file /storage/bfe_reizel/bengst/genomic_reference_data/mm10KnownGenes.bed \
 --meth_difference 25 > dmr_jobs.args

 echo You\'re going to have to manually edit dmr_jobs.args
 echo format: \<name_for_condor_logs\>, \<args for find_dmrs_methylkit.R\>
 echo Note that usually this is run as part of the rrbs pipline, and the name_for_condor_logs is the name of the output_dir.
 echo For other use cases, if name_for_condor_logs is not the name of the output_dir - an additional dir named name_for_condor_logs will be created with a condor_logs sub directory.
}


write_heatmap_jobs_args()
{
  ecgo writing heatmap_jobs.args
  echo heatmap_jobs.args format: \<name_for_condor_logs\>,  \<path to all_samples_100bp_tiles.bed\> \<sample_dir - output dir of dmr_job\> [args for make_heatmap.R]
  echo If all_samples_100bp_tiles.bed include more samples than you want to show in your heatmap you must edit heatmap_jobs.args
  echo  use: --sample_names \<ordered names of samples in all_samples_100bp_tiles.bed\> --include_samples_by_name \<samples to inclide\>
  echo e.g. --sample_names 56n-58n-56p-57n-54n-54p-55n --include_samples_by_name 54n-55n-56n-57n
  all_samp_tiles=$1
  # cat dmr_jobs.args | awk -F , -v all_samp_tiles=$all_samp_tiles 'match($0, /--samp_ids ([^ ]*)/, array)  {print $1",",  all_samp_tiles, $1, "--sample_names " array[1]}' > heatmap_jobs.args
  cat dmr_jobs.args | awk -F , -v all_samp_tiles=$all_samp_tiles 'match($0, /--samp_ids ([^ ]*)/, array)  match($0, /--output_dir ([^ ]*)/, array2) {print $1",",  all_samp_tiles, array2[1], "--sample_names " array[1]}' > heatmap_jobs.args
 #echo edit heatmap_jobs.args to set the path to all_samples_100bp_tiles.bed. And if you want to select part of the samples \(see /scripts_2022/make_heatmap.R\)
}


write_homer_jobs_args()
{
#format: sample_dir
 awk -F , 'match($0, /--output_dir ([^ ]*)/, array2){print array2[1]}' dmr_jobs.args  > homer_jobs.args
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
request_cpus = 10
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
