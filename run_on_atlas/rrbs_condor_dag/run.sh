
#!/bin/bash

# write condor_submition files
submission_file_names=(align_jobs.sub \
                make_tiles_jobs.sub \
                meth_calling_jobs.sub \

                trim_diversity_adaptors_jobs.sub \
                trim_illumina_adaptors_jobs.sub)
executables=()
for i in "${!submission_file_names[@]}"; do

  cat << EOF > ${submission_file_names[$i]}
executable = align_job.sh
Initialdir = /srv01/technion/bengst/storage/analyzed_data/KKTR-TargetingMafAMotifWithTet/test_rrbs_dag_2
request_cpus = 10
RequestMemory = 32GB
Arguments = \$(name) \$(args)
universe = vanilla
log = ./\$(name)/condor_logs/\$(name)_align.log
output = ./\$(name)/condor_logs/\$(name)_align.out
error = ./\$(name)/condor_logs/\$(name)_align.error
queue args from rrbs_jobs.args
EOF

done


ovation_rrbs_presubmit.sh