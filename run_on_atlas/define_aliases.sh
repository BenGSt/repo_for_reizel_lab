#!/bin/bash

set_aliases() {
  alias wgbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/prepare_jobs.sh'
  alias rrbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag/run.sh'
  alias dmrs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/run.sh'
  alias rna_seq='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/run.sh'
}

write_to_bashrc()
{
# Get the user's home directory
USER_HOME=$(eval echo ~${SUDO_USER:-$USER})
cat << EOF >>$USER_HOME/.bashrc

# Reizel Lab Pipeline Aliases
. /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/define_aliases.sh
EOF
echo -e "Custom aliases added to .bashrc\n"
# Reload the .bashrc file
source $USER_HOME/.bashrc
}

#main
if [ "$1" == "-write" ]; then
  source_aliases_in_bashrc
else
  set_aliases
fi


## Add aliases to the user's .bashrc file
#echo >>$USER_HOME/.bashrc
#echo "# Custom Aliases" >>$USER_HOME/.bashrc
#echo "alias wgbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/run_1_dag_per_sample.sh'" >>$USER_HOME/.bashrc
#echo "alias rrbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag/run.sh'" >>$USER_HOME/.bashrc
#echo "alias dmrs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/run.sh'" >>$USER_HOME/.bashrc
#echo "alias rna_seq='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/run.sh'" >>$USER_HOME/.bashrc


