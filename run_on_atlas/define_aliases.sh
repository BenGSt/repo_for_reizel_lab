#!/bin/bash

# USAGE: set aliases in current shell: source define_aliases.sh
#        make .bashrc set the aliases: define_aliases.sh -write

write_to_bashrc()
{
# Get the user's home directory
USER_HOME=$(eval echo ~${SUDO_USER:-$USER})
cat << EOF >>$USER_HOME/.bashrc

# Reizel Lab Pipeline Aliases
. /storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/define_aliases.sh
EOF
echo -e "Custom aliases added to .bashrc\n"
}

#main
if [ "$1" == "-write" ]; then
  write_to_bashrc
else
  alias wgbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/prepare_jobs.sh'
  alias rrbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag/run.sh'
  alias dmrs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/run.sh'
  alias rna_seq='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/run.sh'
fi
