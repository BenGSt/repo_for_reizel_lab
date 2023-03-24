#!/bin/bash

#
#alias wgbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/run_1_dag_per_sample.sh'
#alias rrbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag/run.sh'
#alias dmrs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/run.sh'
#alias rna_seq='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/run.sh'

# Get the user's home directory
USER_HOME=$(eval echo ~${SUDO_USER:-$USER})

# Add aliases to the user's .bashrc file
echo "" >> $USER_HOME/.bashrc
echo "# Custom Aliases" >> $USER_HOME/.bashrc
echo "alias wgbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/bismark_wgbs/run_1_dag_per_sample.sh'" >> $USER_HOME/.bashrc
echo "alias rrbs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rrbs_condor_dag/run.sh'" >> $USER_HOME/.bashrc
echo "alias alias dmrs='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/dmrs_condor_dag/run.sh'" >> $USER_HOME/.bashrc
echo "alias rna_seq='/storage/bfe_reizel/bengst/repo_for_reizel_lab/run_on_atlas/rna_hisat2_htseq_DESeq2_condor_dag/run.sh'" >> $USER_HOME/.bashrc

# Reload the .bashrc file
source $USER_HOME/.bashrc

# Let the user know that the aliases have been added
echo "Custom aliases added to .bashrc!"
