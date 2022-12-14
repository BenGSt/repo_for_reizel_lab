#!/usr/bin/bash

help()
{
  cat << EOF > echo
  Set up sra download jobs. Run from the dir you want the files to be downloaded to
  this scripts uses prefetch and fastq-dump from the sratoolkit. fastqdump is being deprecated So this script needs to be updated to fasterq-dump

EOF
}
  cat << EOF > sra_download.sh
#!/bin/bash
prefetch -p \$1
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip \$1
EOF
  chmod +x sra_download.sh

  cat << EOF > sra_download.sub
getenv = True
executable = ./sra_download.sh
Initialdir = $pwd
request_cpus = 1
RequestMemory = 500MB
Arguments = \$(srr)
universe = vanilla
log = sra_download_\$(srr).log
output = sra_download_\$(srr).out
error = sra_download_\$(srr).out
queue srr from (
$(awk -F , 'NR > 1 {print $1}' SraRunTable.txt)
)
EOF

    echo condor_submit sra_download.sub