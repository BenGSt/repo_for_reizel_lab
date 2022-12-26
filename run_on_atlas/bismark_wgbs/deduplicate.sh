
main()
{
  source /Local/bfe_reizel/anaconda3/bin/activate ovation_rrbs_pipeline_2022
	script_name=$(echo $0 | awk -F / '{print $NF}')

	echo
	echo
	echo \#################################
	echo \#################################
	echo running: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

  time bismark_deduplicate "$@"

	echo
	echo
	echo \#################################
	echo \#################################
	echo finished: $script_name "$@"
	echo date: $(date)
	echo \#################################
	echo \#################################
	echo
	echo
}

