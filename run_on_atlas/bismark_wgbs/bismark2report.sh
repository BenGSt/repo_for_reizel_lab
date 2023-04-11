
main()
{
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

  source /Local/bfe_reizel/anaconda3/bin/activate wgbs_bismark_pipeline_2023
  arg_parse "$@"
	cd "$output_dir" || exit 1


  time bismark2report --splitting_report *splitting_report.txt --mbias_report  *M-bias.txt

  echo
	echo
	echo \#################################
	echo \#################################
	echo finished: $script_name "$@"
	echo date: $(date)
	echo hostname: $(hostname)
	echo pwd: $(pwd)
	echo \#################################
	echo \#################################
	echo
	echo

}

arg_parse()
{
  while [[ $# -gt 0 ]]; do
    case $1 in
      -output-dir)
        output_dir="$2"
        shift
        shift
        ;;
        *)
        echo error: only valid arg is  -output-dir \<path\>
        exit 1
        ;;
    esac
done
}


main "$@"