container=../../containers/2024/dmrseq/dmrseq_28.11.2024.img
script=./find_dmrs_dmrseq.R

sampleheet=$(realpath ./example_samplesheet.csv)
output_dir=$(realpath ./example_output)

#get sample dirs from samplesheet to bind to container
cov_data_paths=$(dirname $(awk 'NR>1' "$sampleheet" | cut -d, -f4) | sort | uniq | sed "s|~|$HOME|")
cov_data_paths=$(echo $cov_data_paths | sed 's/ /,/g')

#write a new samplesheet with continer paths (replace 4th field dirnames with /coverage_data/)
cat $sampleheet | awk -F , '{
  OFS=",";
  if (NR==1) 
    print $0
  else {
    split($4, path, "/");
    $4="/coverage_data/"path[length(path)];
    print $0
  }
}' > ./container_samplesheet.csv

#run the script in the container
singularity exec --no-home \
  -B:./container_samplesheet.csv:/samplesheet.csv \
  -B:$output_dir:/output_dir \
  -B:$cov_data_paths:/coverage_data \
  $container Rscript $script \
    -ncpus 6 \
    -samplesheet /samplesheet.csv \
    -output_dir /output_dir \
    -min_qval 0.05 \
    cutoff 0.1

rm ./container_samplesheet.csv