#!/usr/bin/env bash

#buiild the image
docker build -t  methylkit_dmrs:11.8.2024 .

#save the image
#docker save -o methylkit_dmrs_11.8.2024.tar methylkit_dmrs:11.8.2024

#delete from local
#docker rmi methylkit_dmrs:11.8.2024

#build singularity image
singularity build methylkit_dmrs_11_8_2024.img ./Singularity
