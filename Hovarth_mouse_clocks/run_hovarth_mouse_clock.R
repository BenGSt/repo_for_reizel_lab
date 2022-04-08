#load CpG sites and model coefficients 
all_CpGs_elastic_net=read.csv("Thompson_2018_all_CpGs_elastic_net_coefficients_mm10.bed"
                              ,header = TRUE, sep = "\t")
all_CpGs_elastic_net_intersept = all_CpGs_elastic_net[1,2]
#remove first line contatining intersept
all_CpGs_elastic_net = all_CpGs_elastic_net[2:nrow(all_CpGs_elastic_net),]
