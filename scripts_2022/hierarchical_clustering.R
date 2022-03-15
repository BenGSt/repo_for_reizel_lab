all_samples_100bp_tiles_no_NA <- read.delim("~/extra_storage_1TB/analyzed_data/KKTR-FahRenegeration-rrbs/hierarchical_clustering/all_samples_100bp_tiles_no_NA.bed" )
df=t(all_samples_100bp_tiles_no_NA)

# Original rownames
# "chrom"  "start"  "end"    "AAGAGG" "AAGCCT" "ACCTCA" "AGCATG" "AGTGAG" "GAGTCA" "GCACTA" "GGAGAA" "GTCGTA" "GTGCTT" "TGGTGA"
rownames(df)=c("chrom", "start", "end", "YoungYoung", "YoungYoung", "OldOld", "YoungYoungProlong", "Young", "YoungYoungProlong", "OldOld", "Old", "Young", "OldOld", "Old")
df=df[!(rownames(df) %in% c("YoungYoungProlong")),]

distance_matrix=dist(df[4:nrow(df),])

hc=hclust(distance_matrix,method = "centroid")


png("hc_all_samples_100bp_tiles_no_NA.jpg")
plot(hc)
dev.off()
