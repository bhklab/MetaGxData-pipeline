# Generate pooled.eset.intersecting.genes.RData, to be used in child jobs
Rscript saps_parallel_init.R

for subtype in All Basal Her2 LumB LumA
do
	qsub cwd -b y -e sge_out -o sge_out -N "SAPS_$subtype" -t 1-140 "Rscript saps_parallel.R $subtype"
done
