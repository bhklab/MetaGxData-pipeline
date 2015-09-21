module load R

for subtype in All Basal Her2 LumB LumA
do
	qsub -cwd -b y -e sge_out -o sge_out -q bhklab -N "SAPS_$subtype" -t 1-585 "module load R; Rscript saps_parallel_brca.R $subtype"
done
