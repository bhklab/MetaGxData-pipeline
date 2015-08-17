module load R

for subtype in All nonAngiogenic Angiogenic
do
	qsub cwd -b y -e sge_out -o sge_out -N "SAPS_$subtype" -t 1-443 "Rscript saps_parallel_ovca.R $subtype"
done
