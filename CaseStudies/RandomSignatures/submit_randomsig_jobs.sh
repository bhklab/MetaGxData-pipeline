qsub cwd -b y -e sge_out -o sge_out -N "randomsigs" -t 5-200:5 "module load R; Rscript randomSignatures.R"
