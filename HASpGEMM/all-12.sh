input="../input_data/all.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		export OMP_NUM_THREADS=24
		export GOMP_CPU_AFFINITY="0-23"
		./haspgemm12 $str | tee -a ../plot/DATA/12-all.csv
	done
} < "$input"