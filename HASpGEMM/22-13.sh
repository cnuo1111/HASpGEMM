input="../input_data/22.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		export OMP_NUM_THREADS=32
		export GOMP_CPU_AFFINITY="0-31"
		./haspgemm13 $str | tee -a ../plot/DATA/13-22.csv
	done
} < "$input"