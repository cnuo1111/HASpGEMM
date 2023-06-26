input="../input_data/all.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		./mkl_spgemm $str | tee -a ../plot/DATA/mkl-all.csv
	done
} < "$input"