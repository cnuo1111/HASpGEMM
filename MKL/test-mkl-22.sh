input="../input_data/22.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		./mkl_spgemm $str | tee -a ../plot/DATA/mkl-22.csv
	done
} < "$input"