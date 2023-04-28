input="../input_data/22.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		./spa $str | tee -a ../plot/DATA/12-22-spa.csv
		./hash $str | tee -a ../plot/DATA/12-22-hash.csv
		./esc $str | tee -a ../plot/DATA/12-22-esc.csv
	done
} < "$input"