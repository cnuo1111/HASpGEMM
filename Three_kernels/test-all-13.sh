input="../input_data/all.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
		./spa $str | tee -a ../plot/DATA/13-22-spa.csv
		./hash $str | tee -a ../plot/DATA/13-22-hash.csv
		./esc $str | tee -a ../plot/DATA/13-22-esc.csv
	done
} < "$input"