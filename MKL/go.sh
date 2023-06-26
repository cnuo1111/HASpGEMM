# ./wait
export OMP_NUM_THREADS=16
export GOMP_CPU_AFFINITY="0-15"
input="/home/codeis123/22.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
        ./mkl_spgemm $str | tee -a mkl_P_12_SAM.csv
	done
} < "$input"
export OMP_NUM_THREADS=8
export GOMP_CPU_AFFINITY="16-23"
input="/home/codeis123/22.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
        ./mkl_spgemm $str | tee -a mkl_E_12_SAM.csv
	done
} < "$input"


export OMP_NUM_THREADS=16
export GOMP_CPU_AFFINITY="0-15"
input="/home/codeis123/2717.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
        ./mkl_spgemm $str | tee -a mkl_P_12_ALL.csv
	done
} < "$input"
export OMP_NUM_THREADS=8
export GOMP_CPU_AFFINITY="16-23"
input="/home/codeis123/2717.csv"
{
	read
	i=1
	while IFS=',' read -r num str
	do
        ./mkl_spgemm $str | tee -a mkl_E_12_ALL.csv
	done
} < "$input"

