import csv
matrix_to_max_perf1 = {} 
with open('./DATA/13-spa.csv', 'r') as file1:
    reader = csv.reader(file1)
    for row in reader:
        matrix_name = row[1]
        max_perf = float(row[6])
        if matrix_name in matrix_to_max_perf1:
            matrix_to_max_perf1[matrix_name] = max(matrix_to_max_perf1[matrix_name], max_perf)
        else:
            matrix_to_max_perf1[matrix_name] = max_perf

with open('./DATA/13-spa.csv', 'r') as file2:
    reader = csv.reader(file2)
    rows = []
    for row in reader:
        matrix_name = row[1]
        if matrix_name in matrix_to_max_perf1:
            max_perf1 = matrix_to_max_perf1[matrix_name]
            max_perf2 = float(row[6])
            nnz = row[2]
            cub=row[4]
            if(max_perf1 < 0.1 or max_perf2 < 0.1):
                continue
            speedup = max_perf2 / max_perf1
            rows.append([matrix_name,nnz,cub, max_perf1, max_perf2, speedup])

with open('./DATA/result_SPA.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Matrix Name','nnz','cub', 'Max Perf (Method 1)', 'Max Perf (Method 2)', 'Speedup'])
    writer.writerows(rows)


matrix_to_max_perf1 = {}  
with open('./DATA/13-hash.csv', 'r') as file1:
    reader = csv.reader(file1)
    for row in reader:
        matrix_name = row[1]
        max_perf = float(row[6])
        if matrix_name in matrix_to_max_perf1:
            matrix_to_max_perf1[matrix_name] = max(matrix_to_max_perf1[matrix_name], max_perf)
        else:
            matrix_to_max_perf1[matrix_name] = max_perf

with open('./DATA/13-all.csv', 'r') as file2:
    reader = csv.reader(file2)
    rows = []
    for row in reader:
        matrix_name = row[1]
        if matrix_name in matrix_to_max_perf1:
            max_perf1 = matrix_to_max_perf1[matrix_name]
            max_perf2 = float(row[6])
            nnz = row[2]
            cub=row[4]
            if(max_perf1 < 0.1 or max_perf2 < 0.1):
                continue
            speedup = max_perf2 / max_perf1
            rows.append([matrix_name,nnz,cub, max_perf1, max_perf2, speedup])

with open('./DATA/result_HASH.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Matrix Name','nnz','cub', 'Max Perf (Method 1)', 'Max Perf (Method 2)', 'Speedup'])
    writer.writerows(rows)

matrix_to_max_perf1 = {}  
with open('./DATA/13-esc.csv', 'r') as file1:
    reader = csv.reader(file1)
    for row in reader:
        matrix_name = row[1]
        max_perf = float(row[6])
        if matrix_name in matrix_to_max_perf1:
            matrix_to_max_perf1[matrix_name] = max(matrix_to_max_perf1[matrix_name], max_perf)
        else:
            matrix_to_max_perf1[matrix_name] = max_perf

with open('./DATA/13-all.csv', 'r') as file2:
    reader = csv.reader(file2)
    rows = []
    for row in reader:
        matrix_name = row[1]
        if matrix_name in matrix_to_max_perf1:
            max_perf1 = matrix_to_max_perf1[matrix_name]
            max_perf2 = float(row[6])
            nnz = row[2]
            cub=row[4]
            if(max_perf1 < 0.1 or max_perf2 < 0.1):
                continue
            speedup = max_perf2 / max_perf1
            rows.append([matrix_name,nnz,cub, max_perf1, max_perf2, speedup])

with open('./DATA/result_ESC.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Matrix Name','nnz','cub', 'Max Perf (Method 1)', 'Max Perf (Method 2)', 'Speedup'])
    writer.writerows(rows)

matrix_to_max_perf1 = {}  
with open('./DATA/13-mkl.csv', 'r') as file1:
    reader = csv.reader(file1)
    for row in reader:
        matrix_name = row[1]
        max_perf = float(row[6])
        if matrix_name in matrix_to_max_perf1:
            matrix_to_max_perf1[matrix_name] = max(matrix_to_max_perf1[matrix_name], max_perf)
        else:
            matrix_to_max_perf1[matrix_name] = max_perf

with open('./DATA/13-all.csv', 'r') as file2:
    reader = csv.reader(file2)
    rows = []
    for row in reader:
        matrix_name = row[1]
        if matrix_name in matrix_to_max_perf1:
            max_perf1 = matrix_to_max_perf1[matrix_name]
            max_perf2 = float(row[6])
            nnz = row[2]
            cub=row[4]
            if(max_perf1 < 0.1 or max_perf2 < 0.1):
                continue
            speedup = max_perf2 / max_perf1
            rows.append([matrix_name,nnz,cub, max_perf1, max_perf2, speedup])

with open('./DATA/result_mkl.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Matrix Name','nnz','cub', 'Max Perf (Method 1)', 'Max Perf (Method 2)', 'Speedup'])
    writer.writerows(rows)