import pandas as pd

# Create an array to hold the data for each file
data = []

# Iterate over each file, read the Method, matrix, avg, and max columns, and add them to data
for i in range(1, 6):
    filename = f"13-22-{i}.csv"
    df = pd.read_csv(filename, usecols=[0, 1, 5, 6])
    # df.columns = [f"col{i}_1", f"col{i}_2", f"col{i}_6", f"col{i}_7"]
    data.append(df)

# Merging data
merged_data = pd.concat(data, axis=1)

merged_data.to_csv("./DATA/output.csv", index=False)
