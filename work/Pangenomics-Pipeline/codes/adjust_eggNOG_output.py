import csv
import sys
import pandas as pd

# Define the column names
column_names = ['Query', 'COG Category', 'GO terms']

# Create the DataFrame with column names
df = pd.DataFrame(columns=column_names)

fname = sys.argv[1]

# Read the CSV file and exclude the last three lines
with open(fname, 'r') as file:
    lines = file.readlines()[:-3]

# Write the modified content back to the file
with open(fname, 'w') as file:
    file.writelines(lines)

# Open the CSV file
with open(fname, 'r') as file:
	# Create a CSV reader
	csv_reader = csv.reader(file)

	# Iterate over each row, in our case over each feature
	for row in csv_reader:
		query = row[0].split("\t")[0]
		cog = row[0].split("\t")[1]

		go_array = row[1:]
		go_array.insert(0, row[0].split("\t")[2])
		go = ','.join(go_array)
		
		divide = [query, cog, go]
		divide = pd.DataFrame([divide], columns=column_names)
		df = df.append(divide, ignore_index=True)
df.to_csv(fname, index=False)