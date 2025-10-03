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

# Open the CSV file with tab delimiter
with open(fname, 'r') as file:
    # Create a CSV reader with tab delimiter
    csv_reader = csv.reader(file, delimiter='\t')
    
    # Iterate over each row, in our case over each feature
    for row in csv_reader:
        # Skip header rows (lines starting with #) or empty rows
        if len(row) == 0 or (len(row) > 0 and row[0].startswith('#')):
            continue
            
        # Check if row has at least 3 columns
        if len(row) < 3:
            print(f"Warning: Skipping malformed row with {len(row)} columns: {row}")
            continue
            
        query = row[0]
        cog = row[1]
        
        # Handle GO terms - everything from column 2 onwards
        go_array = row[2:]
        go = ','.join(go_array)
        
        divide = [query, cog, go]
        divide = pd.DataFrame([divide], columns=column_names)
        df = pd.concat([df, divide], ignore_index=True)

df.to_csv(fname, index=False)
print(f"Processed {len(df)} rows successfully")