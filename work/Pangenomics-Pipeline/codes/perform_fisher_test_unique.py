import sys
import csv
import os
import math
import pandas as pd
from scipy.stats import fisher_exact

# File containing query, COG Category, and GO terms
fname = sys.argv[1]

# CD-HIT Sequence output file containing representative sequences 
representatives_file = sys.argv[2]

# File containing info about gene classes (core, unique, accessory)
pangenome_file = sys.argv[3]

# Species Name for reference 
species = sys.argv[4]

# Total number of genes
numb_genes = int(sys.argv[5])

df = pd.read_csv(fname)
#df = df.drop(0)

# First value is number of unique COG X genes, second is number of non-unique COG X genes
COG_counts = {
    "A": [0,0],
    "B": [0,0],
    "C": [0,0],
    "D": [0,0],
    "E": [0,0],
    "F": [0,0],
    "G": [0,0],
    "H": [0,0],
    "I": [0,0],
    "J": [0,0],
    "K": [0,0],
    "L": [0,0],
    "M": [0,0],
    "N": [0,0],
    "O": [0,0],
    "P": [0,0],
    "Q": [0,0],
    "R": [0,0],
    "S": [0,0],
    "T": [0,0],
    "U": [0,0],
    "V": [0,0],
    "W": [0,0],
    # "X": [0,0],
    # "Y": [0,0],
    "Z": [0,0]
}

with open(pangenome_file, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    gene_class = []
    unique_count = 0
    for row in csvreader:
        gene_class.append(row[3])
        if row[3] == 'unique':
            unique_count += 1
csvfile.close()

def get_id (line):
    identifier = line.split(" ")[0]
    return identifier.split(">")[1]

for index, row in df.iterrows():
    category = row['COG Category']
    patric_id = row['Query']
    representatives = open(representatives_file,'r')
    rep_lines = representatives.readlines()

    cluster_count = -1

    if category == '-':
        continue

    else:
        for line in rep_lines: 
            if line.startswith('>'):
                cluster_count += 1
                seq_id = get_id(line)
                if seq_id == patric_id:
                    if gene_class[cluster_count] == 'unique':
                        if len(category) > 1:
                            individual_category = [char for char in category]
                            for char in individual_category:
                                COG_counts[char][0]+= 1
                        else:
                            COG_counts[category][0] += 1
                    else:
                        if len(category) > 1:
                            individual_category = [char for char in category]
                            for char in individual_category:
                                COG_counts[char][1] += 1
                        else:
                            COG_counts[category][1] += 1
            else:
                continue

odds_ratio = {}
p_values = {}

for key, value in COG_counts.items():
    column_names = ['Unique Genes', 'Non-unique Genes']
    row_names = ['With COG '+key, 'Without COG '+key]

    # A is the number of core genes with COG J
    A = value[0]
    # B is the number of non-core genes with COG J
    B = value[1]
    # C is the number of core genes without COG J
    C = unique_count - A
    # D is the number of non-core genes without COG J
    D = (numb_genes - unique_count) - B

    contingency_table = [[A, B], [C, D]]

    # Apply Fisher's exact test
    odds_ratio[key], p_values[key] = fisher_exact(contingency_table)

keys_to_delete = []
for key, value in odds_ratio.items():
    if math.isnan(value):
        keys_to_delete.append(key)

for key in keys_to_delete:
    del odds_ratio[key] 

df_odds_ratio = pd.DataFrame.from_dict(odds_ratio, orient='index', columns=[species])
df_p_values = pd.DataFrame.from_dict(p_values, orient='index', columns=[species])

df_odds_ratio.to_csv(species+"/"+species+"_Unique_COG_LORs.csv", index=True)
df_p_values.to_csv(species+"/"+species+"_Unique_COG_pvals.csv", index=True)
    


            

    
        

    
