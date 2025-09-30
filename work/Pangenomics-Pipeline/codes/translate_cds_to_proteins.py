from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

# Usage: python3 translate_cds_to_proteins.py input_cds.fasta output_proteins.faa

input_file = sys.argv[1]
output_file = sys.argv[2]

protein_records = []

for record in SeqIO.parse(input_file, "fasta"):
    protein_seq = record.seq.translate(to_stop=True)
    protein_record = SeqRecord(protein_seq, id=record.id, description="translated")
    protein_records.append(protein_record)

SeqIO.write(protein_records, output_file, "fasta")
