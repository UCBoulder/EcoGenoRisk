###September 5th, 2024
##This script is a taxnomic similairty search between the input synbio organism and all of the organisms present in the biome 
import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

import pandas as pd 


#________________________________________________________________________________________________________________________________________________________________________________
##1) Input Synbio FASTA
final = []
sixteen_s = SeqIO.parse('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/RiskQ/Vibrio_Natriegens_16S.fasta', 'fasta')
for NA in sixteen_s:
     sequence = str(NA.seq)
     final.append(sequence)
print(final[0]) #This is a list of all the sequences 


#2) Import Biome FASTA and parse for genus and species 
biome = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/RiskQ/biome_data.tsv', sep='\t')
lineage = biome['Bin Lineage'] #We want just the genus and species as that is the BLAST output 
genusspecies = biome['Bin Lineage'].str.rsplit(';', n=2).str[-2:].str.join('')
# print(genusspecies)


##3) BLASTn 
sequence1=final[0]
Bout1 = NCBIWWW.qblast('blastn', 'refseq_rna', sequence1)
Bout2 = NCBIXML.parse(Bout1)

for record in Bout2: 
     query_length = record.query_length 


#4) Parse Blast Output for relevant information 




#5) Merge on scientific name and write results to an external file. 




