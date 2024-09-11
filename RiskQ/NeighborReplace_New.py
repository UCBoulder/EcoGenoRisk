###September 5th, 2024
##This script is a taxnomic similairty search between the input synbio organism and all of the organisms present in the biome 
#This script utilizes the new blast on biopython. https://biopython.org/docs/latest/Tutorial/chapter_blast.html
import Bio
from Bio import Blast
from Bio import SeqIO
import pandas as pd 

Blast.email = 'john.docter@colorado.edu'
#________________________________________________________________________________________________________________________________________________________________________________
##1) Input Synbio FASTA
def input_synbio(input_fasta): 
    final = []
    sixteen_s = SeqIO.parse(input_fasta, 'fasta') #https://www.ncbi.nlm.nih.gov/search/all/?term=Variovorax%20paradoxus%2016s%20rRNA
    for NA in sixteen_s:
        sequence = str(NA.seq)
        final.append(sequence)
        # print(final[0]) #This is a list of all the sequences 

    return final


#2) Import Biome FASTA and parse for genus and species 
def biome_input(input_biome_tsv):
    biome = pd.read_csv(input_biome_tsv, sep='\t')
    lineage = biome['Bin Lineage'] #We want just the genus and species as that is the BLAST output 
    genusspecies = biome['Bin Lineage'].str.rsplit(';', n=2).str[-1:].str.join('')
    genusspecies = genusspecies.reset_index(drop=True)
    genusspecies = genusspecies.to_frame()
    genusspecies.columns = ['tax']
    genusspecies['tax'] = genusspecies['tax'].str.strip()
    print('genusspecies')
    print(genusspecies)

    return genusspecies



##3) BLASTn 
def Bloost(rRNA): 
    sequence1=rRNA[0]
    print('Blasting.....')
    result_stream = Blast.qblast('blastn', 'refseq_rna', sequence1)

    with open("blast_output.xml", "wb") as out_stream: #writes the contents to a file and then reopens to so if there's an issue, don't have to rerun BLAST all over again
        out_stream.write(result_stream.read())
    result_stream.close()

def species_present(): 
    result_stream = open("blast_output.xml", "rb")
    blast_record = Blast.read(result_stream)

    # print(len(blast_record))

    all_records = len(blast_record)
    all_blast_species = []
    for i in range(all_records):
        hit = blast_record[i]
        text = hit.target.description
        # print(text.split())
        text= text.split()[:2]
        text = ' '.join(text)
        all_blast_species.append(text)
    return all_blast_species

        # all_blast_species.append(hit.target.description)

#5) Merge on scientific name and write results to an external file. 
#Turn all blast species into dataframe 
def merge(all_blast_species, genusspecies): 
    all_blast_species_df = pd.DataFrame(all_blast_species, columns=['tax'])
    all_blast_species_df = all_blast_species_df.drop_duplicates()
    all_blast_species_df['tax'] = all_blast_species_df['tax'].str.strip()
    # print(all_blast_species_df)
    remaining = pd.merge(genusspecies, all_blast_species_df,  how='inner', on='tax')
    remaining = remaining.rename(columns={'tax': 'genus species remaining'})
    return remaining, all_blast_species_df

#This accounts for the 
def species_check(all_blast_species_df): 
    all_blast_species_df['tax'] = all_blast_species_df['tax'].str.split().str.get(0)

    remaining2 = pd.merge(biome, all_blast_species_df, how='inner', on=['tax'])
    remaining2 = remaining2.drop_duplicates()
    remaining2 = remaining2.rename(columns={'tax': 'species remaining'})
 
    return remaining2




##Calling Script______________________________________________________________________________________________________
input_fasta = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/RiskQ/Validation Run/Variovorax paradoxus 16s rRNA.fasta'
rRNA = input_synbio(input_fasta) #creates rRNA Sequence to BLAST
input_biome = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/RiskQ/Validation Run/Lake_Biome_Data.tsv'
biome = biome_input(input_biome) #creates list of genus species from the JGI Biome file 
# Bloost(rRNA) ###Turn off here above if you've already blasted. Blasts the input fasta rRNA sequence and writes the results to an xml file 
all_blast_species = species_present() #reads in the xml and parses for all species present in BLAST Output 
# print(all_blast_species)
remaining, all_blast_species_df = merge(all_blast_species, biome)

remaining2 = species_check(all_blast_species_df)
all_remaining = pd.concat([remaining, remaining2])
print(all_remaining)

##A message for future john when writing the RiskQ paper 
#What this script does is takes in a 16s rRNA file from NCBI and a tsv file from JGI (bins by ecosystem, select all, export) and etracts the first 
#16s sequence, then blasts it. Turn on and off the BLAST definition if you are running a new analysis or just working with the output. In other words, the blast definition can be toggled to run 
# a blast or not run a blast. the blast xml gets outputted to a file, which is then read back in and parsed. Then, it's just data processing to compare the species that come out of the blast to the species 
#that are present in the biome. All_remaining is the final dataframe. 