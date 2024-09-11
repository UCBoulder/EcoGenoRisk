Neighbor Replace
================
Checks for likely occupancy of a synbio organism within a given habitat through a taxonomic similarity search. 

What this script does is takes in a synbio 16s rRNA FASTA and a tsv file from JGI (bins by ecosystem, select all, export) and extracts the first 
16s sequence, then blasts it. Turn on and off the BLAST definition if you are running a new analysis or just working with the output. 

In other words, the blast definition can be toggled to run a blast or not run a blast. 

The blast xml gets outputted to a file, which is then read back in and parsed. Then, it's just data processing to compare the species that come out of the blast to the species 
that are present in the biome. All_remaining is the final dataframe. 

