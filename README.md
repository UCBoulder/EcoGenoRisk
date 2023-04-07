# EcoGenoRisk
EcoGenoRisk is a computational risk assessment tool designed specifically for assessing the risk a synthetic organism can introduce into the environment. EcoGenoRisk is composed of HazID, EnvCen, and RiskQ. 
HazID contains three modules: NicheOverlap, CompetitorFind, PoisInhibitor. HazID aims to assess the risk posed to genetically similar organisms. HazID operates on Python and R. All supporting documents and file are currently stored on OneDrive due to the size limitations on GitHub. 

## *NicheOverlap*: 
### Genomic and Enzymatic Library Creation   
**_1_configuring_library.py** : Main file which will call on **_2_genomic_data_download.py** and **_3_genomic_summary.py** <br />
**_2_genomic_data_download**: Creation of genomic and enzymatic libraries. This script will download genomic files for a domain of user's interest from NCBI RefSeq. The genomic files will be protein files of genomes that have complete assembly. DIAMOND sequence aligner will use the protein files to any present EC numbers. <br />
**_3_genomic_summary.py**: Creates the EC number binary summary matrix (EC BSM) for the sampled genomes. This program will run for each of the domains individually <br />
**_0_merged_binary_matrix_merged_lineage.py** : Combines Archaeal and Bacterial EC BSM to create one large file (combined_summary_matrix_YYYY_MM_DD.txt). This code will also compile the taxonomy datasheet from Genome Taxonomy Database to include full lineage for Archaea and Bacteria. Taxonomy files for each domain will be combined into one document. Taxonomy files will be save with the day the script was run (taxonomy_YYYY_MM_DD.txt). <br />

### Synbio Testing 
**_4_new_synbio_analysis**: This code will complete *NicheOverlap* analysis for a synthetic biology (synbio) organism. This will require a .faa file for the synbio organism. This code will run DIAMOND sequence aligner on the given .faa file, create a EC BSM for the organism, and send the results to **_5_synbio_distance_matrix.py**. <br /> 
**_5_synbio_distance_matrix.py**: Imports the synbio's EC BSM and returns a list of top match organisms by comparing the individual EC BSM. This script can also group by taxonomical preference and create a functional distance matrix that includes the synbio organism. <br />
**_7_pca_dendro.R** : This code will conduct EC number frequency analysis, PCA, Robinson-Foulds metric test, and will construct a phenotype-based dendrogram. <br /><br />

## *CompetitorFind*: 
**_1_competitorFind.py**: Completes top match (found in **_5_synbio_distance_matrix.py**) and synbio substrate overlap analysis. Finds the total substrate overlap between the two organism, substrate overlap in modified pathways, and mutualistic compounds that are found as products and reactants. <br />

### *CompetitorFind* Scoring Tests
**CompetitorFind_test_cases.py**: Creates 10 enzymatically similar genome pairs, 10 enzymatically different genome pairs, 1000 random genome pairs to assess. The script sends the genome pairs to **separatefunctions.py** for analysis, and later assigns a *CompetitorFind* score. <br />
**separatefunctions.py**: Contains the **_1_competitorFind.py** machinery to identify the substrate overlap. <br />
**Competitor_Find_Scores.R**: plots the data found in **CompetitorFind_test_cases.py**. <br /><br />
### CompetitorFind vs. Ecological Relationships
**kehe_dataset_testing.R**: Main script that launches Method 1 and Method 2 testing of ecological relationships using the Kehe et al. dataset. This script will subsets the Kehe et al. results dataset using a growth medium of interest  (sucrose, glutamine, etc.), and send the dataframe to **method1_kehe.R** for analysis.<br />
**method1_kehe.R**: Method 1 algorithm for evaluating the substrate overlap between genome groups outlined in Kehe et al. <br />
**CompetitorFind_Method1_plotting_results.R**: plots the output summary document for Method 1. Creates a dot plot for ecological relationships, completes a Wilcoxon test with Competition as the reference.<br />
## *PoisInhibitor* 
**_1_poisinhibitor.py**: Completes *PoisInhibitor* analysis for a user-inputted chemical agent. Compiles a list of organisms that are predicted to be inhibited by the chemical agent and its analogs. 
<br /><br />

**References**: <br />
Kehe, J., Ortiz, A., Kulesa, A., Gore, J., Blainey, P. C., & Friedman, J. (2021). Positive interactions are common among culturable bacteria. Science Advances, 7(45), eabi7159. https://doi.org/10.1126/sciadv.abi7159
