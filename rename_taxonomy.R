### Load in the files that are needed ###
# Need the actually output file
bacteria_binary   <- read.table("bacteria_combined1.csv", header=TRUE, quote = "", sep = ' ')
# Need the assembly summary describing the entire database that was downloaded and other ID values
assembly_summary  <- read.delim("bacteria_assembly_summary.txt",header=TRUE, quote = "", sep = '\t')  #reads as a table, separated by tabs
# Need the lineage names from the taxonomy database
full_lineage      <- read.csv("fullnamelineage.dmp", header=TRUE, quote = "", sep='\t')

# Note, these are imported in as dataframes, using the $ symbol allows you to access each named column
# Take a look at these two
full_lineage$X1
full_lineage$X

# These have funny names because they were imported in oddly, to cleanup things and using only the columns we need
# create a new dataframe and rename the columns with useful descriptors
lineage_only           <- as.data.frame(cbind(full_lineage$X1,full_lineage$X))
colnames(lineage_only) <- c("taxid","id_description")
nrow(lineage_only)

# Now, one thing that is unfortunate is that the initial larger bacteria matrix preserved names that are not
# directly comparable, we just need to process the name column a bit to allow direct merging of rows. First
# we split based on _. Then, we add "GDF_" to the beginning of the ID # to make it directly compatable, finally
# we merge these together and repalce the row column names.
names_of_interest              <- do.call(rbind,strsplit(as.character(bacteria_binary$Name_of_Genome),"_"))
wanted_reference               <- paste("GCF_",names_of_interest[,2],sep='')
bacteria_binary$Name_of_Genome <- wanted_reference


# For the assembly summary, to cleanup things and using only the columns we need
# create a new dataframe and rename the columns with useful descriptors
nrow(assembly_summary)
assembly_only <- data.frame(assembly_summary$assembly_accession, assembly_summary$taxid)
colnames(assembly_only) <- c("Name_of_Genome","taxid")
colnames(lineage_only) <-c("taxid")

# Now we just merge!
assembly_taxid           <- merge(assembly_only, full_lineage, by="taxid")
nrow(assembly_taxid)
nrow(bacteria_binary)
bacteria_binary_taxonomy <- merge(assembly_taxid,bacteria_binary,by="Name_of_Genome",all.y = TRUE)
nrow(bacteria_binary_taxonomy)
sub=subset(bacteria_binary_taxonomy, select=X1.1.1.1:X7.6.2.16)
Output_names_only           <- data.frame(bacteria_binary_taxonomy$root,sub)

