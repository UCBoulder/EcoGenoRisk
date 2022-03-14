a1<-read.delim('archaea_genome_to_genome_complete_matrix.csv',header=TRUE)

# Summary Stats
library(tidyr)
new1<-gather(a1, key="Genome", value="Score", 2:408)
print(new1)
summary(new1)
d<-density(new1$Score)
plot(d, main="Distribution of Difference Based Scoring")

#NMDS
library(vegan)
#make community matrix
com= a1[,2:ncol(a1)]
m_com=as.matrix(com)
set.seed(123)
nmds=metaMDS(m_com, distance="bray", k=2, trymax=100)
plot(nmds)
data.scores=as.data.frame(scores(nmds))
data.scores$Sample=a1$Sample
data.scores$Time=a1$Time
