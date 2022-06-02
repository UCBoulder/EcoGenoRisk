#Load the binary EC matrix 
library(stringr)
library(dplyr)
# bacteria_ec <- read.delim("bacteria_big_matrix_tab.txt", header = TRUE, quote = "",  sep= '\t')
# rownames(bacteria_ec) <- bacteria_ec$Name_of_Genome
# ec_numbers <- ncol(bacteria_ec)
# bacteria_ec_sum <- colSums(bacteria_ec%>%select_if(is.numeric))
# bacteria_ec_sum[bacteria_ec_sum>400]
# sum(bacteria_ec_sum/25899 > 0.75)
# sum(bacteria_ec_sum/25899 > 0.50)
# sum(bacteria_ec_sum/25899 > 0.25)

#hist(bacteria_ec_sum[bacteria_ec_sum>0]/25899)
# 
# # Read in the Processed Distance Matrix
# # Phylum
# d_phylum <-read.delim('DistanceListforSynbio_Phylum.txt',header=FALSE, sep = "\t", row.names = 1)
# names(d_phylum) <- row.names(d_phylum)
# # Add zeroes along the diag back in
# d_phylum[d_phylum<0.001] <- 0
# d_phylum_asdistance <- as.dist(d_phylum)
# fit_phylum <- hclust(d_phylum_asdistance, method="ward")
# pdf("20220601 Phylum Tree.pdf", height = 10, width = 10)
# plot(fit_phylum)
# dev.off()

bacteria_phylum <-read.delim("ec_space_phylum.txt", header = TRUE, quote = "",  sep= ' ')
ec_numbers <- ncol(bacteria_phylum)
bacteria_phylum_sum <- colSums(bacteria_phylum%>%select_if(is.numeric))
bacteria_phylum_sum[bacteria_phylum_sum>400]
sum(bacteria_phylum_sum/25899 > 0.75)
sum(bacteria_phylum_sum/25899 > 0.50)
sum(bacteria_phylum_sum/25899 > 0.25)
fit_prin_phylum <- princomp(t(bacteria_phylum[bacteria_phylum_sum>100]), cor=TRUE)
summary(fit_prin_phylum) # print variance accounted for
loadings(fit_prin_phylum) # pc loadings
plot(fit_prin_phylum,type="lines") # scree plot
write.table(fit_prin_phylum$scores,"phylum_scores.txt") # the principal components
biplot(fit_prin_phylum)

# # Read in the Processed Distance Matrix
# # Class
# d_class <-read.delim('DistanceListforSynbio_Class.txt',header=FALSE, sep = "\t", row.names = 1)
# names(d_class) <- row.names(d_class)
# # Add zeroes along the diag back in
# d_class[d_class<0.001] <- 0
# d_class_asdistance <- as.dist(d_class)
# fit_class<- hclust(d_class_asdistance, method="ward")
# pdf("20220601 Class Tree.pdf", height = 10, width = 10)
# plot(fit_class)
# dev.off()

bacteria_class <-read.delim("ec_space_class.txt", header = TRUE, quote = "",  sep= ' ')
ec_numbers <- ncol(bacteria_class)
bacteria_class_sum <- colSums(bacteria_class%>%select_if(is.numeric))
bacteria_class_sum[bacteria_class_sum>400]
sum(bacteria_class_sum/25899 > 0.75)
sum(bacteria_class_sum/25899 > 0.50)
sum(bacteria_class_sum/25899 > 0.25)
fit_prin_class <- princomp(t(bacteria_class[bacteria_class_sum>100]), cor=TRUE)
summary(fit_prin_class) # print variance accounted for
loadings(fit_prin_class) # pc loadings
plot(fit_prin_class,type="lines") # scree plot
write.table(fit_prin_class$scores,"class_scores.txt") # the principal components
biplot(fit_prin_class)
# 
# # Read in the Processed Distance Matrix
# # Order
# d_order<-read.delim('DistanceListforSynbio_Order.txt',header=FALSE, sep = "\t", row.names = 1)
# names(d_order) <- row.names(d_order)
# # Add zeroes along the diag back in
# d_order[d_order<0.001] <- 0
# d_order_asdistance <- as.dist(d_order)
# fit_order<- hclust(d_order_asdistance, method="ward")
# 
# pdf("20220601 Order Tree.pdf", height = 10, width = 10)
# plot(fit_order)
# dev.off()

bacteria_order <-read.delim("ec_space_order.txt", header = TRUE, quote = "",  sep= ' ')
ec_numbers <- ncol(bacteria_order)
bacteria_order_sum <- colSums(bacteria_order%>%select_if(is.numeric))
bacteria_order_sum[bacteria_order_sum>400]
sum(bacteria_order_sum/25899 > 0.75)
sum(bacteria_order_sum/25899 > 0.50)
sum(bacteria_order_sum/25899 > 0.25)
fit_prin_order <- princomp(t(bacteria_class[bacteria_order_sum>100]), cor=TRUE)
summary(fit_prin_order) # print variance accounted for
loadings(fit_prin_order) # pc loadings
plot(fit_prin_order,type="lines") # scree plot
write.table(fit_prin_order$scores,"order_scores.txt") # the principal components
biplot(fit_prin_order)

# # Read in the Processed Distance Matrix
# # Family
# d_family<-read.delim('DistanceListforSynbio_Order.txt',header=FALSE, sep = "\t", row.names = 1)
# names(d_family) <- row.names(d_family)
# # Add zeroes along the diag back in
# d_family[d_family<0.001] <- 0
# d_family_asdistance <- as.dist(d_family)
# fit_family<- hclust(d_family_asdistance, method="ward")
# 
# pdf("20220601 Family Tree.pdf", height = 10, width = 10)
# plot(fit_family)
# dev.off() 

bacteria_family <-read.delim("ec_space_family.txt", header = TRUE, quote = "",  sep= ' ')
ec_numbers <- ncol(bacteria_family)
bacteria_family_sum <- colSums(bacteria_family%>%select_if(is.numeric))
bacteria_family_sum[bacteria_family_sum>400]
sum(bacteria_family_sum/25899 > 0.75)
sum(bacteria_family_sum/25899 > 0.50)
sum(bacteria_family_sum/25899 > 0.25)
fit_prin_family <- princomp(t(bacteria_class[bacteria_family_sum>100]), cor=TRUE)
summary(fit_prin_family) # print variance accounted for
loadings(fit_prin_family) # pc loadings
plot(fit_prin_family,type="lines") # scree plot
write.table(fit_prin_family$scores,"family_scores.txt") # the principal components
biplot(fit_prin_family)


# 
# # Genera
# d_genera<-read.delim('DistanceListforSynbio.txt',header=FALSE, sep = "\t", row.names = 1)
# names(d_genera) <- row.names(d_genera)
# # Add zeroes along the diag back in
# d_genera[d_genera<0.001] <- 0
# d_genera_asdistance <- as.dist(d_genera)
# fit_genera <- hclust(d_genera_asdistance, method="ward")
# 
# pdf("20220601 Genera Tree.pdf", height = 20, width = 225)
# plot(fit_genera)
# dev.off()

bacteria_genus <-read.delim("ec_space_genus.txt", header = TRUE, quote = "",  sep= ' ')
ec_numbers <- ncol(bacteria_genus)
bacteria_genus_sum <- colSums(bacteria_genus%>%select_if(is.numeric))
bacteria_genus_sum[bacteria_genus_sum>400]
sum(bacteria_genus_sum/25899 > 0.75)
sum(bacteria_genus_sum/25899 > 0.50)
sum(bacteria_genus_sum/25899 > 0.25)
fit_prin_genus <- princomp(t(bacteria_class[bacteria_genus_sum>100]), cor=TRUE)
summary(fit_prin_genus) # print variance accounted for
loadings(fit_prin_genus) # pc loadings
plot(fit_prin_genus,type="lines") # scree plot
write.table(fit_prin_genus$scores,"family_scores.txt") # the principal components
biplot(fit_prin_genus)


#d <- dist(bacteria_ec[bacteria_ec_sum>1], method = "euclidean")
#fit <- hclust(d, method="ward")
#plot(fit)
#d2 <- dist(t(bacteria_ec[bacteria_ec_sums>100]), method = "euclidean")
#fit2 <- hclust(d2, method="ward")
#plot(fit2)

#fit_prin <- princomp(t(bacteria_ec[bacteria_ec_sums>100]), cor=TRUE)

#summary(fit_prin) # print variance accounted for
#loadings(fit_prin) # pc loadings
#plot(fit_prin,type="lines") # scree plot
#write.table(fit_prin$scores,"scores.txt") # the principal components
#biplot(fit_prin)

#help(biplot)