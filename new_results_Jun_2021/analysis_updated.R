setwd("~/Desktop/metabarcode_cerrado/barcode_cerrado/new_results_Jun_2021")
# rm(list=ls())

library(phyloseq)
library("data.table")
library(bipartite)
library(ggplot2)
#library(gplots)
library(vegan)
library("devtools")
devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)

#setwd("/Volumes/MBD-Storage/Project_global/data_ITS2/ITS2_Brazil_Meliponini.import")

tax <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T,sep=",")))
otu <- otu_table(read.table("asv_table.merge.txt"), taxa_are_rows=T)
map <- sample_data(read.table("sample_list.csv", sep=";", header=T,row.names=1))

sum(!(sample_names(map) %in% sample_names(otu)))
sum(!(sample_names(otu) %in% sample_names(map)))

# combine datasets
(data <- merge_phyloseq(otu,tax,map))
dataset.comp <- data

# renaming unresolved taxa
tax_table(dataset.comp)[tax_table(dataset.comp)[,"phylum"]=="","phylum"]<-paste(tax_table(dataset.comp)[tax_table(dataset.comp)[,"phylum"]=="","kingdom"],"_spc",sep="")
tax_table(dataset.comp)[tax_table(dataset.comp)[,"order"]=="","order"]<-paste(tax_table(dataset.comp)[tax_table(dataset.comp)[,"order"]=="","phylum"],"_spc",sep="")
tax_table(dataset.comp)[tax_table(dataset.comp)[,"family"]=="","family"]<-paste(tax_table(dataset.comp)[tax_table(dataset.comp)[,"family"]=="","order"],"_spc",sep="")
tax_table(dataset.comp)[tax_table(dataset.comp)[,"genus"]=="","genus"]<-paste(tax_table(dataset.comp)[tax_table(dataset.comp)[,"genus"]=="","family"],"_spc",sep="")
tax_table(dataset.comp)[tax_table(dataset.comp)[,"species"]=="","species"]<-paste(tax_table(dataset.comp)[tax_table(dataset.comp)[,"species"]=="","genus"],"_spc",sep="")

#Remove chlorophytes and off target products
dataset.comp.filter = dataset.comp
(dataset.comp.filter = subset_taxa(dataset.comp.filter, phylum=="p:Streptophyta"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, family!="p:Streptophyta_spc_spc"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, species!="k:Viridiplantae_spc_spc_spc_spc_spc"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, species!="p:Streptophyta_spc_spc_spc_spc"))

#Make taxa labels nice for plots
replace_tax_prefixes <- function(phyloseq){
  tmp_tax_table <- apply(tax_table(dataset.comp.filter), c(1, 2),function(y) gsub("^\\w:","",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_spc_.*","_spc",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_"," ",y))
  tax_table(phyloseq)<- tmp_tax_table
  return(phyloseq)
}

tail(tax_table(dataset.comp.filter))

dataset.comp.filter <- replace_tax_prefixes(dataset.comp.filter)

# collating zotus to species
(data.species <- tax_glom(dataset.comp.filter,taxrank="species"))
taxa_names(data.species) <- tax_table(data.species)[,"species"]

tail(tax_table(data.species))



# check throughput
sum(sample_sums(data.species) < 1000)
sum(sample_sums(data.species) < 500)
# it seems 24 samples have low throughput (<1000), 23 of those very low, I suggest to remove them, these are:

sample_data(data.species)[sample_sums(data.species) < 1000,]
# as expected, a lot of the honey, but also some pollen samples

# good ones are: 
sample_data(data.species)[sample_sums(data.species) >= 1000,]

data.species.bak<-data.species
#data.species <- prune_samples(sample_sums(data.species) >= 1000, data.species) # here they are removed. if you still want them included, just change the >= to 0 and then everything should be kept. 
# i would say, for networks and community analyses, e.g. diversity etc, samples should have at least 500-1000 reads. Sometimes, if it is known to be a specialist, can reduce to 200-300. 
# for all below that, you can look at them, to just see what was found inside, but it will probably not be close to complete foraging spectra and abundances might be very problematic

# transformation to relative data
data.species.rel = transform_sample_counts(data.species, function(x) x/sum(x))

# see count distribution of taxa before filtering
par(mar=c(4,15,1,1))
barplot(t(as.data.frame(sort(taxa_sums(data.species)[1:50], decreasing=T))), las=2, horiz = T)

# relative data and removal of rare taxa
otu_table(data.species.rel)[otu_table(data.species.rel)<0.01 ]<-0
otu_table(data.species)[otu_table(data.species.rel)<0.01 ]<-0
data.filter		= prune_taxa(taxa_sums(data.species)>0, data.species)
(data.rel.filter = prune_taxa(taxa_sums(data.species.rel)>0, data.species.rel))

# see count distribution of taxa after low-abundance filtering
par(mar=c(4,15,1,1))
barplot(t(as.data.frame(sort(taxa_sums(data.rel.filter), decreasing=T)[1:50])), las=2, horiz = T)
write.table(as.data.frame(sort(taxa_sums(data.rel.filter), decreasing=T)[1:100]), file="most_abundant_taxa2.csv",sep=";")


#splitting dataset
data.meli.rel <- subset_samples(data.rel.filter, include=="yes") 
data.other.rel <- subset_samples(data.rel.filter, include=="no") 
data.meli <- subset_samples(data.filter, include=="yes") 
data.other <- subset_samples(data.filter, include=="no") 


data.rel.filter.all <- data.rel.filter
data.filter.all <- data.filter

data.rel.filter <-data.meli.rel
data.filter <-data.meli

(data.rel.filter = prune_taxa(taxa_sums(data.rel.filter)>0, data.rel.filter))
taxa_sums(subset_taxa(data.rel.filter, genus=="Galinsoga") )


# Diversity and Species Richness (Observed)
plot_richness(data.filter,x="BeeSpecies" , measures=c("Observed","Shannon"), col="Sample_type")+geom_boxplot() +theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))# to group with metadata
plotdiv2 <- plot_richness(data.filter,x="BeeSpecies", color="Sample_type" , measures=c("Observed","Shannon"))+geom_boxplot() +theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))# to group with metadata

plotdiv2$layers <- plotdiv2$layers[-1]
plotdiv2
barplot(t(as.data.frame(sort(taxa_sums(data.rel.filter), decreasing=T)[1:50])), las=2, horiz = T)

# NMDS ordinations
plot_ordination(data.rel.filter,ordinate(data.rel.filter, method="NMDS", trymax=500, k=2), color="BeeSpecies", shape="Sample_type")+geom_point(size=6)+theme_bw() # , color="Variable1",shape="Variable2"

# Barplots
# needs some dataset preparation for the plot to be ordered by bee species
sample_data(data.rel.filter)$IndName <- interaction(sample_data(data.rel.filter)$BeeSpecies, sample_data(data.rel.filter)$Sample_type,  sample_names(data.rel.filter))
sample_data(data.rel.filter)$IndName <- factor(sample_data(data.rel.filter)$IndName, levels= sort(as.character(sample_data(data.rel.filter)$IndName)))
melted_data <- psmelt(data.rel.filter)
ggplot(melted_data, aes(y=Abundance, x=IndName,fill=family)) + geom_bar( stat= "identity")+theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

### AK: a better visualization is probably to use facets, means and SE?: 

library(reshape2)
library(plyr)

means <- ddply(melted_data, c("BeeSpecies", "family","Sample_type"), summarise,
           mean=mean(Abundance))

plot <- ggplot(data= means, aes(x= family, y=mean, fill= family)) + 
facet_grid(BeeSpecies ~ Sample_type, scales = "free") +
  geom_bar(stat="identity",
           width = 0.8) +                           
  xlab(" ") + ylab("Percentage (%)") + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.y=element_text(size=16, face="bold")) + 
  theme(axis.title.y=element_text(size=16, face="bold", vjust=1)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=11)) +
  theme(legend.position="none")

means.sem <- ddply(melted_data, c("family", "BeeSpecies","Sample_type"), summarise,
                   mean=mean(Abundance), sem=sd(Abundance)/sqrt(length(Abundance)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)

# Add SEM & change appearance of barplot
plotSE <- plot + geom_errorbar(data=means.sem, aes(ymax=upper,  ymin=lower), 
                                width=0.15)
plotSE

# Networks
sample_data(data.rel.filter)$SpcType <- factor(interaction(sample_data(data.rel.filter)$BeeSpecies, sample_data(data.rel.filter)$Sample_type))
merge_spctype <- merge_samples(data.rel.filter, "SpcType")
merge_spctype = transform_sample_counts(merge_spctype, function(x) x/sum(x))

plotweb(t(data.frame(otu_table(merge_spctype))),text.rot=90)

# Subsetting
# if you want to look at individual species of bees, or plants, you can subset your data and then run the plots above again, perhaps on deeper levels for plant species
# e.g. 

Melipona_rufiventris <- subset_samples(data.rel.filter, (BeeSpecies=="Melipona rufiventris" & Sample_type =="pollen" ))
(Melipona_rufiventris = prune_taxa(taxa_sums(Melipona_rufiventris)>0, Melipona_rufiventris))

plot_bar(Melipona_rufiventris, x="IndName", fill="species") +theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25))


Melipona_rufiventris <- subset_samples(data.rel.filter, (BeeSpecies=="Melipona rufiventris" & Sample_type =="honey" ))
(Melipona_rufiventris = prune_taxa(taxa_sums(Melipona_rufiventris)>0, Melipona_rufiventris))

plot_bar(Melipona_rufiventris, x="IndName", fill="species") +theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25))





# carol bee x plant table
data.rel.filter
sample_data(data.rel.filter)$IndName2 <- interaction(sample_data(data.rel.filter)$BeeSpecies, sample_data(data.rel.filter)$Sample_type)
merged_export <- merge_samples(data.rel.filter, "IndName2")
merged_export = transform_sample_counts(merged_export, function(x) x/sum(x))
export <- (cbind(tax_table(merged_export), t(otu_table(merged_export))))
write.table(export,file="export_carol_rm-zero-taxa.txt",sep=";", row.names=F)
 

#---------------------------------
#---------------------------------
#---------------------------------
View(export)


