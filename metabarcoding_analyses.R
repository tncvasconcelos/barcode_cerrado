library(phyloseq)
library("data.table")
library(bipartite)
library(ggplot2)
library(vegan)
library(reshape2)
library(plyr)
library(lemon)
library('ggthemes')
library(circlize)
library(viridis)


# load functions available at https://github.com/chiras/metabarcoding_pipeline
source('./metabarcoding_tools_0-1a.R')

# load in data
setwd("/Users/ra39huv/TMP/Basespace/ITS2/2023/ITS2_Aline_Martins/old_binf")
tax <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T,sep=",")))
otu <- otu_table(read.table("asv_table.merge.txt"), taxa_are_rows=T)
map <- sample_data(read.table("sample_list.csv", sep=";", header=T,row.names=1))

# check that all samples have metadata
sum(!(sample_names(map) %in% sample_names(otu)))
sum(!(sample_names(otu) %in% sample_names(map)))

# combine datasets
(data <- merge_phyloseq(otu,tax,map))
data.comp <- data

# check incomplete taxonomic assignments
data.comp <- propagate_incomplete_taxonomy(data.comp)
tail(tax_table(data.comp))

#Remove taxa that received no taxnomic classification and Marchantiales
dataset.comp.filter = data.comp
(dataset.comp.filter = subset_taxa(dataset.comp.filter, phylum=="p:Streptophyta"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, family!="p:Streptophyta_spc_spc"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, species!="k:Viridiplantae_spc_spc_spc_spc_spc"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, species!="p:Streptophyta_spc_spc_spc_spc"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, order!="o:Marchantiales"))
(dataset.comp.filter = subset_taxa(dataset.comp.filter, kingdom!=""))

#Make taxa labels nice for plots
tail(tax_table(dataset.comp.filter))
dataset.comp.filter <- replace_tax_prefixes(dataset.comp.filter)

# collating ASVs to species
(data.species <- tax_glom(dataset.comp.filter,taxrank="species"))
taxa_names(data.species) <- tax_table(data.species)[,"species"]

tail(tax_table(data.species))

# transformation to relative read abundance data
data.species.rel = transform_sample_counts(data.species, function(x) x/sum(x))

# low abundance filtering
otu_table(data.species.rel)[otu_table(data.species.rel)<0.01 ]<-0
otu_table(data.species)[otu_table(data.species.rel)<0.01 ]<-0
data.filter		= prune_taxa(taxa_sums(data.species)>0, data.species)
(data.rel.filter = prune_taxa(taxa_sums(data.species.rel)>0, data.species.rel))

# reduce dataset to Meliponini only and split to datasets 
data.rel.filter <- subset_samples(data.rel.filter, include=="yes") 
data.filter <- subset_samples(data.filter, include=="yes") 

(data.filter = prune_taxa(taxa_sums(data.rel.filter)>0, data.filter))
(data.rel.filter = prune_taxa(taxa_sums(data.rel.filter)>0, data.rel.filter))

pollen <- subset_samples(data.rel.filter, Sample_type=="pollen")
honey <- subset_samples(data.rel.filter, Sample_type=="honey")
dir.create("plots_meliponini")

# Diversity and Species Richness (Observed)
plotdiv2 <- plot_richness(data.filter,x="BeeSpecies", color="Sample_type" , measures=c("Shannon"))+
  geom_boxplot() + 
  theme_bw()+
  scale_colour_grey(start = 0.4,end = 0.6)

plotdiv2$layers <- plotdiv2$layers[-1]

plotdiv2 = plotdiv2 + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position=c(0.875,0.1), 
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank()	)+
  theme(axis.text.x = element_text( family = "sans", face = "italic"))+
  ylab("Shannon Diversity")+
  xlab("")

plotdiv2

pdf("plots_meliponini/diversity.pdf", width=4.5, height=5)
plotdiv2
dev.off()

# Shannon and Spec Richness Kruskal-Wallis tests
pollen_test <- cbind(sample_data(pollen),diversity(t(otu_table(pollen))),specnumber(t(otu_table(pollen))))

kruskal.test(pollen_test[,9],g=pollen_test[,2])
kruskal.test(pollen_test[,8],g=pollen_test[,2])

honey_test <- cbind(sample_data(honey),diversity(t(otu_table(honey))),specnumber(t(otu_table(honey))))
kruskal.test(honey_test[,9],g=honey_test[,2])
kruskal.test(honey_test[,8],g=honey_test[,2])

### PERMANOVAs

# overall
meta <- as.data.frame(cbind(data.frame(sample_data(data.rel.filter)[,c("BeeSpecies")])[,1],data.frame(sample_data(data.rel.filter)[,"Sample_type"])[,1]))
colnames(meta) <- c("BeeSpecies","Sample_type")
(adonis2( vegdist(t(otu_table(data.rel.filter))) ~ meta$BeeSpecies * meta$Sample_type))

# pollen
meta <- as.data.frame(cbind(data.frame(sample_data(pollen)[,c("BeeSpecies")])[,1],data.frame(sample_data(pollen)[,"Sample_type"])[,1]))
colnames(meta) <- c("BeeSpecies","Sample_type")
(adonis2( vegdist(t(otu_table(pollen))) ~ meta$BeeSpecies))

# honey
meta <- as.data.frame(cbind(data.frame(sample_data(honey)[,c("BeeSpecies")])[,1],data.frame(sample_data(honey)[,"Sample_type"])[,1]))
colnames(meta) <- c("BeeSpecies","Sample_type")
(adonis2( vegdist(t(otu_table(honey))) ~ meta$BeeSpecies))

###  Ordination plots

pdf("plots_meliponini/ordination.pdf", width=7, height=7)

# pollen
plot_ordination(data.rel.filter,ordinate(pollen, method="NMDS", trymax=500, k=2), color="BeeSpecies", shape="BeeSpecies")+
  geom_point(size=6)+
  theme_bw()+
  ggtitle("Pollen") +
  scale_colour_grey(start = 0.7,end = 0)+
  theme(legend.position=c(0.82,0.9), 
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        plot.title=element_text(face="bold", size=16, hjust = 0.02, margin=margin(t=20,b=-25)))+
  theme(legend.text = element_text( family = "sans", face = "italic"))

# honey
plot_ordination(data.rel.filter,ordinate(honey, method="NMDS", trymax=500, k=2), color="BeeSpecies", shape="BeeSpecies")+
  geom_point(size=6)+
  theme_bw()+
  ggtitle("Honey") +
  scale_colour_grey(start = 0.7,end = 0)+
  theme(legend.position=c(0.82,0.9), 
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        plot.title=element_text(face="bold", size=16, hjust = 0.02, margin=margin(t=20,b=-25)))+
  theme(legend.text = element_text( family = "sans", face = "italic"))

dev.off()

### Creating the habitat barplots
toptaxa = names(tail(sort(rowSums(otu_table(data.rel.filter))),n=30))
data.rel.filter.top  <- subset_taxa(data.rel.filter, species %in% toptaxa)
melted_data <- psmelt(data.rel.filter.top)


plant_meta <- read.csv("plant_traits.csv", header=T, sep=";")

# not elegant, but working mapping. takes a bit
melted_data$habitat = c()
for (i in 1:length(melted_data$species)){
  
  print(paste(i, " ",melted_data$species[i], "-- >" ,plant_meta[plant_meta$species == melted_data$species[i],2]))
  melted_data$habitat[i] <- plant_meta[plant_meta$species == melted_data$species[i],2][1]
}

means <- ddply(melted_data, c("BeeSpecies", "species","Sample_type","habitat"), summarise,
               mean=mean(Abundance))

plot <- ggplot(data= means, aes(x= species, y=mean, fill= habitat)) + 
  facet_grid(BeeSpecies ~ Sample_type, scales='fixed') + 
  geom_bar(stat="identity",
           width = 0.8) + 
  xlab(" ") + 
  ylab("Relative Read Abundance (%)") + 
  theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.y=element_text(size=10)) + 
  theme(axis.title.y=element_text(size=13, vjust=1)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=11, face="italic")) +
  theme(legend.position=c(0.5,0.93), legend.direction = "horizontal", legend.title= element_blank(),legend.box.background = element_rect(colour = "black")) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,2), "cm"))+
  theme(panel.spacing.x = unit(15, "pt"))+
  theme(panel.border=element_blank(), axis.line=element_line(), strip.text.y=element_text(face="italic"))

# Add SEM & change appearance of barplot
means.sem <- ddply(melted_data, c("species", "BeeSpecies","Sample_type","habitat"), summarise,
                   mean=mean(Abundance), sem=sd(Abundance)/sqrt(length(Abundance)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)

plotSE <- plot + geom_errorbar(data=means.sem, aes(ymax=upper,  ymin=lower), 
                               width=0.15)
pdf("plots_meliponini/byhabitat_top.pdf", width=13, height=8)

plotSE
dev.off()


### Network stats
sample_data(pollen)$SpcType <- factor(interaction(sample_data(pollen)$BeeSpecies))
merge_spctype <- merge_samples(pollen, "BeeSpecies")
d_circ 		= prune_taxa(taxa_sums(merge_spctype)>0.05, merge_spctype)
mat=t(otu_table(d_circ))

dfun(t(mat))
H2fun(t(mat),H2_integer=F)

### Network plotting

bees=c("Melipona rufiventris","Scaptotrigona sp","Tetragonisca angustula")

{
  families_sorted <- names(sort(table(tax_table(d_circ)[,c("family")]),decreasing=T))
  families_colors <- data.frame(cbind(families_sorted,viridis(length(families_sorted))))
  #families_colors <- data.frame(cbind(families_sorted,plasma(length(families_sorted))))
  #families_colors <- data.frame(cbind(families_sorted,cividis(length(families_sorted))))
  families_colors <- data.frame(cbind(families_sorted,turbo(length(families_sorted))))
  plant_spc = (tax_table(d_circ)[,c("family","species")])
  grid.col.pollen = c(rand_color(length(plant_spc)))
  
  rownames(families_colors) <- families_colors[,1]
  group = structure(c(t(plant_spc[,"family"])), names =c(t(plant_spc[,"species"])))
  
  for (i in 1:length(plant_spc[,1])){
    plant_spc[i,1] <- families_colors[c(plant_spc[i,1]),2]
  } # plant species to family assignments
  
  grid.col.bees = c(gray.colors(length(bees), start = 0.7, end = 0))
  grid.col = c(plant_spc[,1], grid.col.bees)
  group = c(group,structure(c("Bees","Bees","Bees"), names= bees))
} # color and group assignments

mat = t(otu_table(d_circ))

# top taxa naming
{
  # only plot 20 most abundant, rest in numbers
  toptaxa <- data.frame(rra = sort(rowSums(mat), decreasing=T), display = c(rep(0,20), seq(1,length(rowSums(mat))-20)), track=NA, family=NA)
  
  for (i in 1:20){
    toptaxa[i,2]<- rownames(toptaxa)[i]
  }
  
  for (i in 1:length(toptaxa[,1])){
    toptaxa[i,4]<- c(group[rownames(toptaxa)[i]])
  }
} 

# track assignments for species
{
  head(toptaxa, n=25)
  toptaxa$track = rep(9,length(toptaxa[,1]))
  toptaxa[toptaxa[,2]=="Hedyosmum brasiliense",3] <- 8
  toptaxa[toptaxa[,2]=="Maprounea guianensis",3] <- 7
  toptaxa[toptaxa[,2]=="Stryphnodendron spc",3] <- 8
  toptaxa[toptaxa[,2]=="Struthanthus/Psittacanthus spc",3] <- 6
  toptaxa[toptaxa[,2]=="15",3] <- 7
  toptaxa[toptaxa[,2]=="1",3] <- 6
  toptaxa[toptaxa[,2]=="14",3] <- 6
  toptaxa[toptaxa[,2]=="33",3] <- 5
  toptaxa[toptaxa[,2]=="Miconia stenostachya",3] <- 8
  toptaxa[toptaxa[,2]=="Miconia leucocarpa",3] <- 7
  toptaxa[toptaxa[,2]=="21",3] <- 8
  toptaxa[toptaxa[,2]=="8",3] <- 4
  toptaxa[toptaxa[,2]=="9",3] <- 6
  toptaxa[toptaxa[,2]=="29",3] <- 5
  toptaxa[toptaxa[,2]=="31",3] <- 4
  toptaxa[toptaxa[,2]=="10",3] <- 6
  toptaxa[toptaxa[,2]=="6",3] <- 5
  toptaxa[toptaxa[,2]=="32",3] <- 5
  
  toptaxa[toptaxa[,2]=="23",3] <- 7
  toptaxa[toptaxa[,2]=="19",3] <- 6
  
  toptaxa[toptaxa[,2]=="Myrcia linearifolia",3] <- 8
  toptaxa[toptaxa[,2]=="Syzygium cumini",3] <- 7
  toptaxa[toptaxa[,2]=="18",3] <- 8
  toptaxa[toptaxa[,2]=="Myrcia pinifolia",3] <- 6
  toptaxa[toptaxa[,2]=="5",3] <- 8
  toptaxa[toptaxa[,2]=="27",3] <- 8
  toptaxa[toptaxa[,2]=="16",3] <- 7
  toptaxa[toptaxa[,2]=="25",3] <- 7
  
  toptaxa[toptaxa[,2]=="Guapira graciliflora",3] <- 6
  
  toptaxa[toptaxa[,2]=="28",3] <- 7
  
  toptaxa[toptaxa[,2]=="Piper aduncum",3] <- 8
  
  toptaxa[toptaxa[,2]=="11",3] <- 7
  toptaxa[toptaxa[,2]=="Rubus urticifolius",3] <- 6
  toptaxa[toptaxa[,2]=="35",3] <- 7
  toptaxa[toptaxa[,2]=="Matayba guianensis",3] <- 8
  toptaxa[toptaxa[,2]=="20",3] <- 6
  
  toptaxa[toptaxa[,2]=="Tapirira guianensis",3] <- 7
  
  toptaxa[toptaxa[,2]=="17",3] <- 6
  toptaxa[toptaxa[,2]=="34",3] <- 7
  toptaxa[toptaxa[,2]=="Baccharis dracunculifolia",3] <- 8
  toptaxa[toptaxa[,2]=="22",3] <- 6
} 

pdf("plots_meliponini/network_numbered.pdf", width=13, height=8)

# plotting basic network plot
{
  circos.clear()
  chordDiagram(mat, annotationTrack = "grid",
               grid.col = grid.col,
               transparency = 0.4, 
               big.gap = 5,
               group=group,
               preAllocateTracks = list(list(track.height = 0.02),list(track.height = 0.02),list(track.height = 0.02),list(track.height = 0.02),
                                        list(track.height = 0.02),list(track.height = 0.02),list(track.height = 0.02),list(track.height = 0.02),list(track.height = 0.02)))
} 

# plotting bee names
trackid=c(9,8,9)
for (i in seq_len(ncol(mat))) { # use for loop to label each sector
  myFactor <- colnames(mat)[i] # assuming this defines the sectors
  myCol <- grid.col.bees[i] # defined in the question
  circos.trackPlotRegion(track.index = trackid[i], factor = myFactor,
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           sector.name = get.cell.meta.data("sector.index")
                           niceFacing = T
                           circos.text(mean(xlim), 
                                       ylim[1] + .1, 
                                       cex=1,
                                       sector.name,
                                       font=4,
                                       col = myCol,
                                       facing = "bending",
                                       niceFacing = T)
                         },
                         bg.border = NA)
} 

# plotting plant species names/numbers
for (i in seq_len(nrow(mat))) { # use for loop to label each sector
  myFactor <- rownames(mat)[i] # assuming this defines the sectors
  myCol <- c(plant_spc[myFactor,"family"]) # defined in the question
  #trackid=trackid-1
  #if (trackid < 1) {trackid=9}
  circos.trackPlotRegion(track.index = toptaxa[myFactor,3], factor = myFactor,
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           sector.name =  toptaxa[get.cell.meta.data("sector.index"),2] 
                           niceFacing = T
                           circos.text(mean(xlim), 
                                       ylim[1] + .1, 
                                       cex=0.7,
                                       sector.name,
                                       col = myCol,
                                       font = 3,
                                       facing = "bending",
                                       niceFacing = F)
                           #   circos.axis(h = "top", 
                           #     labels.cex = 0.5,
                           #     major.tick.percentage = 0.2, 
                           #    sector.index = sector.name, 
                           #    track.index = 1)
                         },
                         bg.border = NA)
} 

# plotting families
{families_colors$track <- rep(1,length(families_colors[,1]))
  
  families_colors["Myrtaceae",]$track <- 3
  families_colors["Nyctaginaceae",]$track <- 3
  families_colors["Moraceae",]$track <- 3
  families_colors["Boraginaceae",]$track <- 3
  families_colors["Piperaceae",]$track <- 3
  families_colors["Malpighiaceae",]$track <- 3
  families_colors["Fabaceae",]$track <- 3
  families_colors["Cannabaceae",]$track <- 4
  families_colors["Pinaceae",]$track <- 4
  families_colors["Rosaceae",]$track <- 4
  families_colors["Aquifoliaceae",]$track <- 3
  families_colors["Sapindaceae",]$track <- 3
  families_colors["Anacardiaceae",]$track <- 4
  families_colors["Asteraceae",]$track <- 4
  families_colors["Chrysobalanaceae",]$track <- 5
  families_colors["Clusiaceae",]$track <- 3
  families_colors["Combretaceae",]$track <- 4} # track assignments for families
for (i in 1:length(families_colors[,1])){
  
  highlight.sector(row.names(toptaxa)[toptaxa[,4]==families_colors[i,1]], track.index = families_colors$track[i], col="white",
                   text = families_colors[i,1], cex = 1.1, font = 2, text.col = families_colors[i,2], niceFacing = F, facing = "bending")
  
  highlight.sector(row.names(toptaxa)[toptaxa[,4]==families_colors[i,1]], track.index = 2,  col = families_colors[i,2],lwd = 4)
  
} 
dev.off()


#### Export data
write.table(file="export.pollen.csv",cbind(tax_table(pollen),otu_table(pollen)), sep=";")
write.table(file="export.honey.csv",cbind(tax_table(honey),otu_table(honey)), sep=";")
write.table(file="export.diversity-pollen.csv",cbind(specnumber(t(otu_table(pollen))),round(diversity(t(otu_table(pollen))), digits=3),sample_data(pollen)), sep=";")
write.table(file="export.diversity-honey.csv",cbind(specnumber(t(otu_table(honey))),round(diversity(t(otu_table(honey))), digits=3),sample_data(honey)), sep=";")
