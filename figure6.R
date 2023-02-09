rm(list=ls())


setwd("~/Desktop/barcode_cerrado")
library(ggplot2)
library(gridExtra)

#table_for_figure <- read.csv("new_results_Jun_2021/abundance_table_Melipona rufiventris.csv")
#table_for_figure <- table_for_figure[1:20,]

tables_for_figures <- list.files("datasets", full.names = T)[grep("abundance_table", list.files("datasets"))]
bee_species <- gsub(paste0(c("datasets/abundance_table_",".csv"),collapse="|"),"",tables_for_figures)
tables_for_figures <- lapply(tables_for_figures, read.csv)
names(tables_for_figures) <- bee_species
metadata <- read.csv("datasets/For Fig 6 8 Feb 2023.csv")
metadata$habitat[which(is.na(metadata$habitat))] <- "indet species"

###################
# Bee species 1:
one_table <- tables_for_figures[[1]]
rownames(one_table) <- one_table$X
one_table <- merge(one_table, metadata, by.x="X", by.y="species")

pal <- c("#DFC27D","#B8E186","#80CDC1","#BABABA","#B2ABD2")
bee_species[1]
bee1 <- ggplot(data=one_table, aes(x=reorder(X, -abundance), y=abundance, fill=habitat)) +
  geom_bar(stat="identity")+theme_minimal() +
  ggtitle(paste0("a) ", bee_species[1])) +
  theme(plot.title = element_text(size = 0.2)) +
  #theme_bw() +
  #theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10),
  #      axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(text = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title=element_text(size=15)) +
  theme(plot.subtitle=element_text(size=15)) +
  theme(axis.text=element_text(size=8)) +
  scale_fill_manual(values=pal)


#pal <- hcl.colors(8, palette = "Viridis", alpha = 1)

#################
# Bee species 2:
one_table <- tables_for_figures[[2]]
rownames(one_table) <- one_table$X
one_table <- merge(one_table, metadata, by.x="X", by.y="species")

pal <- c("#DFC27D","#B8E186","#80CDC1","#BABABA","#B2ABD2")

bee_species[2]
bee2 <- ggplot(data=one_table, aes(x=reorder(X, -abundance), y=abundance, fill=habitat)) +
  geom_bar(stat="identity")+theme_minimal() +
  ggtitle(paste0("a) ", bee_species[2])) +
  theme(plot.title = element_text(size = 0.2)) +
  #theme_bw() +
  #theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10),
  #      axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(text = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title=element_text(size=15)) +
  theme(plot.subtitle=element_text(size=15)) +
  theme(axis.text=element_text(size=8)) +
  scale_fill_manual(values=pal)


#pal <- hcl.colors(8, palette = "Viridis", alpha = 1)

################
# Bee species 3:
one_table <- tables_for_figures[[3]]
rownames(one_table) <- one_table$X
one_table <- merge(one_table, metadata, by.x="X", by.y="species")

pal <- c("#DFC27D","#B8E186","#80CDC1","#BABABA","#B2ABD2")

bee_species[3]
bee3 <- ggplot(data=one_table, aes(x=reorder(X, -abundance), y=abundance, fill=habitat)) +
  geom_bar(stat="identity")+theme_minimal() +
  ggtitle(paste0("a) ", bee_species[3])) +
  theme(plot.title = element_text(size = 0.2)) +
  #theme_bw() +
  #theme(legend.position="none", text = element_text(size = 15), axis.text.y = element_text(size = 10),
  #      axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(text = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title=element_text(size=15)) +
  theme(plot.subtitle=element_text(size=15)) +
  theme(axis.text=element_text(size=8)) +
  scale_fill_manual(values=pal)


#pal <- hcl.colors(8, palette = "Viridis", alpha = 1)
#pal <- c("#DFC27D","#B8E186","#80CDC1","#B2ABD2")


pdf("test_figure6.pdf",height=12, width=8)
grid.arrange(bee1,bee2,bee3,
             nrow=3, ncol=1)
dev.off()


