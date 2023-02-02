
setwd("~/Desktop/metabarcode_cerrado/barcode_cerrado")

table_for_figure <- read.csv("For Fig 6 25 abr 2022.csv")
table_for_figure <- table_for_figure[1:20,]

rownames(table_for_figure) <- table_for_figure$species
barplot(table_for_figure$X, names.arg=table_for_figure$species, las=2)

