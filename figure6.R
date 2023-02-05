
setwd("~/Desktop/barcode_cerrado")

table_for_figure <- read.csv("new_results_Jun_2021/abundance_table_Melipona rufiventris.csv")
#table_for_figure <- table_for_figure[1:20,]

rownames(table_for_figure) <- table_for_figure$X
barplot(table_for_figure$abundance, names.arg=table_for_figure$X, las=2)

