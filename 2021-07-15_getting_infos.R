
setwd("~/Desktop/metabarcode_cerrado/new_results_Jun_2021")
results <- read.csv("final_otu_vs_bee.csv")
unique_genera <- unique(results$genus)[-grep("spc", unique(results$genus))]
unique_family <- unique(results$family)[-grep("spc", unique(results$family))]
