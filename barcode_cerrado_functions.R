# 
FindMarkerGenBank <- function (taxa, intervals, markers) {
  results = as.data.frame(matrix(nrow=length(taxa), ncol=length(intervals)))
  colnames(results) <- intervals
  
  for(tip_index in 1:length(taxa)) {
    taxon = taxa[tip_index]
    rownames(results)[tip_index] <- taxon
    print(tip_index)
    for(int_index in sequence(ncol(results))) {
      ###################### accessing NCBI ##########################
      id0 <- NA
      while(is.na(id0)[1]) {
        Sys.sleep(1) # workaround to avoid errors when accessing GenBank
        skip_to_next <- FALSE
        tryCatch(id0 <- entrez_search(db="nucleotide", term = paste0(taxon,"[ORGN] AND ",colnames(results)[int_index],"[PDAT]"),retmax = 5000),error = function(e) { skip_to_next <<- TRUE})
        if(is.na(id0)[1]) {
          if(skip_to_next) {
            Sys.sleep(5)
          }
        }
      }  
      if (id0$count == 0) {
        results[tip_index,int_index] <- NA
      } else {
        chuncks <- split(id0$ids, ceiling(seq_along(id0$ids)/300))
        tmp_result <- c()
        for(chunck_index in 1:length(chuncks)) {
          tmp_chuncks <- chuncks[[chunck_index]]
          
          sequences <- NA
          while(is.na(sequences)[1]) {
            Sys.sleep(1)
            skip_to_next <- FALSE
            tryCatch(sequences <- extract_from_esummary(entrez_summary(db="nucleotide", id=tmp_chuncks), "title"),error = function(e) { skip_to_next <<- TRUE})
            if(is.na(sequences)[1]) {
              if(skip_to_next) {
                Sys.sleep(5)
              }
            }
          }
          ##############################################################
          markers = markers # match sequence names with one of these
          matches <- grepl(paste(markers,collapse="|"), sequences, ignore.case=FALSE)
          tmp_result[chunck_index] <- table(matches)["TRUE"]
        }
        results[tip_index,int_index] <- sum(tmp_result[which(!is.na(tmp_result))])
      }  
    }
  }
  return(results)
}


# Function to taxize tip names accoding to NCBI's taxonomy
taxizeNcbi <- function(species, return_species_only = TRUE) {
  gnr_resolve_x <- function(x, return_species_only) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "NCBI"], best_match_only=FALSE)$matched_name)
    ########################################################
    # workaround to find plant names in case of conflict: #
    if(length(tmp.name) > 1) { # If two or more names are matched...
      splitted_name <- strsplit(tmp.name, " ")
      year <- tail(strsplit(tmp.name, " ")[[1]], 1)
      if(!is.na(as.numeric(year))){ # ... and the last character string of the name is a year (NCBI uses author and year of publication to identify conflicting names)
        possible.names <- sapply(1:length(splitted_name), function(x) paste0(splitted_name[[x]][-length(splitted_name[[x]])], collapse = " "))
        new.name.plant <- taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "The International Plant Names Index"], best_match_only=TRUE)$matched_name
        new.name <- tmp.name[which(possible.names %in% new.name.plant)]
      } else { new.name <- tmp.name[1] } # else uses best match from first search
    } else { new.name <- tmp.name[1] } # else uses best match from first search
    ########################################################
    if(return_species_only) {
      if (length(unlist(strsplit(new.name, " "))) == 1) {
        new.name <- paste0("UNMATCHED_", x)
      }
    }
    if(is.null(tmp.name)) {
      new.name <- paste0("UNMATCHED_",x)
    }
    return(new.name)
  }
  new.names <- suppressWarnings(pbapply::pbsapply(species, gnr_resolve_x, return_species_only))
  return(as.character(new.names))
}


setwd("~/Desktop/barcode_cerrado/")
# rm(list=ls())


#-------------------------------
# Checking data availability on GenBank for a list of species
#-------------------------------
# First, let's load a list of species for a quick example
# This is the list of species that occur in the IBGE, an area of cerrado preservation in Brasilia/Brazil.
# We wanted to check which species already had ITS sequences on GenBank and which would have to be prioritized
# for extractions and sequencing.
flora.dir <- paste0(getwd())
species <-  paste(lapply(strsplit(as.character(read.csv(paste0(flora.dir,"/Barcode Lista IBGE Prioridades 1 e 2 revisar.csv"))[,3]), " "), "[[", 1), 
                  lapply(strsplit(as.character(read.csv(paste0(flora.dir,"/Barcode Lista IBGE Prioridades 1 e 2 revisar.csv"))[,3]), " "), "[[", 2), sep=" ")

# We then have to taxize the names according to the NCBI taxonomic backbone
library(taxize)
taxized_species <- taxizeNcbi(species[1:5], return_species_only=TRUE) # test with 5 species
# If return_species_only=TRUE (the default), species that are not found on Genbank are returned marked with "UNMATCHED".
# We can remove those with one line:
taxized_species <- taxized_species[!grepl("UNMATCHED", taxized_species)]

# Next, to check how many entries of a certain marker are on GenBank, we need to specify a temporal interval to search:
intervals = c("1990:2020") # Retrieving entries from 1990 to 2020 (present)

# We also have to specify the names of the markers to search:
markers = c("ITS1", "ITS2", "internal transcribed spacer","internal transcribed spacers", "ITS")

library(rentrez)
# It may take a while to run for many species
results = FindMarkerGenBank(taxa=taxized_species, intervals, markers) 

# For the pollen metabarcode project, we were also interested in classifying priorities for DNA extraction based
# on the number of GenBank entries: 
results$prioridade <- NA
for(i in sequence(nrow(results))) {
  if(is.na(results$`1990:2020`)[i]) {
    results$prioridade[i] <- "prioridade1"
    results$`1990:2020`[i] <- 0
  }
  else if(results$`1990:2020`[i] == 0) {
    results$prioridade[i] <- "prioridade1"
  }  
  else if(results$`1990:2020`[i] == 1) {
    results$prioridade[i] <- "prioridade2"
  }
  else if(results$`1990:2020`[i] == 2) {
    results$prioridade[i] <- "prioridade3"
  }
  else if(results$`1990:2020`[i] > 2) {
    results$prioridade[i] <- "prioridade4"
  }
}

# Done.
# write.csv(results, file="prioridades_revisadas_06.12.2020.csv")
