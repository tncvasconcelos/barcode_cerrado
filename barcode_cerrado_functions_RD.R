# 
taxa=taxized_species
FindMarkerGenBank <- function (taxa, intervals, markers) {
  results = as.data.frame(matrix(nrow=length(taxa), ncol=length(intervals)))
  colnames(results) <- intervals

  #tip_index = 3
  for (tip_index in 1:length(taxa)) {
    taxon = taxa[tip_index]
    rownames(results)[tip_index] <- taxon
    print(paste(tip_index,length(taxa),sep='/'))
    flush.console()
    
    #int_index = 1
    for (int_index in 1:ncol(results)) {
      ###################### accessing NCBI ##########################
      id0 <- NA
      while(is.na(id0)[1]) {
        Sys.sleep(1) # workaround to avoid errors when accessing GenBank
        skip_to_next <- FALSE
        tryCatch(
          id0 <- entrez_search(db="nucleotide", term = paste0(taxon,"[ORGN] AND ",intervals[int_index],"[PDAT]"),retmax = 5000),
          error = function(e) { skip_to_next <- TRUE}
        )
        if(is.na(id0)[1]) {
          if(skip_to_next) {
            Sys.sleep(5)
          }
        }
      }  
      if (id0$count == 0) {
        results[tip_index,int_index] <- NA
      } else {
        print(id0$ids)
        flush.console()
        chuncks <- split(id0$ids, ceiling(seq_along(id0$ids)/300))
        tmp_result <- c()
        #chunck_index = 1
        for (chunck_index in 1:length(chuncks)) {
          tmp_chuncks <- chuncks[[chunck_index]]
          sequences <- NA
          while(is.na(sequences)[1]) {
            print(paste0('tentou ',tip_index))
            flush.console()
            Sys.sleep(1)
            skip_to_next <- FALSE
            tryCatch(
              sequences = extract_from_esummary(entrez_summary(db='nucleotide',id=tmp_chuncks),'title'),
              error = function(e) { skip_to_next <<- TRUE}
            )
            if(is.na(sequences)[1]) {
              if(skip_to_next) {
                Sys.sleep(5)
              }
            }
          }
          #markers = markers # match sequence names with one of these
          matches <- grepl(paste(markers,collapse="|"), sequences, ignore.case=FALSE)
          tmp_result[chunck_index] <- paste0(sequences[matches],collapse='道')
          print(tmp_result[chunck_index])
          flush.console()
        }
        results[tip_index,int_index] <- paste0(tmp_result,collapse='道')
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


#setwd("~/Desktop/barcode_cerrado/")
setwd("~/Documents/LabFer/Serrapilheira/2021-04-30")
# rm(list=ls())


#-------------------------------
# Checking data availability on GenBank for a list of species
#-------------------------------
# First, let's load a list of species for a quick example
# This is the list of species that occur in the IBGE, an area of cerrado preservation in Brasilia/Brazil.
# We wanted to check which species already had ITS sequences on GenBank and which would have to be prioritized
# for extractions and sequencing.
d = read.csv('IBGE_specieslink_RD.csv',sep='\t')
#d = read.csv('Barcode Lista IBGE Prioridades 1 e 2 revisar.csv')
names(d)[3] <- 'full.name'
names(d)[7] <- 'scientific.name'
head(d)

# método 1, com expressão regular (regex)
re = regexpr('\\S+\\s+\\S+\\K ',d$scientific.name,perl=T) # posição do primeiro espaço após o segundo nome
re[re < 0] = nchar(d$scientific.name)[re < 0] # onde só tem dois nomes, pega o comprimento do texto
sp1 = trimws(substr(d$scientific.name,1,re)) # retorna só os dois primeiros nomes

# método 2, com a biblioteca stringr
library(stringr)
sp2 = word(d$scientific.name,1,2)

# método 3, com lapply
sp3 = paste(lapply(strsplit(d$scientific.name, ' '), '[[', 1),
            lapply(strsplit(d$scientific.name, ' '), '[[', 2))

# provando que os 3 métodos são equivalentes
any(sp1 != sp2) # false
any(sp1 != sp3) # false
any(sp2 != sp3) # false
d[which(d$scientific.name != sp1),]

species = sp1

# We then have to taxize the names according to the NCBI taxonomic backbone
library(taxize)
taxized_species = species
#taxized_species <- taxizeNcbi(species[1:5], return_species_only=TRUE) # test with 5 species
taxized_species <- taxizeNcbi(species, return_species_only=TRUE)
# If return_species_only=TRUE (the default), species that are not found on Genbank are returned marked with "UNMATCHED".
species[1:5]
taxized_species
# We can remove those with one line:
taxized_species1 <- taxized_species[!grepl("UNMATCHED", taxized_species)]
table(taxized_species != species)
species[taxized_species != species]
taxized_species[taxized_species != species]

# Next, to check how many entries of a certain marker are on GenBank, we need to specify a temporal interval to search:
intervals = c("1990:2021") # Retrieving entries from 1990 to 2020 (present)

# We also have to specify the names of the markers to search:
#markers = c("ITS1", "ITS2", "internal transcribed spacer","internal transcribed spacers", "ITS")
markers = c('internal transcribed spacer','internal transcribed spacers', 'ITS')

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
