# This script creates a file entitled "MDT_IntList.csv" which describes all interactions happening in each
# habitat of each site, and their frequencies.

rm(list=ls(all=TRUE))

# PC_PATH <- "/media/asauve/DATAPART1/landscape_project"
PC_PATH <- "/home/alix/Documents/Recherche/LandscapeProject"

DAT_DIR <- paste(PC_PATH, "/MDT/data/30JAN21", sep = "")

setwd(DAT_DIR)
MDT_AllObs <- read.csv('all_web_interactions_20210130.csv', header = TRUE, stringsAsFactors = FALSE)

# remove lines with a species to remove (marked as "Unknown_spREMOVEFROMWEB")
RmLines <- which((MDT_AllObs$Lower_Taxon == 'Unknown_sp_REMOVEFROMWEB') | (MDT_AllObs$Upper_Taxon == 'Unknown_sp_REMOVEFROMWEB'))
if (length(RmLines) > 0){
  MDT_AllObs <- MDT_AllObs[-RmLines, ]
}
# remove lines with lower/upper taxa as "NA"
RmLines <- which((is.na(MDT_AllObs$Lower_Taxon)) | (is.na(MDT_AllObs$Upper_Taxon)))
if (length(RmLines) > 0){
  MDT_AllObs <- MDT_AllObs[-RmLines, ]
}

MDT_IntList <- subset(MDT_AllObs, select = c(Web, Site, M_D_T, Habitat, Lower_Taxon, Upper_Taxon)) # removing some columns
MDT_IntList <- unique(MDT_IntList) # remove duplicates of interactions
MDT_IntList$Freq <- NA # create a column for interactions frequencies
# calculate the interaction frequencies
ListIndexInt <- vector("list", nrow(MDT_IntList))
AggIndexInt <- c()
for (int in 1:nrow(MDT_IntList)){
  index_int <- which((MDT_AllObs$Web == MDT_IntList$Web[int]) &
                     (MDT_AllObs$Site == MDT_IntList$Site[int]) &
                     (MDT_AllObs$Habitat == MDT_IntList$Habitat[int]) &
                     (MDT_AllObs$Lower_Taxon == MDT_IntList$Lower_Taxon[int]) &
                     (MDT_AllObs$Upper_Taxon == MDT_IntList$Upper_Taxon[int]))

  MDT_IntList$Freq[int] <- length(index_int)
  if (length(index_int) == 0){
    stop("Program does not find the interaction.")
  }
}

# check that for each site, all habitats are in the data
Sites <- sort(unique(MDT_IntList$Site))
NSites <- length(Sites)
for (s in 1:NSites){
  Hab_s <- unique(MDT_IntList$Habitat[MDT_IntList$Site == Sites[s]])
  MDT <- unique(MDT_IntList$M_D_T[MDT_IntList$Site == Sites[s]])
  MDT <- ifelse((MDT == "monad"), 1, ifelse((MDT == "dyad"), 2, 3))
  if (length(Hab_s) != MDT){
    stop("Problem with the habitats of sites s: Not all sites are present in the dataset.")
  }
}

# add two columns to describe the guild to which each taxon belongs
MDT_IntList$Upper_Taxon <- as.character(MDT_IntList$Upper_Taxon); MDT_IntList$Lower_Taxon <- as.character(MDT_IntList$Lower_Taxon)
MDT_IntList$Lower_Guild <- sapply(strsplit(MDT_IntList$Web, "_"), "[[", 1)
MDT_IntList$Upper_Guild <- sapply(strsplit(MDT_IntList$Web, "_"), "[[", 2)

# adding two columns (one for the upper guild and one for the lower guild) for taxon's trophic level (hereafter TL)
# if PL, then TL = 1; if LM or FV or SF or CP, then TL = 2; if P, then TL = 3
MDT_IntList$TL_lower_taxon <- NA; MDT_IntList$TL_upper_taxon <- NA
MDT_IntList$TL_lower_taxon[MDT_IntList$Lower_Guild == "PL"] <- 1 # plants (PL)
MDT_IntList$TL_upper_taxon[(MDT_IntList$Upper_Guild == "FV") | (MDT_IntList$Upper_Guild == "LM") | (MDT_IntList$Upper_Guild == "SF") | (MDT_IntList$Upper_Guild == "CP")] <- 2 # flower visitors (FV), leaf miners (LM), seed-feeding insect (SF) and caterpillars (CP) as upper guild
MDT_IntList$TL_lower_taxon[(MDT_IntList$Lower_Guild == "LM") | (MDT_IntList$Lower_Guild == "CP") | (MDT_IntList$Lower_Guild == "SF")] <- 2 # leaf miners as lower guild (LM)
MDT_IntList$TL_upper_taxon[(MDT_IntList$Upper_Guild == "SFP") | (MDT_IntList$Upper_Guild == "P")] <- 3 # parasitoids (SFP and P)

# adding two columns (one for the upper guild, and one for the lower guild) for taxon's life stage (hereafter LS)
# if {CP | LM | SF | P | SFP}, then LS = 'larva'; if {FV}, then LS = 'adult', NA for plant species
MDT_IntList$LS_lower_taxon <- NA; MDT_IntList$LS_upper_taxon <- NA
MDT_IntList$LS_upper_taxon[MDT_IntList$Upper_Guild == "FV"] <- 'adult'
MDT_IntList$LS_upper_taxon[(MDT_IntList$Upper_Guild == "P") | (MDT_IntList$Upper_Guild == "SFP") | (MDT_IntList$Upper_Guild == "LM") | (MDT_IntList$Upper_Guild == "SF") | (MDT_IntList$Upper_Guild == "CP")] <- 'larva'
MDT_IntList$LS_lower_taxon[(MDT_IntList$Lower_Guild == "LM") | (MDT_IntList$Lower_Guild == "CP")] <- 'larva'

write.table(MDT_IntList, file = "MDT_IntList.csv", sep = ",", row.names = FALSE, quote = FALSE)
