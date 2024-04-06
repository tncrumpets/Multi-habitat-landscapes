# setting a file for plant species abundances in each habitat, retaining only species involved in interactions
# FU = 0 if NA, FU = 1 is found interacting
# clear workspace
rm(list=ls())

PC_PATH <- "/home/alix/Documents/Recherche/LandscapeProject"

DAT_DIR <- paste(PC_PATH, "/MDT/data/30JAN21", sep = "")

# all observations
setwd(DAT_DIR)
Webs <- read.csv("MDT_IntList.csv", header = TRUE, stringsAsFactors = FALSE)

# plant abundances
PlSurvey <- read.csv("MDT_PlDiv.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = NA)
PlSurvey <- subset(PlSurvey, select = -c(Comments, Name, FlowersperFU_Bud, FlowersperFU_Open, FlowersperFU_Wilt, NA.))
PlSurvey$Veg_Cover <- as.numeric(PlSurvey$Veg_Cover)

Sites <- sort(unique(Webs$Site))
NSites <- length(Sites)

for (s in 1:NSites){
  Web_s <- subset(Webs, Site == Sites[s])
  Hab_s <- sort(unique(Web_s$Habitat))
  Plants_s <- sort(unique(c(Web_s$Lower_Taxon[Web_s$Lower_Guild == "PL"], unique(PlSurvey$Plant_Species[PlSurvey$Site == Sites[s]]))))
  
  CombHabPl <- expand.grid(Hab_s, Plants_s)
  PlData_s <- data.frame(Site = Sites[s], Habitat = CombHabPl$Var1, PlSpecies = CombHabPl$Var2, CP = 0, FU = 0)
  
  for (ph in 1:nrow(PlData_s)){
    p <- PlData_s$PlSpecies[ph]; h <- PlData_s$Habitat[ph]
    Data_sph <- subset(PlSurvey, (Site == Sites[s]) & (Plant_Species == p) & (Habitat == h))
    TestPlSurvey <- (nrow(Data_sph) != 0) # is plant species observed during the plant survey
    IndexPlInt <- which((Web_s$Habitat == h) & (Web_s$Lower_Taxon == p))
    TestPlInt <- (length(IndexPlInt) != 0) # is plant species observed during network sampling
    
    IndexPlData <- which((PlData_s$Habitat == h) & (PlData_s$PlSpecies == p))
    # if species p is observed interacting and during the plant survey
    if ((TestPlInt) & (TestPlSurvey)){
      # cross points
      PlData_s$CP[IndexPlData] <- sum(Data_sph$Cross_Points, na.rm = TRUE)
      if (PlData_s$CP[IndexPlData] != 0){
        PlData_s$CP[IndexPlData] <- PlData_s$CP[IndexPlData]/nrow(Data_sph[!is.na(Data_sph$Cross_Points), ]) # take the average
      }
      else {
        PlData_s$CP[IndexPlData] <- 1
      }
      # floral units
      PlData_s$FU[IndexPlData] <- sum(Data_sph$FU_Open, na.rm = TRUE)
      if (PlData_s$FU[IndexPlData] != 0){
        PlData_s$FU[IndexPlData] <- PlData_s$FU[IndexPlData]/nrow(Data_sph[!is.na(Data_sph$FU_Open), ]) # take the average
      }
    }
    
    # if species p is observed interacting, but not in the plant survey
    if ((TestPlInt) & (!TestPlSurvey)){
      # indicate the minimum number for cross points and floral units
      PlData_s$CP[IndexPlData] <- 1
      PlData_s$FU[IndexPlData] <- 1
    }
    
    # if species p is not observed interacting, but during the plant survey
    if ((!TestPlInt) & (TestPlSurvey)){
      # cross points
      PlData_s$CP[IndexPlData] <- sum(Data_sph$Cross_Points, na.rm = TRUE)
      if (PlData_s$CP[IndexPlData] != 0){
        PlData_s$CP[IndexPlData] <- PlData_s$CP[IndexPlData]/nrow(Data_sph[!is.na(Data_sph$Cross_Points), ]) # take the average
      }
      else {
        if (!is.na(mean(Data_sph$Veg_Cover, na.rm = TRUE))){
          PlData_s$CP[IndexPlData] <- 1
        }
      }
      # floral units
      PlData_s$FU[IndexPlData] <- sum(Data_sph$FU_Open, na.rm = TRUE)
      if (PlData_s$FU[IndexPlData] != 0){
        PlData_s$FU[IndexPlData] <- PlData_s$FU[IndexPlData]/nrow(Data_sph[!is.na(Data_sph$FU_Open), ]) # take the average
      }
    }
  }
  write.table(PlData_s, file = paste0('plant_abundances', Sites[s], "_PlAbundInt.csv"), sep = ",", quote = FALSE, row.names = FALSE)
}
