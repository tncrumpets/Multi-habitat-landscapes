# This script creates a file entitled "MDT_IntList.csv" which describes
# all interactions happening in each site, and their frequencies.

rm(list=ls(all=TRUE))

PC_PATH <- "/home/alix/Documents/Recherche/LandscapeProject"

DAT_DIR <- paste(PC_PATH, "/MDT/data/30JAN21", sep = "")

setwd(DAT_DIR)
MDT_IntList <- read.csv('MDT_IntList.csv', header = TRUE, stringsAsFactors = FALSE)

MDT_IntListAgg <- subset(MDT_IntList, select = -c(Habitat, Freq))
MDT_IntListAgg <- unique(MDT_IntListAgg)
MDT_IntListAgg$Freq <- NA

# calculate the frequencies of each interaction in each site
for (Int in 1:nrow(MDT_IntListAgg)){
  IndexInt <- which((MDT_IntList$Web == MDT_IntListAgg$Web[Int]) &
                    (MDT_IntList$Site == MDT_IntListAgg$Site[Int]) &
                    (MDT_IntList$Lower_Taxon == MDT_IntListAgg$Lower_Taxon[Int]) &
                    (MDT_IntList$Upper_Taxon == MDT_IntListAgg$Upper_Taxon[Int]))
  MDT_IntListAgg$Freq[Int] <- sum(MDT_IntList$Freq[IndexInt])
  if (length(IndexInt) == 0){
    stop("Program does not find the interaction.")
  }
}

write.table(MDT_IntListAgg, file = "MDT_IntListAgg.csv", sep = ",", row.names = FALSE, quote = FALSE)
