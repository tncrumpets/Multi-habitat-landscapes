
setwd("C:/Users/user/Dropbox/Kate/Talya_FD")

# fill in your gbif.org credentials 
user <- "katemaia" # your gbif.org username 
pwd <- "kate0021" # your gbif.org password
email <- "kate.spirogyra@gmail.com" # your email

library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_

# match the names 
gbif_taxon_keys <- 
  readr::read_csv("likely_plant_species_search.csv") %>%
  pull("Taxon name") %>% # use fewer names if you want to just test 
  taxize::get_gbifid_(method = "backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT") %>% # get only accepted (REMOVED: & status == "ACCEPTED") and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys
# gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!
occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred("country", "GB"),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV",
  user = user,pwd = pwd,email = email
)
