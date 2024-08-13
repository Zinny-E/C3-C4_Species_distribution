##########################################################
# This script contain nutnet cover 2007-2023
# removes labels that are not needed for the analysis and 
#fixes some of trait  and category issues.
#####cleaning and setting dataset for anaylsis of percentage cover
##########################################################
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)

setwd("/Users/zinny/git/TXArXSN_nutnet/scripts")
# Load full_cover data
nutnet.cover <- read.csv("../data/full-cover_2023-11-07.csv")
site_cordinates <- read.csv("../data/sites_2023-11-07.csv")
#load TERN. dataset
TERN_ps.pathway<- 
  read_xlsx("../Photosynthetic_pathways_of_plant_species_in_TERN_plots.xlsx")

#convert to csv
write.csv(TERN_ps.pathway, "../tern.pspathway.csv", row.names = FALSE)
tern.pspathway<- read.csv("../tern.pspathway.csv")

#load China trait database
china.pspathway <- read.csv("../china_plant_trait_database/Photo Pathway.csv")
china.speciesname <- read.csv("../china_plant_trait_database/Species translations.csv")
china.taxonomicname <- read.csv("../china_plant_trait_database/Taxonomic standardisation.csv")



####converting taxon in nutnet cover to lowercase and underscore
# Convert only the first letter of each entry in the species_name column to 
#uppercase
nutnet.cover$Taxon <- sapply(nutnet.cover$Taxon, function(x) {
  x <- tolower(x) # Ensure the string is in lowercase first
  # Convert first letter to uppercase
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))) 
})

#convert the uppercase rows in the Family
# convert a string to sentence case
toSentenceCase <- function(s) {
  # Convert the entire string to lowercase first
  s <- tolower(s)
  # Then capitalize the first letter of each sentence
  s <- sub("^(.)", "\\U\\1", s, perl = TRUE)
  return(s)
}
nutnet.cover$Family <- sapply(nutnet.cover$Family, toSentenceCase)



##removes "_" for species name
tern.pspathway$Species.name <- gsub("_", " ", tern.pspathway$Species.name)

##Finding common species
#common_species <- intersect(nutnet.cover$Taxon, tern.pspathway$Species.name)

##removing leading and trailing whitespaces in species names. 
nutnet.cover$Taxon<- trimws(tolower(nutnet.cover$Taxon))
tern.pspathway$Species.name <- trimws(tolower(tern.pspathway$Species.name))

##find/checking for  matching species
#for(i in 1:nrow(nutnet.cover)) {
  # Debugging: Print current species name and if a match exists in B
 # print(paste("Checking:", nutnet.cover$Taxon[i]))
 # species_in_B <- tern.pspathway$Species.name == nutnet.cover$Taxon[i]
 # if(any(species_in_B)) {
 #   print("Match found")
#  } else {
 #   print("No match found")
 # }
#}

#for TERN Australia Database
## create a new column to track updated pathways
nutnet.cover$pspathway.updated <- FALSE
# Convert "NULL" strings to actual NA values for consistency, if necessary
nutnet.cover$ps_path[nutnet.cover$ps_path == "NULL"] <- NA
##looping through datsets and matching ps_pathways
# Iterate over each row in nutnet.cover
for(i in 1:nrow(nutnet.cover)) {
  # Check if pathway is NA (which now represents the original "NULL" values)
  if(is.na(nutnet.cover$ps_path[i])) {
    # Find the matching species in tern.pspathway dataset
    species_match_index <- which(tern.pspathway$Species.name == nutnet.cover$Taxon[i])
    
    # If a matching species is found and it has a pathway in tern.pspathway
    if(length(species_match_index) > 0 ) {
      
      # Print "Match found" to the console
      print("Match found")
      
      # Update the pathway in nutnet.cover with the pathway from tern.pspathway
      nutnet.cover$ps_path[i] <- 
        tern.pspathway$Photosynthetic.pathway..combined.[species_match_index]
      # Mark the row as updated... replaces false with true when a match is made. 
    nutnet.cover$pspathway.updated[i] <- TRUE
    }
  }
}

###checks for values that were updated. 
nutnet.cover[nutnet.cover$pspathway.updated == TRUE, ]


#For China trait database
#Transform column names for both dataframes
names(china.pspathway) <- tolower(trimws(names(china.pspathway)))
names(china.speciesname) <- tolower(trimws(names(china.speciesname)))
names(china.taxonomicname) <- tolower(trimws(names(china.taxonomicname)))
china.taxonomicname$species.id <- tolower(china.taxonomicname$species.id)


# Transform data within columns for china.pspathway
china.pspathway <- china.pspathway %>%
  mutate(across(everything(), ~tolower(trimws(.))))
china.pspathway$photo_path <- toupper(china.pspathway$photo_path)

# Transform data within columns for china.speciesname
china.speciesname <- china.speciesname %>%
  mutate(across(everything(), ~tolower(trimws(.))))

china_ps.pathway.combine <- china.pspathway %>%
  left_join(china.speciesname, by = "species.id") %>%
    select(-c(site.id, field.identified.genus, field.identified.species))

china_ps.pathway.combine <- left_join(china_ps.pathway.combine, 
            select(china.taxonomicname, species.id, family),  by = "species.id")

# Combine the 'genus' and 'species' columns to create a new 'taxon' column
china_ps.pathway.combine$taxon <- paste(china_ps.pathway.combine$accepted.genus,
                                      china_ps.pathway.combine$accepted.species)


##China database
# Loop through each row in nutnet
for (i in 1:nrow(nutnet.cover)) {
  # Check if photopath is NA
  if (is.na(nutnet.cover$ps_path[i])) {
    # Use which() to find indices of matching family in china
    #match_indices <- which(china_ps.pathway.combine$family == 
    #nutnet.cover$Family[i] & !is.na(china_ps.pathway.combine$photo_path))
   match_indices <- which(china_ps.pathway.combine$taxon == nutnet.cover$Taxon[i])
    # Check if there is at least one match and use the first one
    if (length(match_indices) > 0) {
      
      # Print "Match found" to the console
      print("Match found")
      
      # Update the photopath in nutnet from the first matching entry in china
      nutnet.cover$ps_path[i] <- china_ps.pathway.combine$photo_path[match_indices[1]]
# nutnet.cover$ps_path[i] <- tern.pspathway$Photosynthetic.pathway..combined.[species_match_index]
      # Optionally, mark as updated or perform some action
      # For example, setting an update flag (if you have one)
       nutnet.cover$pspathway.updated[i] <- TRUE
    }
  }
}



# remove non-living groups and non-control plots
nutnet.cover <- subset(nutnet.cover, functional_group != "NON-LIVE" & 
                    functional_group != "NULL" & functional_group != "LICHEN") 


###assign pathway
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus aggregatus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus retrorsus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus grayi"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus filiculmis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus filiculmis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus esculentus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus sp. (cedr.us)"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus eragrostis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus reflexus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus corymbosus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus schweinitzii"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus mollipes"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus plukenetii"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus difformis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus tenax"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus niveus var. leucocephalus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "eleocharis vivipara"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "leocharis sp."] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "fimbristylis autumnalis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "mollugo verticillata"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "allionia incarnata"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "portulaca amilis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "aristida sp. (rams.us)"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "sehima sp. (kidman.au)"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "andropogon sp."] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "kallstroemia angustifolia"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "blysmus sinocompressus"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "bolboschoenus robustus"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "cyperus niveus var. leucocephalu"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "pycreus filicinus"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "pycreus polystachyos"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$functional_group == "BRYOPHYTE"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$functional_group == "WOODY"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "ambrosia psilostachya"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "lysimachia europaea"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "solidago virgaurea"] <- "C3"
nutnet.cover$ps_path[nutnet.cover$Taxon == "bulbostylis"] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "bulbostylis sp."] <- "C4"
nutnet.cover$ps_path[nutnet.cover$Taxon == "bidens pilosas"] <- "C4"

nutnet.cover$ps_path[is.na(nutnet.cover$ps_path)] <- "C3"




#join the site lon,lat,country,elevation,experiment_type
nutnet.cover <- nutnet.cover %>%
  left_join(site_cordinates %>% 
              dplyr::select(-site_name), by = "site_code")

# removing unwanted columns
nutnet.cover <- nutnet.cover %>% 
  select(-'subplot', -year_trt, -first_year, -final_year, -number_experiment_years, 
         -slope, -aspect, -experiment_type)


# calculate total C4 percentage cover per year per site per plot
C4_nutnet.percent.cover <- nutnet.cover %>% 
  group_by(year, site_code, plot, block, longitude, latitude, elevation, trt) %>% 
  summarise("c4_percent" = (sum(max_cover[ps_path == "C4"], na.rm = TRUE) * 100 / 
                              sum(max_cover)))

# sites with 0% c4 is set to 100% c3
C4_nutnet.percent.cover$c3_percent <- 100 - C4_nutnet.percent.cover$c4_percent

#nutnet_percent.all <- nutnet.cover %>%
 # select(year, site_code, plot,  block, trt, Family, Taxon, functional_group, 
  #       local_lifeform,local_lifespan, longitude,latitude, elevation, ps_path, 
   #      max_cover) %>%
  #distinct() %>%
 # left_join(C4_nutnet.percent.cover, by = c("year", "site_code", "plot"))

























#C4_cover <- nutnet.cover %>% 
#  group_by(year, site_code, plot) %>% 
 # summarise("c4_percent" = sum(max_cover[ps_path == "C4"], na.rm = TRUE * 100 / sum(max_cover)))


 
#C4_nutnet_ <- nutnet.cover %>%
  # Ensure you filter rows before grouping if you only want to include specific conditions
 # filter(ps_path == "C4") %>%
  #group_by(year, site_code, plot) %>%
  # Use summarise() to calculate c4_percent
  #summarise(c4_percent = sum(max_cover, na.rm = TRUE) * 100 / 
    #          sum(nutnet.cover$max_cover[nutnet.cover$site_code == site_code & nutnet.cover$year == year & 
   #                                        nutnet.cover$plot == plot], na.rm = TRUE),
     #       .groups = 'drop') # This line controls the grouping structure of the output


#write.csv(nutnet.cover, "../c4percent.csv")
