

library(raster)

####################################
# 1. INPUT ENVIRONMENT LAYERS
env <- stack("/Users/au529793/Projects/Thesis/DATA/occurrences/GIS DATA/rasters/PR_env_layers.grd")
env <- env[[c(1,2,5,7,8)]]


####################################
# 2 : GET RECORDS
focsp <- read.csv("/Users/au529793/Desktop/Hyraulics_Focal_Species.csv")

pts <- readRDS('/Users/au529793/Projects/Postdoc/CLCC/DATA/occurrences/PR_PLANT_OCCS_8.3.16.RDS')

pts$scrubbed.name[pts$scrubbed.name=='Prestoea acuminata'] <- 'Prestoea montana'
pts$scrubbed.name[pts$scrubbed.name=='Myrcia amazonica'] <- 'Myrcia leptoclada'
pts$scrubbed.name[pts$scrubbed.name=='Nectandra turbacensis'] <- 'Ocotea sintenisii'

pts <- unique(pts)

focpts <- droplevels(pts[pts$scrubbed.name %in% 
                           c(as.character(focsp$NAME), as.character(focsp$NAME2)),])

focpts$PPT <- extract(exp(env[[2]]), focpts[,c('LON','LAT')])

splist <- split(focpts, focpts$scrubbed.name)

df <- data.frame()
for(i in 1:length(splist)){
  path <- "/Users/au529793/Projects/Thesis/DATA/Chapter 4/allAICc_unscaled_mods/_grids/"
  
  if(file.exists(paste0(path, focsp$CODE[i], '.grd'))){
    r <- raster(paste0(path, focsp$CODE[i], '.grd'))
    opt.ppt <- cellStats((r == cellStats(r, max)) * exp(env$log.annual.ppt), max)
  } else {
    opt.ppt <- NA
  } 
  df <- rbind(df, cbind(as.character(focsp$CODE[i]), 
                        names(splist)[i], 
                        median(splist[[i]]$PPT, na.rm=T), 
                        min(splist[[i]]$PPT, na.rm=T), 
                        max(splist[[i]]$PPT, na.rm=T), 
                        opt.ppt,
                        quantile(splist[[i]]$PPT, 0.05, na.rm=T),
                        quantile(splist[[i]]$PPT, 0.95, na.rm=T)))
}

for(j in 3:8) { df[,j] <- round(as.numeric(as.character(df[,j])), 0) }

names(df) <- c('SPCODE', 'SpeciesName', 'MedianPPT', 'MinPPT', 'MaxPPT',
               'OptimalPPT','5%quantile','95%quantile')

rownames(df) <- NULL

df


# write.csv(df, '/Users/au529793/Desktop/PPT_at_Occurrences.csv', row.names = F)


####################################
## GBIF
####################################

library(rgbif)
library(CoordinateCleaner)
occ_gbif <- occ_search(scientificName=as.character(focsp$NAME), 
                     limit=5000,
                     hasCoordinate=TRUE)

# Now we clean the coordinates for each species in turn
occ_gbif_cc <- list()
for(i in 1:length(occ_gbif)){
  tmpdata <- data.frame(species=focsp$NAME[i],
                        decimallongitude=occ_gbif[[i]]$data$decimalLongitude,
                        decimallatitude=occ_gbif[[i]]$data$decimalLatitude)
  
  # Run the coordinate cleaner function
  occ_gbif_cc[[i]] <- clean_coordinates(tmpdata, value="clean")
}

# Clean up and make them a single data frame
occ_gbif_cc_clean <- do.call(rbind, occ_gbif_cc)


# Read the gridded climate data just downloaded
path <- "/Users/au529793/Projects/Manuscripts/Global Change Biology/Uriarte Seedlings/CODE/chelsa/"
env <- stack(paste0(path, list.files(path)[1:19]))

# Extract the climate data at occurrences
occ_gbif_cc_clean <- cbind(occ_gbif_cc_clean, extract(env, occ_gbif_cc_clean[,2:3]))

for(i in which(grepl("CHELSA", names(occ_gbif_cc_clean)))){
  tmp <- round(tapply(occ_gbif_cc_clean[,i], occ_gbif_cc_clean$species, mean, na.rm=T), 0)
  df[,gsub("CHELSA_", "", names(occ_gbif_cc_clean)[i])] <- tmp[match(df$SpeciesName, names(tmp))]
  
  tmp95 <- round(tapply(occ_gbif_cc_clean[,i], 
                        occ_gbif_cc_clean$species, quantile, 0.95, na.rm=T), 0)
  df[,paste0(gsub("CHELSA_", "", names(occ_gbif_cc_clean)[i]), "_95%")] <- 
    tmp[match(df$SpeciesName, names(tmp))]
  
  tmp05 <- round(tapply(occ_gbif_cc_clean[,i], 
                        occ_gbif_cc_clean$species, quantile, 0.05, na.rm=T), 0)
  df[,paste0(gsub("CHELSA_", "", names(occ_gbif_cc_clean)[i]), "_05%")] <- 
    tmp[match(df$SpeciesName, names(tmp))]
}

df



write.csv(df, '/Users/au529793/Desktop/PPT_at_Occurrences.csv', row.names = F)





####################################
## R BIEN
####################################


# Get occurrences from BIEN
library(BIEN)
occ <- BIEN_occurrence_species(as.character(focsp$NAME),
                               only.new.world=T)
occ$species <- occ$scrubbed_species_binomial
occ$decimallongitude <- occ$longitude
occ$decimallatitude <- occ$latitude

# Remove rows with no coordinates
occ <- occ[!is.na(occ$latitude),]

# Run the coordinate cleaner function
occ_cc <- CoordinateCleaner::clean_coordinates(occ, value="clean")

### Get the Chelsa climate data (bioclimatic variables)
# library(chelsaDL)
# library(stringr)
# ch_data <- ch_queries(variables = c("bio"),
#                       layers = 18:19,
#                       timeframes = "1979-2013")
# 
# ### Read the gridded climate data just downloaded
# temp <- stack(ch_dl(ch_data)$path)

path <- "/Users/au529793/Projects/Manuscripts/Global Change Biology/Uriarte Seedlings/CODE/chelsa/"
env <- stack(paste0(path, list.files(path)[1:12]))

# Extract the climate data at occurrences
occ_cc <- cbind(occ_cc, extract(env, occ_cc[,3:2]))

for(i in which(grepl("CHELSA",names(occ_cc)))){
  tmp <- round(tapply(occ_cc[,i], occ_cc$species, mean, na.rm=T), 0)
  df[,gsub("CHELSA_", "", names(occ_cc)[i])] <- tmp[match(df$SpeciesName, names(tmp))]
}

