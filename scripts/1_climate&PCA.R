########################################################################
# LIBRARIES                                                            #
########################################################################

library(tidyverse)

library(ncdf4)
library(raster)
library(parzer)

library(ggplot2)
library(patchwork)
library(leaflet)

library(factoextra)

#also run script "1a_coalesce_join_function.R" > function is used later on


########################################################################
# DATA                                                                 #
########################################################################

# POPULATION & CULTIVAR LOCATIONS
pops <- read.csv("./data/climate/populations_locations.csv", header = FALSE, stringsAsFactors=FALSE)
str(pops)

colnames(pops) <- pops[1,]
pops <- pops[-1,]
str(pops)

## Keep only populations & cultivars that are used in the analysis
pc.keep <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", # L. bienne Spain/Med
             "10", "11", "12", "13", "14", "15", "CGa1",  "19", 
             "BH", "Bro", "CR", "IOW1", "IOW2", # L. bienne Atlantic
             "Lil", "LJLb1", "Lla", "Man", "Mat", "Roc", 
             "Saf", "Sut", "Tal", "Tor", "Tym", "Vil", "Dor", 
             "Ara", "Ari", "Ble", "Bol", "Ede", "Gis", "Lir", # L. usitatissimum
             "Mar", "Mon", "Olg", "Ome", "Pri", "Rab", "Suz", 
             "Tin", "Vol")

pops <- pops %>% filter(Pop %in% pc.keep) %>% droplevels() %>% as.data.frame()
str(pops)

## calculate decimal coordinates and check they are correct
pops$Lat_dms <- gsub("\"", "s", as.character(pops$Lat_dms))
pops$Lon_dms <- gsub("\"", "s", as.character(pops$Lon_dms))

pops$Lat_dms <- gsub("d|m|s", " ", as.character(pops$Lat_dms))
pops$Lon_dms <- gsub("d|m|s", " ", as.character(pops$Lon_dms))

coords_correct <- data.frame(
  x = parse_lon(pops$Lon_dms),
  y = parse_lat(pops$Lat_dms),
  ID = pops$Pop)

str(coords_correct)

leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng = coords_correct$x, lat = coords_correct$y, popup = coords_correct$ID)

## add decimal coords to l. bienne populations
pops$Lat <- ifelse(is.na(pops$Lat) & match(pops$Pop, coords_correct$ID), coords_correct$y, pops$Lat)
pops$Lon <- ifelse(is.na(pops$Lon) & match(pops$Pop, coords_correct$ID), coords_correct$x, pops$Lon)

## check again everything is ok with map
sppCol <- colorFactor(palette = 'RdYlGn', pops$spp)

leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addCircleMarkers(
    data = pops,
    lng = ~as.numeric(Lon),
    lat = ~as.numeric(Lat),
    popup = ~Pop,
    color = ~sppCol(spp))

## turn locations into spatial dataframe
pops$Lat <- as.numeric(pops$Lat)
pops$Lon <- as.numeric(pops$Lon)

coordinates(pops) = ~Lon+Lat
proj4string(pops) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# CLIMATE

common_source <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/"

vars <- c("wc2.1_30s_tmin.zip",
          "wc2.1_30s_tmax.zip",
          "wc2.1_30s_tavg.zip",
          "wc2.1_30s_prec.zip",
          "wc2.1_30s_srad.zip",
          "wc2.1_30s_wind.zip",
          "wc2.1_30s_vapr.zip")

common_dest <- paste0(getwd(), "/data/climate/")

# ## download climate data (WorldClim 2.1, historical, 30sec - see https://www.worldclim.org/data/worldclim21.html)
# ### this section is commented out because it has to be run only once!
# for(i in 1:length(vars)){
# download.file(paste0(common_source, vars[i]), paste0(common_dest, vars[i]))
#   }
# 
# ### unzip
# for(i in 1:length(vars)){
#   unzip(paste0(common_dest, vars[i]), 
#         overwrite = FALSE, 
#         exdir = paste0(common_dest, str_remove(vars[i], ".zip")))
#   }
# 
# ### create stack for each variable and save to list, then change layers names to months names (1 = Jan, 12 = Dec)
# files.list <- NULL #make list of file names for each variable
# 
# for(i in 1:length(vars)){
#   files.list[[i]] <- list.files(paste0(common_dest, str_remove(vars[i], ".zip")), full.names = TRUE)
#   names(files.list)[i] <- str_remove(str_remove(vars[i], ".zip"), "wc2.1_30s_")
#   }
# 
# layers.names <- as.character(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# 
# stacks.list <- list()
# 
# for(i in 1:length(vars)){
#   stk <- stack(files.list[[i]])
#   stacks.list[[i]] <- stk
#   names(stacks.list)[i] <- names(files.list)[i]
#   names(stacks.list[[i]]) <- layers.names
#   }
# 
# ### crop stacks and save files
# ext <- extent(-25, 75, 25, 71) #ca. L. bienne distribution
# 
# cropped.list <- list()
# 
# for(i in 1:length(vars)){
#   cropped.list[[i]] <- stack(crop(stacks.list[[i]], ext))
#   names(cropped.list)[i] <- names(files.list)[i]
#   }
# 
# writeRaster(cropped.list[[7]], #here I save raster as layers from 1 to 7 - change number to select rasters
#             filename = paste0(common_dest, paste0(names(files.list)[7], ".tif")), 
#             options = "INTERLEAVE=BAND", 
#             overwrite = FALSE)

## import climate data (can be run only if data have already been downloaded)
### remember to unzip climate data folder before running the following
common_dest <- paste0(getwd(), "/data/climate/climate_data/") #update path
ext <- extent(-25, 75, 25, 71) #ca. L. bienne distribution
layers.names <- as.character(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

files.list <- list()

for (i in 1:length(vars)){
  nm <- paste0(paste0(common_dest, str_remove(str_remove(vars[i], ".zip"), "wc2.1_30s_")), ".tif")
  files.list[[i]] <- nm
  names(files.list)[i] <- str_remove(str_remove(vars[i], ".zip"), "wc2.1_30s_")
  }

stacks.list <- list()

for(i in 1:length(vars)){
  stk <- stack(files.list[[i]])
  stacks.list[[i]] <- stk
  names(stacks.list)[i] <- names(files.list)[i]
  names(stacks.list[[i]]) <- layers.names
  }

### extract variables
rasterOptions(maxmemory = 1e+10) #check this is ok with your computer
climdf.list <- NULL

for(i in 1:length(vars)){
  climdf.list[[i]] <- raster::extract(stacks.list[[i]], pops, method = "simple", cellnumbers = TRUE, df = TRUE)
  climdf.list[[i]]$Pop <- pops$Pop
  names(climdf.list)[i] <- names(files.list)[i]
  }

climdf.list
#### some locations are very close to coast and I get NA when extracting raster values.
#### this is because some raster cells for points close to the sea, in fact represent the sea

### I can find closest non-NA cell to those points and extract values
### pick any list item (any clim. variable df in the list) and extract cell numbers where rows are all NA

na.cells.nr <- climdf.list[[1]]$cells[rowSums(is.na(climdf.list[[1]])) == 12] #one cultivar (Mar, NA = 13), will stay NA because clim. raster covers europe, and Mar is from USA, but I don't need clim. data for cultivars
na.pops.row <- climdf.list[[1]]$ID[rowSums(is.na(climdf.list[[1]])) == 12]

na.cells <- data.frame(xyFromCell(stacks.list[[1]]@layers[[1]], na.cells.nr, spatial = FALSE)) #find xy coordinates for NA cells
na.cells$Pop <- pops$Pop[na.pops.row]  

coordinates(na.cells) = ~x+y
proj4string(na.cells)=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #create spatial object

climdf.list.na <- NULL

for(i in 1:length(vars)){
  climdf.list.na[[i]] <- raster::extract(stacks.list[[i]], 
                                         na.cells, 
                                         buffer = 3000, fun = mean, na.rm = TRUE,
                                         df = TRUE)
  climdf.list.na[[i]]$Pop <- na.cells$Pop
  names(climdf.list.na)[i] <- names(files.list)[i]
  }

### fill climate dataframes NAs with values found using 3km buffer
for(i in 1:length(vars)){
  climdf.list[[i]] <- coalesce_join(mutate_all(climdf.list[[i]], as.character), mutate_all(climdf.list.na[[i]], as.character), by = 'Pop')
  }

### create unique climate dataframe for all clim. variables
climdf <- bind_rows(climdf.list, .id = "column_label")

str(climdf)

colnames(climdf)[1] <- "clim_var"
climdf <- climdf[, -2]
climdf <- climdf[, c(15, 1:14)]

str(climdf)

climdf <- climdf %>% #long format
  mutate_at(layers.names, as.numeric) %>%
  gather(month, value, Jan:Dec)

str(climdf)

climdf <- climdf %>% #summarise climate data by season
  mutate(season = case_when(
    month %in% c("Dec", "Jan", "Feb") ~ "DJF",
    month %in% c("Mar", "Apr", "May") ~ "MAM",
    month %in% c("Jun", "Jul", "Aug") ~ "JJA",
    month %in% c("Sep", "Oct", "Nov") ~ "SON"
    )) %>% #View()
  group_by(Pop, clim_var, season) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() #%>% View()

# write.csv(climdf, "./data/climate/seasonalClim_populations.csv")


# ADDITIONAL GEOGRAPHICAL VARIABLES

## Import elevation tiff (30s resolution, downloaded directly from worldclim website)
# elev <- stack(paste0(common_dest, "wc2.1_30s_elev.tif"))
# names(elev) <- "elevation"
# 
# elev.cropped <- crop(elev, ext)
# writeRaster(elev.cropped, 
#             filename = paste0(common_dest, "elev.tif"), 
#             options = "INTERLEAVE=BAND", 
#             overwrite = FALSE)

elev <- stack(paste0(common_dest, "elev.tif"))
pops$elev_wc <- raster::extract(elev, pops)

## Import distance from coast raster, created from E-obs data
dst <- stack(paste0(common_dest, "distance.from.coast.grd"))
names(dst) <- "Distance.from.Coast"
pops$dCst <- raster::extract(dst, pops)

## Save population locations with additional geographical variables
str(pops)
pops_df <- as.data.frame(pops)
names(pops_df)[14] <- "elev_wc"
str(pops_df)

# write.csv(pops_df, paste0(common_dest, "populations_locations_1.csv"))



########################################################################
# DESCRIPTION OF POPULATIONS CLIMATE WITH PCA APPROACH                 #
########################################################################

# Prepare input data frame
pop.clim <- climdf %>%
  pivot_wider(id_cols = Pop,
              names_from = c(clim_var, season),
              values_from = value) %>%
  left_join(pops_df %>%
              mutate(Alt = if_else(is.na(as.numeric(Alt)), elev_wc, as.numeric(Alt))) %>%
              dplyr::select(c(spp, Pop, Country, Lat, Lon, Alt, Distance.from.Coast)) %>%
              as.data.frame(),
            .,
            by = c("Pop")) %>%
  filter(spp == "bienne") %>% # !!!NB: PCA is done for L. bienne only, not for L. usitatissimum
  as.data.frame() #%>% str()

# Run PCA
clim.pca <- prcomp(pop.clim[,c(8:35)], center = TRUE, scale. = TRUE) #8:35 all climatic variables per season

screeplot(clim.pca, type = "l", npcs = 12, main = "Screeplot of the first 12 PCs")
abline(h = 1, col="red", lty = 5)
legend("topright", legend = c("Eigenvalue = 1"),
       col = c("red"), lty = 5, cex = 0.6)

cumpro <- cumsum(clim.pca$sdev^2 / sum(clim.pca$sdev^2))
plot(cumpro[0:12], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 3, col="blue", lty = 5)
abline(h = 0.9, col = "blue", lty = 5)
legend("topleft", legend = c("Cut-off @ PC3"),
       col = c("blue"), lty = 3, cex = 0.6)

pc12.clim <- fviz_pca_biplot(clim.pca,
                             axes = c(1, 2),
                             col.ind = pop.clim$Lat, # Color by the quality of representation
                             fill.ind = pop.clim$Lat,
                             geom = "point",
                             pointsize = 3,
                             stroke = 1.5,
                             labelsize = 3,
                             alpha = 0.5,
                             #addEllipses = TRUE, # Concentration ellipses
                             #ellipse.type = "confidence",
                             legend.title = "Latitude",
                             #palette = pal.pc[c(1:12)], #pal.pc[c(3,10)],
                             col.var = "gray42",
                             label = c("var"),
                             repel = TRUE,
                             title = NULL) +     # Avoid text overlapping
  scale_x_continuous(position = "top") +
  theme_bw(base_size = 15)

pc13.clim <- fviz_pca_biplot(clim.pca,
                             axes = c(1, 3),
                             col.ind = pop.clim$Lat, # Color by the quality of representation
                             fill.ind = pop.clim$Lat,
                             geom = "point",
                             pointsize = 3,
                             stroke = 1.5,
                             labelsize = 3,
                             alpha = 0.5,
                             #addEllipses = TRUE, # Concentration ellipses
                             #ellipse.type = "confidence",
                             legend.title = "Latitude",
                             #palette = pal.pc[c(1:12)], #pal.pc[c(3,10)],
                             col.var = "gray42",
                             label = c("var"),
                             repel = TRUE,
                             title = NULL) +     # Avoid text overlapping
  scale_x_continuous(position = "top") +
  theme_bw(base_size = 15)


# Get variables contributions to PC 1-3 and plot
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
  }

loadings <- clim.pca$rotation
sdev <- clim.pca$sdev
var.coord <- t(apply(loadings, MARGIN = 1, var_coord_func, sdev)) 
head(var.coord[, 1:3])

var.cos2 <- var.coord^2
head(var.cos2[, 1:3])

comp.cos2 <- apply(var.cos2, MARGIN = 2, FUN = sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- as.data.frame(t(apply(var.cos2, MARGIN = 1, contrib, comp.cos2)))[1:3]

var.contrib$vars <- rownames(var.contrib)
rownames(var.contrib) <- seq(1:nrow(var.contrib))

var.clim <- var.contrib %>%
  gather(PC, value, PC1:PC3) %>% #str()
  ggplot(., aes(x = reorder(factor(vars), -value), y = value, fill = PC)) +
  geom_col(width = 0.85, position = position_dodge(), col = "black") +
  #scale_fill_manual(values = c("#002261", "darksalmon", "#F4E2D4", "#F6FBEE", "ivory3")) +
  scale_fill_grey(start = 0, end = 1) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(position = "right") +
  xlab("Climatic Variables") + ylab("Contribution %") +
  theme_bw(base_size = 15) +
  #theme(legend.position = "none") +
  coord_flip()

# Make and save Plots & PCs file
p <- pc12.clim + pc13.clim + var.clim  +
  plot_layout(guides = 'collect')

# ggsave("./output/clim_pca.pdf", p, width = 24, height = 6, units = "in")

pop.clim <- pop.clim %>%
  bind_cols(.,
            as.data.frame(clim.pca$x[,1:3])) %>%
  dplyr::select(c(Pop, PC1:PC3)) %>%
  as.data.frame() #%>% View()

# write.csv(pop.clim, "./data/climate/climPCs.csv")
