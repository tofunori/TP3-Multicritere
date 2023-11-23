#Welcome! TP3_PIF By Thierry Laurent
# dont forget to change the path of your files

# Step 1:

# load libraries
library(terra)
library(tmap)
library(tmaptools)
library(sf)
library(dplyr)
library(randomForest)
library(ggplot2)
library(RColorBrewer)
library(yardstick)
library(raster)


# RcolorBrewer colors. Many choices!

display.brewer.all()


################################################################################
# First load the shp and reprojected them with with the desired EPSG. In this case is EPSG: 2951

# Load the raster
raster_file <- "canopy.tif"
raster_data <- raster(raster_file)

#
shapefile_paths <- c("shp/Batiments.shp", 
                     "shp/Lacs.shp", 
                     "Routes_chemins_repro.shp",
                     "shp/Tourbieres.shp",
                     "shp/Rivieres.shp")

desired_crs <- "EPSG:2951"  # Choose the desired projection

shapefiles <- lapply(shapefile_paths, function(shp) {
  shp_data <- st_read(shp)
  st_transform(shp_data, crs = desired_crs)
})

# Export the shp to the desired folder
output_paths <- c("Batiments_repro2.shp", 
                  "Lacs_repro2.shp", 
                  "Routes_chemins_repro2.shp",
                  "Tourbieres_repro2.shp",
                  "Rivieres_repro2.shp")

Map(function(shp, out_path) {
  st_write(shp, out_path)
}, shapefiles, output_paths)


# Here its optional. First you need to extract the non paved roads from the shp

shp <- st_read("Routes_repro.shp")
print(shp)

head(shp)

filtered_shp <- shp[shp$DESCRIPTIO == "Chemin carrossable non pavé", ]

#Export it
st_write(filtered_shp, "routes_chemins_repro.shp")

################################################################################

#This process will create a buffer for each shapefile, then rasterize and reclassify them accordingly

################################################################################


# Load the reference raster that serve as a extend
ref_raster <- raster("canopy.tif")

# Load the reprojected shp
shapefiles <- c("Batiments_repro2.shp", 
                "Lacs_repro2.shp", 
                "Routes_chemins_repro2.shp",
                "Tourbieres_repro2.shp",
                "Rivieres_repro2.shp")

# Buffer 
buffer_distances <- c(100, 20, 20, 40, 10)  # Modify these values as needed

buffers <- mapply(function(shp, dist) {
  shp_data <- st_read(shp)
  st_buffer(shp_data, dist = dist)
}, shp=shapefiles, dist=buffer_distances, SIMPLIFY=FALSE)

# Rasterize each buffer using the extent of the reference raster and save each raster
rasters <- lapply(seq_along(buffers), function(i) {
  rast <- raster(ref_raster)
  rasterized <- rasterize(buffers[[i]], rast, field=1, background=1, fun=function(x, ...) 0) #here you need to determine the values of the buffers and the background
  # Save the rasterized buffer
  output_filename <- paste0("reclass_", i, ".tif")
  writeRaster(rasterized, output_filename, format="GTiff", overwrite=TRUE)
  return(rasterized)
})


raster_paths <- c("batiments_reclass.tif", 
                  "lacs_reclass.tif", 
                  "routes_chemins_reclass.tif",
                  "tourbieres_reclass.tif",
                  "rivieres_reclass.tif")

raster_stack <- stack(raster_paths)

# Multiply the 5 raster together 

rasters_sum <- calc(raster_stack, fun=function(x) prod(x, na.rm=TRUE))

#Optional plot

#plot(rasters_sum)

# Export the binary raster of the sum of the 5 rasters 
writeRaster(rasters_sum, "rasters_sum.tif", format="GTiff", overwrite=TRUE)


#################################################################################

# Part 2: Now its time to produce categorital rasters to give certain weight for a range of value.
#The Eucledean distance fonction is use for the routes raster

################################################################################

# LOAD the Raster # Now Here we have to put a raster to calculate the distance. 
# In this case I use the DEM from the same source of the canopy and slope one


raster_file <- "dem.tif"  
r <- rast(raster_file)

raster_file <- "routes_eucle.tif"  
routes_eucle <- rast(raster_file)

#Eucledean distance

routes <- vect("routes_chemins_repro.shp")
routes2 <- st_read("routes_chemins_repro.shp")

dist <- distance(r, routes)

tm_shape(dist)+
  tm_raster(style= "cont")+
  tm_shape(routes2)+
  tm_lines(col="black")+
  tm_layout(legend.outside=TRUE)


#######Classify the eucledeane distance raster

# Load the routes eucledean raster
raster <- rast("routes_eucle.tif")

# Print the minimum and maximum values 
min_max_values <- range(raster, na.rm = TRUE)
print(min_max_values)

# Choose the classification 
reclass_func <- function(x) {
  ifelse(x >= 0 & x <= 500, 10,
         ifelse(x > 500 & x <= 1500, 8,
                ifelse(x > 1500 & x <= 2500,6,
                       ifelse(x > 2500 & x <= 3500, 4,
                              ifelse(x > 3500 & x <= 6000, 2, NA)
                       )
                )
         )
  )
}

reclassified_raster <- app(raster, reclass_func)


# Choose the color palette using ColorBrewer library
colors <- brewer.pal(5, "Set1")  # Choose a palette with 5 colors
plot(reclassified_raster, col=colors, main="Reclassified Raster")


# Save the raster
output_file <- "routes_eucle_reclass2.tif"
writeRaster(reclassified_raster, filename=output_file, overwrite=TRUE)

# Create frequency table
attribute_table <- freq(reclassified_raster)
print(attribute_table)


#########################

# Now with the slope and canopy raster, we will aggregate, low filter and resample them before doing a reclassification

#### Raster Canopy #######

#Load raster
file_path <- "canopy.tif"
raster <- rast(file_path)


#Aggregate the raster. Choose the factor you want
aggregated_raster <- aggregate(raster, fact=10, fun=max)

#Low pass filter
filtered_raster <- focal(aggregated_raster, w=matrix(1, 3, 3), fun=max)  # Better to you the max for the canopy

#Resample the raster
resampled_raster <- resample(filtered_raster, raster, method="bilinear")

# Print the minimum and maximum 
min_max_values <- range(resampled_raster, na.rm = TRUE)
print(min_max_values)


# 4. Reclassify using your custom function
reclass_func <- function(x) {
  ifelse(x >= 0 & x <= 3, 0,
         ifelse(x > 3 & x <= 10, 4,
                ifelse(x > 10 & x <= 20, 8,
                       ifelse(x > 20 & x <= 25, 12,
                              ifelse(x > 25 & x <= 30, 16,
                                     ifelse(x > 30 & x <= 41, 20, NA)
                              )
                       )
                )
         )
  )
}
reclassified_raster <- app(resampled_raster, reclass_func)

# View the final raster
plot(reclassified_raster)

# Save the canopy reclass raster
output_file <- "canopy_reclass2.tif"  
writeRaster(reclassified_raster, filename=output_file, overwrite=TRUE)


##### Raster Slope ##########################################################

# Load the raster
file_path <- "pentes.tif"
raster <- rast(file_path)

#Aggregate the raster 
aggregated_raster <- aggregate(raster, fact=10, fun=mean)

#Apply a low-pass filter 
filtered_raster <- focal(aggregated_raster, w=matrix(1, 3, 3), fun=mean)

#Resample the raster 
resampled_raster <- resample(filtered_raster, raster, method="bilinear")

# Print the minimum and maximum 
min_max_values <- range(resampled_raster, na.rm = TRUE)
print(min_max_values)

# 4. Reclassify the slope raster
reclass_func <- function(x) {
  ifelse(x >= 0 & x <= 10, 10,
         ifelse(x > 10 & x <= 20, 8,
                ifelse(x > 20 & x <= 30, 6,
                       ifelse(x > 30 & x <= 40, 4,
                              ifelse(x > 40 & x <= 50, 2,
                                     ifelse(x > 50 & x <= 103, 0, NA)
                              )
                       )
                )
         )
  )
}
reclassified_raster <- app(resampled_raster, reclass_func)

# View the final raster
plot(reclassified_raster)

#Save the raster
output_file <- "pentes_reclass.tif"  
writeRaster(reclassified_raster, filename=output_file, overwrite=TRUE)

############################## Final step ######################################

#Raster algebra :) Do the sum of the canopy, slope and routes rasters and then multiply by the raster sum binaire

# Load the all the rasters
raster1 <- rast("routes_eucle_reclass2.tif")
raster2 <- rast("pentes_reclass.tif")
raster3 <- rast("canopy_reclass2.tif")
raster4 <- rast("rasters_sum.tif")
raster5 <- rast("rivieres_reclass.tif")
raster6 <- rast("routes_reclass.tif")
raster7 <- rast("lacs_reclass.tif")
raster8 <- rast("marais_reclass.tif")
raster9 <- rast("batiments_reclass.tif")


# Align the rasters
aligned_raster2 <- resample(raster2, raster1, method="bilinear")
aligned_raster3 <- resample(raster3, raster1, method="bilinear")
aligned_raster4 <- resample(raster4, raster1, method="bilinear")

# Sum the first three rasters
sum_raster <- raster1 + aligned_raster2 + aligned_raster3


# Multiply the sum with the fourth raster
result_raster <- sum_raster * aligned_raster4

# View the result
plot(result_raster)

# Save the final raster!
output_file <- "raster_final.tif"  
writeRaster(result_raster, filename=output_file, overwrite=TRUE)

########## Here you can reclassify in 5 classes the final raster  ##############

# Load the final raster
raster_final <- "raster_final.tif"  
r <- rast(raster_final)

# Print the minimum and maximum 
min_max_values <- range(r, na.rm = TRUE)
print(min_max_values)
plot(r)

# Load the raster
raster_final <- "raster_final.tif"  
r <- rast(raster_final)
reclass_func <- function(x) {
  ifelse(x == 0, NA,  # Set 0 to NA
         ifelse(x >= 1 & x <= 15, 1,  # Classify 1-15 as 1
                ifelse(x > 15 & x <= 20, 2,  # 16-20 as 2
                       ifelse(x > 20 & x <= 25, 3,  # 21-25 as 3
                              ifelse(x > 25 & x <= 30, 4,  # 26-30 as 4
                                     ifelse(x > 30 & x <= 40, 5, NA)  # 31-40 as 5, others as NA
                              )
                       )
                )
         )
  )
}


# Apply the manual reclassification
raster_final_reclass <- app(r, manual_reclass_func)


# Define colors for the plot
num_classes <- 5
color_palette <- brewer.pal(num_classes, "RdYlGn") # Replace "YlGnBu" with your chosen palette

# Plot the reclassified raster
plot(raster_final_reclass, col=color_palette, main="Potentiel de coupe forestiere")

### Create a histogram and plot it ####

# Extract the values 
values <- values(raster_final_reclass)

# Calculate the area in km² 
freq_data <- freq(raster_final_reclass)
freq_data$area_km2 <- freq_data$count / 1e6  # Convert area to km²

# Plot the histogram
data_for_plot <- data.frame(category = freq_data[, "value"], area_km2 = freq_data[, "area_km2"])

ggplot(data_for_plot, aes(x = category, y = area_km2, fill = as.factor(category))) +
  geom_bar(stat = "identity", color="black") +
  geom_text(aes(label = round(area_km2, 2)), vjust = -0.5, color = "black") +  
  scale_fill_brewer(palette = "Set3", name = "Category") +  # Custom legend title
  theme_light() +
  labs(title = "Area per Category in km²",
       x = "Category",
       y = "Area (km²)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Save the final raster! Good job ;)
output_file <- "raster_final_5class.tif"  
writeRaster(raster_final_reclass, filename=output_file, overwrite=TRUE)


###### Plot all the rasters in one plot ###########

# Plotting setup: 2 rows, 2 columns
par(mfrow = c(3, 3))

# Plot each raster
plot(raster1, main = "Routes_eucle")
plot(raster2, main = "Pentes")
plot(raster3, main = "Canopy")
plot(raster4, main = "Rasters_binaire")
plot(raster5, main = "Rivieres")
plot(raster6, main = "Routes")
plot(raster7, main = "Lacs")
plot(raster8, main = "Marais")
plot(raster9, main = "Batiments")

# Reset plotting layout
par(mfrow = c(1, 1))



