# Masterarbeit Alexander Julian Koreneef
# Moderne Fernerkundungsmethoden zur Ermittlung der Vitalität von Baumbeständen wiedervernässter Feuchtwälder
# M.Sc. Geographie VT: Stadt- und Landschaftsökologie, Ruhr-Universität Bochum
# R-Script

# === 0 Packages ===============================================================================================================
library(sp)
library(raster)
library(caret)
library(sf)
library(rgdal)
library(randomForest)
library(RStoolbox)

# === 1 Einladen der Daten (Import/stack) ======================================================================================
# Datengrundlage
#setwd("C:/Users/Alexander/Documents/Uni/M.Sc.Landschaftsökologie/2019_SS_Masterarbeit/IV_R/RStudio")

# ===== 1.1) UAV 2016 ==========================================================================================================
AU3_2016_G <- raster('data/Area1MSUtm_201608261_proj_clip.tif', band=1) # G
AU3_2016_R <- raster('data/Area1MSUtm_201608261_proj_clip.tif', band=2) # R
AU3_2016_N1 <- raster('data/Area1MSUtm_201608261_proj_clip.tif', band=3) # NIR1
AU3_2016_N2 <- raster('data/Area1MSUtm_201608261_proj_clip.tif', band=4) # NIR2
AU3_2016_N3 <- raster('data/Area1MSUtm_201608261_proj_clip.tif', band=5) # NIR3

AU3_2016 <- stack(list(AU3_2016_G, AU3_2016_R, AU3_2016_N1, AU3_2016_N2, AU3_2016_N3)) 

# ===== 1.2) UAV 2019 ==========================================================================================================
AU3_2019 <- stack("data/AU3_2019_proj_clip.tif")

# ===== 1.3) UG/AOI  ===========================================================================================================
eUG <- readOGR("data/eUG.shp")
AU3 <- readOGR("data/AU3.shp")
UG <- readOGR("data/AU3_Kreise.shp")

# ===== 1.4) Trainingsdaten/Groundtrouth =======================================================================================
trainData2016 <- st_read("data/AU3_2016_Baumkronen_join_2020.shp")
train2016 <- as_Spatial(trainData2016)
trainData2019 <- st_read("data/AU3_2019_Baumkronen_join.shp")
train2019 <- as_Spatial(trainData2019)

# ===== 1.5) Sentinel 2016 =====================================================================================================
S_2016_B <- raster('data/T32ULC_20160826T104022_B02.jp2') # Blau
S_2016_G <- raster('data/T32ULC_20160826T104022_B03.jp2') # Grün
S_2016_R <- raster('data/T32ULC_20160826T104022_B04.jp2') # Rot
S_2016_N <- raster('data/T32ULC_20160826T104022_B08.jp2') # NIR

S_2016_s <- stack(list(S_2016_B, S_2016_G, S_2016_R, S_2016_N))

# ===== 1.6) Sentinel 2019 =====
S_2019a_B <- raster('data/T32ULC_20190826T104029_B02_10m.jp2') # Blau
S_2019a_G <- raster('data/T32ULC_20190826T104029_B03_10m.jp2') # Grün
S_2019a_R <- raster('data/T32ULC_20190826T104029_B04_10m.jp2') # Rot
S_2019a_N <- raster('data/T32ULC_20190826T104029_B08_10m.jp2') # NIR

S_2019a_s <- stack(list(S_2019a_B, S_2019a_G, S_2019a_R, S_2019a_N))

S_2019b_B <- raster('data/T32ULC_20190915T104019_B02_10m.jp2') # Blau
S_2019b_G <- raster('data/T32ULC_20190915T104019_B03_10m.jp2') # Grün
S_2019b_R <- raster('data/T32ULC_20190915T104019_B04_10m.jp2') # Rot
S_2019b_N <- raster('data/T32ULC_20190915T104019_B08_10m.jp2') # NIR

S_2019b_s <- stack(list(S_2019b_B, S_2019b_G, S_2019b_R, S_2019b_N))

# ===== 1.7) Mask/Crop =====
AU3_2019_c <- crop(AU3_2019, AU3)
AU3_2016_c <- crop(AU3_2016, AU3)

AU3_2019_m <- mask(AU3_2019_c, AU3)
AU3_2016_m <- mask(AU3_2016_c, AU3)

S_2016_c <- crop(S_2016_s, eUG)
S_2019a_c <- crop(S_2019a_s, eUG)
S_2019b_c <- crop(S_2019b_s, eUG)

# === 2 Aufbereitung (Projektion/Normalisieren) ====
# ===== 2.1 Koordinatensystem anpassen Auflösung modifizieren
sr <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"     # Projektion von AU3_2016_s
AU3_2019_proj <- projectRaster(AU3_2019_m, crs = sr)

S_2016_c1 <- projectRaster(S_2016_c, crs = sr)
S_2019a_c1 <- projectRaster(S_2019a_c, crs = sr)
S_2019b_c1 <- projectRaster(S_2019b_c, crs = sr)

# ===== 2.2 Auflösung modifizieren
  AU3_2016_proj_ngb_10 <- projectRaster(AU3_2016_m, crs = sr, res = 10, method = "ngb")
  AU3_2016_proj_bil_10 <- projectRaster(AU3_2016_m, crs = sr, res = 10, method = "bilinear")
  
  AU3_2019_proj_ngb_10 <- projectRaster(AU3_2019_m, crs = sr, res = 10, method = "ngb")
  AU3_2019_proj_bil_10 <- projectRaster(AU3_2019_m, crs = sr, res = 10, method = "bilinear")
  
# ===== 2.3 Normalisierung [Wertebereich 0-1] 
normalize_f <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  nrm <- ((x - min) / (max - min))
  return(nrm)
}

AU3_2016_s_norm <- normalize_f(AU3_2016_m)

  AU3_2016_proj_ngb_10_norm <- normalize_f(AU3_2016_proj_ngb_10)
  AU3_2016_proj_bil_10_norm <- normalize_f(AU3_2016_proj_bil_10)

AU3_2019_proj_norm <- normalize_f(AU3_2019_proj)
 
  AU3_2019_proj_ngb_10_norm <- normalize_f(AU3_2019_proj_ngb_10)
  AU3_2019_proj_bil_10_norm <- normalize_f(AU3_2019_proj_bil_10)

S_2016_c_n <- normalize_f(S_2016_c1)
S_2019a_c_n <- normalize_f(S_2019a_c1)
S_2019b_c_n <- normalize_f(S_2019b_c1)

# === 3 Vegetationsindizes =====================================================================================================
# Definition der angewandten Vegetationsindizes

dvi_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- (bn - br)
  return(vi)
}

evi_f <- function(img, n, r, b) {
  bn <- img[[n]]
  br <- img[[r]]
  bb <- img[[b]]
  vi <- 2.5*((bn - br) / (bn + 6 * br - 7.5 * bb +1))
  return(vi)
}

gdvi_f <- function(img, n, g) {
  bn <- img[[n]]
  bg <- img[[g]]
  vi <- (bn - bg)
  return(vi)
}

gndvi_f <- function(img, n, g) {
  bn <- img[[n]]
  bg <- img[[g]]
  vi <- (bn - bg) / (bn + bg)
  return(vi)
}

lai_f <- function(img, n, r, b) {
  bn <- img[[n]]
  br <- img[[r]]
  bb <- img[[b]]
  vi <- (3.618*(2.5*((bn - br) / (bn + 6 * br - 7.5 * bb +1)))-0.118)
  return(vi)
}

msavi2_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- (2 * bn + 1 - sqrt (( 2 * bn + 1)^2 - 8 * (bn-br)))/2
  return(vi)
}

ndvi_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- (bn - br) / (bn + br)
  return(vi)
}

osavi_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- ((bn - br) / (bn + br + 0.16))   
  return(vi)
}

rdvi_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- (bn - br) / sqrt(bn + br)
  return(vi)
}

wdrvi_f <- function(img, n, r) {
  bn <- img[[n]]
  br <- img[[r]]
  vi <- ((0.1*bn - br) / (0.1*bn + br))
  return(vi)
}

# === 4 Auswertung VI ==========================================================================================================
# Berechnung und statistische Betrachtung der Verteilungen je Aufnahmejahr und Vegetationsindizes
# ===== 4.1) UAV ==============================================================================================================
# ===== NDVI =====
AU3_2016_ndvi <- ndvi_f(AU3_2016_s_norm, 4, 2)
AU3_2016_proj_ngb_10_ndvi <- ndvi_f(AU3_2016_proj_ngb_10_norm, 4, 2)
AU3_2016_proj_bil_10_ndvi <- ndvi_f(AU3_2016_proj_bil_10_norm, 4, 2)

AU3_2019_ndvi <- ndvi_f(AU3_2019_proj_norm, 5, 3)
AU3_2019_proj_ngb_10_ndvi <- ndvi_f(AU3_2019_proj_ngb_10_norm, 4, 2)
AU3_2019_proj_bil_10_ndvi <- ndvi_f(AU3_2019_proj_bil_10_norm, 4, 2)

#--Save Raster       
r <- writeRaster(x = AU3_2016_ndvi,
                 filename = "A1/UAV_ndvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_ndvi,
                 filename = "A1/UAV_ndvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------GNDVI--------------------------------------
AU3_2016_gndvi <- gndvi_f(AU3_2016_s_norm, 4, 1)
AU3_2019_gndvi <- gndvi_f(AU3_2019_proj_norm, 5, 2)

#--Save Raster       
r <- writeRaster(x = AU3_2016_gndvi,
                 filename = "A1/UAV_gndvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_gndvi,
                 filename = "A1/UAV_gndvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------MSAVI2----
AU3_2016_msavi2 <- msavi2_f(AU3_2016_s_norm, 4, 2)
AU3_2019_msavi2 <- msavi2_f(AU3_2019_proj_norm, 5, 3)

#--Save Raster       
r <- writeRaster(x = AU3_2016_msavi2,
                 filename = "A1/UAV_msavi2_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_msavi2,
                 filename = "A1/UAV_msavi2_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------DVI----------------------------------------
AU3_2016_dvi <- dvi_f(AU3_2016_s_norm, 4, 2)
AU3_2019_dvi <- dvi_f(AU3_2019_proj_norm, 5, 3)

#--Save Raster       
r <- writeRaster(x = AU3_2016_dvi,
                 filename = "A1/UAV_dvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_dvi,
                 filename = "A1/UAV_dvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------GDVI---------------------------------------
AU3_2016_gdvi <- gdvi_f(AU3_2016_s_norm, 4, 1)
AU3_2019_gdvi <- gdvi_f(AU3_2019_proj_norm, 5, 2)

#--Save Raster       
r <- writeRaster(x = AU3_2016_gdvi,
                 filename = "A1/UAV_gdvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_gdvi,
                 filename = "A1/UAV_gdvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------OSAVI----
AU3_2016_osavi <- osavi_f(AU3_2016_s_norm, 4, 2)
AU3_2019_osavi <- osavi_f(AU3_2019_proj_norm, 5, 3)

#--Save Raster       
r <- writeRaster(x = AU3_2016_osavi,
                 filename = "A1/UAV_osavi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_osavi,
                 filename = "A1/UAV_osavi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------RDVI (for high healthy vegetation)----------
AU3_2016_rdvi <- rdvi_f(AU3_2016_s_norm, 4, 2)
AU3_2019_rdvi <- rdvi_f(AU3_2019_proj_norm, 5, 3)

#--Save Raster       
r <- writeRaster(x = AU3_2016_rdvi,
                 filename = "A1/UAV_rdvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_rdvi,
                 filename = "A1/UAV_rdvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------WDRVI (for wide range of vegetation)--------
AU3_2016_wdrvi <- wdrvi_f(AU3_2016_s_norm, 4, 2)
AU3_2019_wdrvi <- wdrvi_f(AU3_2019_proj_norm, 5, 3)

#--Save Raster       
r <- writeRaster(x = AU3_2016_wdrvi,
                 filename = "A1/UAV_wdrvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = AU3_2019_wdrvi,
                 filename = "A1/UAV_wdrvi_2019.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

# ===== 4.2) Sentinel ==========================================================================================================
#------------NDVI-----------------------
S_2016_ndvi <- ndvi_f(S_2016_c_n, 4, 3)
S_2019a_ndvi <- ndvi_f(S_2019a_c_n, 4, 3)
S_2019b_ndvi <- ndvi_f(S_2019b_c_n, 4, 3)

#--Save Raster       
r <- writeRaster(x = S_2016_ndvi,
                 filename = "A1/Sentinel_ndvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019a_ndvi,
                 filename = "A1/Sentinel_ndvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019b_ndvi,
                 filename = "A1/Sentinel_ndvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------GNDVI--------------------------------------
S_2016_gndvi <- gndvi_f(S_2016_c_n, 4, 2)
S_2019a_gndvi <- gndvi_f(S_2019a_c_n, 4, 2)
S_2019b_gndvi <- gndvi_f(S_2019b_c_n, 4, 2)

#--Save Raster       
r <- writeRaster(x = S_2016_gndvi,
                 filename = "A1/Sentinel_gndvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019a_gndvi,
                 filename = "A1/Sentinel_gndvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019b_gndvi,
                 filename = "A1/Sentinel_gndvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------msavi2----
S_2016_msavi2 <- msavi2_f(S_2016_c_n, 4, 3)
S_2019a_msavi2 <- msavi2_f(S_2019a_c_n, 4, 3)
S_2019b_msavi2 <- msavi2_f(S_2019b_c_n, 4, 3)

#--Save Raster       
r <- writeRaster(x = S_2016_msavi2,
                 filename = "A1/S_msavi2_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019a_msavi2,
                 filename = "A1/S_msavi2_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019b_msavi2,
                 filename = "A1/S_msavi2_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------DVI----------------------------------------
S_2016_dvi <- dvi_f(S_2016_c_n, 4, 3)
S_2019a_dvi <- dvi_f(S_2019a_c_n, 4, 3)
S_2019b_dvi <- dvi_f(S_2019b_c_n, 4, 3)

#--Save Raster       
r <- writeRaster(x = S_2016_dvi,
                 filename = "A1/Sentinel_dvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019a_dvi,
                 filename = "A1/Sentinel_dvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019b_dvi,
                 filename = "A1/Sentinel_dvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------GDVI---------------------------------------
S_2016_gdvi <- gdvi_f(S_2016_c_n, 4, 2)
S_2019a_gdvi <- gdvi_f(S_2019a_c_n, 4, 2)
S_2019b_gdvi <- gdvi_f(S_2019b_c_n, 4, 2)

#--Save Raster       
r <- writeRaster(x = S_2016_gdvi,
                 filename = "A1/Sentinel_gdvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019a_gdvi,
                 filename = "A1/Sentinel_gdvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

r <- writeRaster(x = S_2019b_gdvi,
                 filename = "A1/Sentinel_gdvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------OSAVI----
S_2016_osavi <- osavi_f(S_2016_c_n, 4, 3)
S_2019a_osavi <- osavi_f(S_2019a_c_n, 4, 3)
S_2019b_osavi <- osavi_f(S_2019b_c_n, 4, 3)


#--Save Raster       
r <- writeRaster(x = S_2016_osavi,
                 filename = "A1/S_osavi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019a_osavi,
                 filename = "A1/S_osavi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019b_osavi,
                 filename = "A1/S_osavi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------RDVI (for high healthy vegetation)----------
S_2016_rdvi <- rdvi_f(S_2016_c_n, 4, 3)
S_2019a_rdvi <- rdvi_f(S_2019a_c_n, 4, 3)
S_2019b_rdvi <- rdvi_f(S_2019b_c_n, 4, 3)

#--Save Raster       
r <- writeRaster(x = S_2016_rdvi,
                 filename = "A1/Sentinel_rdvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019a_rdvi,
                 filename = "A1/Sentinel_rdvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019b_rdvi,
                 filename = "A1/Sentinel_rdvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

#------------WDRVI (for wide range of vegetation)--------
S_2016_wdrvi <- wdrvi_f(S_2016_c_n, 4, 3)
S_2019a_wdrvi <- wdrvi_f(S_2019a_c_n, 4, 3)
S_2019b_wdrvi <- wdrvi_f(S_2019b_c_n, 4, 3)

#--Save Raster       
r <- writeRaster(x = S_2016_wdrvi,
                 filename = "A1/Sentinel_wdrvi_2016.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019a_wdrvi,
                 filename = "A1/Sentinel_wdrvi_2019a.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)
r <- writeRaster(x = S_2019b_wdrvi,
                 filename = "A1/Sentinel_wdrvi_2019b.tif",
                 format = "GTiff",
                 datatype = 'FLT4S',
                 overwrite = TRUE)

# === 5.1) UAV Change-Detection 2016-2019==============================================================================================
# ===== 5.1) UAV
# Auflösung anpassen: AU3_2016_... auf Pixelgröße von AU3_2019_... bringen

AU3_2016_ndvi_rs <- projectRaster(AU3_2016_ndvi, AU3_2019_ndvi, method = 'ngb')
AU3_2016_gdvi_rs <- projectRaster(AU3_2016_gdvi, AU3_2019_gdvi, method = 'ngb')
AU3_2016_msavi2_rs <- projectRaster(AU3_2016_msavi2, AU3_2019_msavi2, method = 'ngb')

AU3_2016_dvi_rs <- projectRaster(AU3_2016_dvi, AU3_2019_dvi, method = 'ngb') 
AU3_2016_gndvi_rs <- projectRaster(AU3_2016_gndvi, AU3_2019_gndvi, method = 'ngb')
AU3_2016_osavi_rs <- projectRaster(AU3_2016_osavi, AU3_2019_osavi, method = 'ngb')
AU3_2016_rdvi_rs <- projectRaster(AU3_2016_rdvi, AU3_2019_rdvi, method = 'ngb')
AU3_2016_wdrvi_rs <- projectRaster(AU3_2016_wdrvi, AU3_2019_wdrvi, method = 'ngb')

# Veränderung berechnen (2019-2016)
AU3_ndvi_change <- AU3_2019_ndvi - AU3_2016_ndvi_rs
AU3_gdvi_change <- AU3_2019_gdvi - AU3_2016_gdvi_rs
AU3_msavi2_change <- AU3_2019_msavi2 - AU3_2016_msavi2_rs

AU3_dvi_change <- AU3_2019_dvi - AU3_2016_dvi_rs
AU3_gndvi_change <- AU3_2019_gndvi - AU3_2016_gndvi_rs
AU3_gi2_change <- AU3_2019_gi2 - AU3_2016_gi2_rs
AU3_osavi_change <- AU3_2019_osavi - AU3_2016_osavi_rs
AU3_rdvi_change <- AU3_2019_rdvi - AU3_2016_rdvi_rs
AU3_wdrvi_change <- AU3_2019_wdrvi - AU3_2016_wdrvi_rs

# ===== Ergebnisdarstellung Text =====
# Definition der Breaks (Wertebereich Symbologie)
breakpoints <- c( -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.0,  0.1,  0.2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) 
breakpoints1 <- c( 0.0, 0.05,  0.1, 0.15,  0.2, 0.25,  0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1) 
breakpoints2 <- c( -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2) 
breakpoints3 <- c(-Inf, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf)
breakpoints4 <- c(-30, -27, -24, -21, -18, -15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30)
color <- rev(terrain.colors(21))
color1 <- rev(terrain.colors(23))

# ===== NDVI =====
par(mfrow = c(1,2))
plot(AU3_2016_ndvi, main = "UAV NDVI 2016", breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_ndvi, main = "UAV NDVI 2019", breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
dev.off()

#--Statistic
summary(AU3_2016_ndvi)
summary(AU3_2019_ndvi)
summary(AU3_ndvi_change)

par(mfrow = c(1,2))
hist(AU3_2016_ndvi,
     main = "UAV 2016",
     xlab = "NDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-0.5, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.25), labels = seq(-0.5,1, 0.25))
abline(v = 0.6529, col = "red")
legend("topleft", legend = c("Median"), text.font = 2, lty = 1, lwd = 1.5, col = "red", cex = 0.5)

hist(AU3_2019_ndvi,
     main = "UAV 2019",
     xlab = "NDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-0.5, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.25), labels = seq(-0.5,1, 0.25))
abline(v = 0.2792, col = "red")
legend("topright", legend = c("Median"), text.font = 2, lty = 1, lwd = 1.5, col = "red", cex = 0.5)
dev.off()

par(mfrow = c(1,2))
plot(AU3_ndvi_change)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
hist(AU3_ndvi_change,
     main = "",
     xlab = "NDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-2, 2),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-2,2, 0.25), labels = seq(-2,2, 0.25))
abline(v = -0.3718, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topright", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)
dev.off()

# ===== GDVI =====
par(mfrow = c(1,2))
plot(AU3_2016_gdvi, main = "UAV GDVI 2016", breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_gdvi, main = "UAV GDVI 2019", breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

#--Statistic
summary(AU3_2016_gdvi)
summary(AU3_2019_gdvi)
summary(AU3_gdvi_change)

hist(AU3_2016_gdvi,
     main = "2016",
     xlab = "GDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 40,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.3403, col = "red")
legend("topleft", legend = c("Median"),lty = 1, col = "red", cex = 0.5)

hist(AU3_2019_gdvi,
     main = "2019",
     xlab = "GDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.1050, col = "red")
legend("topright", legend = c("Median"),lty = 1, col = "red", cex = 0.5)

plot(AU3_gdvi_change)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
hist(AU3_gdvi_change,
     main = "",
     xlab = "GDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 40,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.2246, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topright", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)
dev.off()

# ===== OSAVI =====
par(mfrow = c(1,2))
plot(AU3_2016_osavi, breaks = breakpoints, col = color, main = "UAV OSAVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_osavi, breaks = breakpoints, col = color, main = "UAV OSAVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

#--Statistic
summary(AU3_2016_osavi)
summary(AU3_2019_osavi)
summary(AU3_osavi_change)

par(mfrow = c(1,2))
hist(AU3_2016_osavi,
     main = "UAV OSAVI 2016",
     xlab = "OSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.5199, col = "red")
legend("topleft", legend = c("Median"),lty = 1, col = "red", cex = 0.5)

hist(AU3_2019_osavi,
     main = "UAV OSAVI 2019",
     xlab = "OSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.2028, col = "red")
legend("topleft", legend = c("Median"),lty = 1, col = "red", cex = 0.5)
dev.off()

par(mfrow = c(1,2))
plot(AU3_osavi_change)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
hist(AU3_osavi_change,
     main = "",
     xlab = "OSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1.5, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1.5,1, 0.25), labels = seq(-1.5,1, 0.25))
abline(v = -0.3138, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)
dev.off()

# ===== Ergebnisdarstellung Anhang =====
# ===== DVI =====
par(mfrow = c(2,3))
plot(AU3_2016_dvi, breaks = breakpoints, col = color, main = "UAV DVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_dvi, breaks = breakpoints, col = color, main = "UAV DVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_dvi_change, col = color, breaks = breakpoints2, main = "UAV DVI \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(AU3_2016_dvi)
summary(AU3_2019_dvi)
summary(AU3_dvi_change)

hist(AU3_2016_dvi,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.25), labels = seq(-0.5,1, 0.25))
abline(v = 0.4698, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(AU3_2019_dvi,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.25), labels = seq(-0.5,1, 0.25))
abline(v = 0.1238, col = "red", lty = 2)

hist(AU3_dvi_change,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.3275, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== GNDVI =====
par(mfrow = c(2,3))
plot(AU3_2016_gndvi,breaks = breakpoints, col = color, main = "UAV GNDVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_gndvi, breaks = breakpoints, col = color, main = "UAV GNDVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_gndvi_change, col = color, main = "UAV GNDVI \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(AU3_2016_gndvi)
summary(AU3_2019_gndvi)
summary(AU3_gndvi_change)

hist(AU3_2016_gndvi,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 20,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.4393, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(AU3_2019_gndvi,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 20,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.2311, col = "red", lty = 2)

hist(AU3_gndvi_change,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-2, 2),
     breaks = 20,
     xaxt = 'n')
axis(side = 1, at = seq(-2,2, 0.25), labels = seq(-2,2, 0.25))
abline(v = -0.2050, col = "red", lty = 2)
abline(v = 0, col = "red")

dev.off()

# ===== MSAVI2 =====
par(mfrow = c(2,3))
plot(AU3_2016_msavi2, breaks = breakpoints, col = color, main = "UAV MSAVI2 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_msavi2, breaks = breakpoints, col = color, main = "UAV MSAVI2 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_msavi2_change, col = color, main = "UAV MSAVI2 \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(AU3_2016_msavi2)
summary(AU3_2019_msavi2)
summary(AU3_msavi2_change)

hist(AU3_2016_msavi2,
     main = "",
     xlab = "MSAVI2",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.5603, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(AU3_2019_msavi2,
     main = "",
     xlab = "MSAVI2",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-0.5, 1),
     breaks = 40,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.1), labels = seq(-0.5,1, 0.1))
abline(v = 0.1743, col = "red", lty = 2)

hist(AU3_msavi2_change,
     main = "",
     xlab = "MSAVI2",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.3709, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== RDVI ====
par(mfrow = c(2,3))
plot(AU3_2016_rdvi, breaks = breakpoints, col = color, main = "RDVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_rdvi, breaks = breakpoints, col = color, main = "RDVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_rdvi_change, col = color, main = "RDVI Bilanz \n (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(AU3_2016_rdvi)
summary(AU3_2019_rdvi)
summary(AU3_rdvi_change)

hist(AU3_2016_rdvi,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.5423, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(AU3_2019_rdvi,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.1909, col = "red", lty = 2)

hist(AU3_rdvi_change,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.3470, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== WDRVI ====
par(mfrow = c(2,3))
plot(AU3_2016_wdrvi, breaks = breakpoints, col = color, main = "WDRVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_wdrvi, breaks = breakpoints, col = color, main = "WDRVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_wdrvi_change, col = color, main = "WDRVI Bilanz \n (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(AU3_2016_wdrvi)
summary(AU3_2019_wdrvi)
summary(AU3_wdrvi_change)

hist(AU3_2016_wdrvi,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.3548, col = "red", lty = 2)
legend("topright", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(AU3_2019_wdrvi,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = -0.6986, col = "red", lty = 2)

hist(AU3_wdrvi_change,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-2, 2),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-2,2, 0.25), labels = seq(-2,2, 0.25))
abline(v = -0.3398, col = "red", lty = 2)
abline(v = 0, col = "red")

dev.off()

# ===== 5.2) Sentinel change detection =====
S_ndvi_change_a <- S_2019a_ndvi - S_2016_ndvi
S_gdvi_change_a <- S_2019a_gdvi - S_2016_gdvi
S_osavi_change_a <- S_2019a_osavi - S_2016_osavi

S_ndvi_change_b <- S_2019b_ndvi - S_2016_ndvi
S_gdvi_change_b <- S_2019b_gdvi - S_2016_gdvi
S_osavi_change_b <- S_2019b_osavi - S_2016_osavi

S_dvi_change <- S_2019a_dvi - S_2016_dvi
S_gndvi_change <- S_2019a_gndvi - S_2016_gndvi
S_lai_change <- S_2019a_lai - S_2016_lai
S_msavi2_change <- S_2019a_msavi2 - S_2016_msavi2
S_rdvi_change <- S_2019a_rdvi - S_2016_rdvi
S_wdrvi_change <- S_2019a_wdrvi - S_2016_wdrvi

# ===== Sentinel Ergebnisdarstellung Text =====
# ===== S NDVI =====
par(mfrow = c(1,3))#, mai = c(0, 0.5, 0.5, 0.5))
plot(S_2016_ndvi, breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_ndvi, breaks = breakpoints, col = color, main = "NDVI Sentinel-2")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019b_ndvi, breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
mtext("NDVI (Sentinel-2)", outer = TRUE, cex = 1.2)
dev.off()

#--Statistic
summary(S_2016_ndvi)
summary(S_2019a_ndvi)
summary(S_2019b_ndvi)
summary(S_ndvi_change_a)
summary(S_ndvi_change_b)

par(mfrow = c(3,1), mai = c(0.6, 0.8, 0.5, 0.9))
hist(S_2016_ndvi,
     main = "",
     xlab = "NDVI Sentinel-2 26/08/2016",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0222, col = "red")
legend("topright", legend = c("Median"), lty = 1, lwd = 1.5, col = "red", cex = 0.5)

hist(S_2019a_ndvi,
     main = "",
     xlab = "NDVI Sentinel-2 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.5038, col = "red")

hist(S_2019b_ndvi,
     main = "",
     xlab = "NDVI Sentinel-2 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.2407, col = "red")
dev.off()

par(mfrow = c(2,2), mai = c(0.5, 0.5, 0.5, 0.5))
plot(S_ndvi_change_a)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

plot(S_ndvi_change_b)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

hist(S_ndvi_change_a,
     main = "",
     xlab = "NDVI 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-2, 2),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-2, 2, 0.25), labels = seq(-2, 2, 0.25))
abline(v = 0.4760, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_ndvi_change_b,
     main = "",
     xlab = "NDVI 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-2, 2),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-2, 2, 0.25), labels = seq(-2, 2, 0.25))
abline(v = 0.1949, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S GDVI =====
par(mfrow = c(1,3))
plot(S_2016_gdvi,  breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_gdvi, main = "GDVI Sentinel-2", breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019b_gdvi, breaks = breakpoints, col = color)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
dev.off()

#--Statistic
summary(S_2016_gdvi)
summary(S_2019a_gdvi)
summary(S_2019b_gdvi)
summary(S_gdvi_change_a)
summary(S_gdvi_change_b)

par(mfrow = c(3,1), mai = c(0.6, 0.5, 0.2, 0.5))
hist(S_2016_gdvi,
     main = "",
     xlab = "GDVI Sentinel-2 26/08/2016",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.2796, col = "red")
legend("topright", legend = c("Median"), lty = 1, lwd = 1.5, col = "red", cex = 0.8)

hist(S_2019a_gdvi,
     main = "",
     xlab = "GDVI Sentinel-2 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.2818, col = "red")

hist(S_2019b_gdvi,
     main = "",
     xlab = "GDVI Sentinel-2 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.1216, col = "red")
dev.off()

par(mfrow = c(2,2), mai = c(0.7, 0.5, 0.5, 0.5))
plot(S_gdvi_change_a)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

plot(S_gdvi_change_b)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

hist(S_gdvi_change_a,
     main = "",
     xlab = "GDVI 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0313, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_gdvi_change_b,
     main = "",
     xlab = "GDVI 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = -0.1356, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S OSAVI =====
par(mfrow = c(1,3))
plot(S_2016_osavi, breaks = breakpoints, col = color, main = "Sentinel-2 OSAVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_osavi, breaks = breakpoints, col = color, main = "Sentinel-2 OSAVI 08/2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019b_osavi, breaks = breakpoints, col = color, main = "Sentinel-2 OSAVI 09/2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

#--Statistic
summary(S_2016_osavi)
summary(S_2019a_osavi)
summary(S_2019b_osavi)
summary(S_osavi_change_a)
summary(S_osavi_change_b)

par(mfrow = c(3,1), mai = c(0.6, 0.5, 0.2, 0.5))
hist(S_2016_osavi,
     main = "",
     xlab = "OSAVI Sentinel-2 26/08/2016",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0193, col = "red")
legend("topright", legend = c("Median"), lty = 1, lwd = 1.5, col = "red", cex = 0.8)

hist(S_2019a_osavi,
     main = "",
     xlab = "OSAVI Sentinel-2 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.4190, col = "red")

hist(S_2019b_osavi,
     main = "",
     xlab = "OSAVI Sentinel-2 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 50,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.2072, col = "red")
dev.off()

par(mfrow = c(2,2), mai = c(0.7, 0.5, 0.5, 0.5))
plot(S_osavi_change_a)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

plot(S_osavi_change_b)
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

hist(S_osavi_change_a,
     main = "",
     xlab = "OSAVI 26/08/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.3963, col = "red", lty = 2)
abline(v = 0, col = "red")
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_osavi_change_b,
     main = "",
     xlab = "OSAVI 12/09/2019",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 15,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.1645, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== Sentinel Ergebnisdarstellung Anhang =====
# ===== S DVI =====
par(mfrow = c(2,3))
plot(S_2016_dvi, col = color, breaks = breakpoints, main = "Sentinel-2 DVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_dvi, breaks = breakpoints, col = color, main = "Sentinel-2 DVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_dvi_change, col = color, breaks = breakpoints, main = "Sentinel-2 DVI \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_dvi)
summary(S_2019a_dvi)
summary(S_dvi_change)

hist(S_2016_dvi,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 20,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0200, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_dvi,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = 0.4292, col = "red", lty = 2)

hist(S_dvi_change,
     main = "",
     xlab = "DVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.4068, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S GNDVI =====
par(mfrow = c(2,3))
plot(S_2016_gndvi,breaks = breakpoints, col = color, main = "Sentinel-2 GNDVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_gndvi, breaks = breakpoints, col = color, main = "Sentinel-2 GNDVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_gndvi_change, breaks = breakpoints, col = color, main = "Sentinel-2 GNDVI \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_gndvi)
summary(S_2019a_gndvi)
summary(S_gndvi_change)

hist(S_2016_gndvi,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.3033, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_gndvi,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = 0.2914, col = "red", lty = 2)

hist(S_gndvi_change,
     main = "",
     xlab = "GNDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.0139, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S LAI =====
par(mfrow = c(2,3))
plot(S_2016_lai, breaks = breakpoints4, col = color, main = "Sentinel-2 LAI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_lai, breaks = breakpoints4, col = color, main = "Sentinel-2 LAI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_lai_change, col = color, main = "Sentinel-2 LAI \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_lai)
summary(S_2019a_lai)
summary(S_lai_change)

hist(S_2016_lai,
     main = "",
     xlab = "LAI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.3033, col = "red")
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_lai,
     main = "",
     xlab = "LAI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = 0.2914, col = "red")

hist(S_lai_change,
     main = "",
     xlab = "LAI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.0139, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S MSAVI2 =====
par(mfrow = c(2,3))
plot(S_2016_msavi2, breaks = breakpoints, col = color, main = "Sentinel-2 MSAVI2 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_msavi2, breaks = breakpoints, col = color, main = "Sentinel-2 MSAVI2 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_msavi2_change, breaks = breakpoints, col = color, main = "Sentinel-2 MSAVI2 \n Bilanz (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_msavi2)
summary(S_2019a_msavi2)
summary(S_msavi2_change)

hist(S_2016_msavi2,
     main = "",
     xlab = "MSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0202, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_msavi2,
     main = "",
     xlab = "MSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = 0.4761, col = "red", lty = 2)

hist(S_msavi2_change,
     main = "",
     xlab = "MSAVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.4311, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()
# ===== S RDVI ====
par(mfrow = c(2,3))
plot(S_2016_rdvi, breaks = breakpoints, col = color, main = "RDVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_rdvi,breaks = breakpoints, col = color, main = "RDVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_rdvi_change, breaks = breakpoints, col = color, main = "RDVI Bilanz \n (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_rdvi)
summary(S_2019a_rdvi)
summary(S_rdvi_change)


hist(S_2016_rdvi,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = 0.0221, col = "red", lty = 2)
legend("topleft", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_rdvi,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = 0.4719, col = "red", lty = 2)

hist(S_rdvi_change,
     main = "",
     xlab = "RDVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.4279, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== S WDRVI ====
par(mfrow = c(2,3))
plot(S_2016_wdrvi, breaks = breakpoints, col = color, main = "WDRVI 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_wdrvi, breaks = breakpoints, col = color, main = "WDRVI 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_wdrvi_change, breaks = breakpoints, col = color, main = "WDRVI Bilanz \n (2019-2016)")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

summary(S_2016_wdrvi)
summary(S_2019a_wdrvi)
summary(S_wdrvi_change)


hist(S_2016_wdrvi,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq(-1, 1, 0.25))
abline(v = -0.7289, col = "red", lty = 2)
legend("topright", legend = c("Median", "Wert 0"),lty = c(2, 1), col = c("red", "red"), cex = 0.5)

hist(S_2019a_wdrvi,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1, 1, 0.25), labels = seq( -1, 1, 0.25))
abline(v = -0.3750, col = "red", lty = 2)

hist(S_wdrvi_change,
     main = "",
     xlab = "WDRVI",
     ylab= "Häufigkeit",
     col = "darkgreen",
     xlim = c(-1, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-1,1, 0.25), labels = seq(-1,1, 0.25))
abline(v = 0.3436, col = "red", lty = 2)
abline(v = 0, col = "red")
dev.off()

# ===== 6 Überwachte Klassifikation===============================================================================================
# ===== 6.1) 2016
sink("Endergebnisse/Klassifikation_2016_rf_1.txt")

UAV_ndvi_2016_rf <- superClass(AU3_2016_ndvi, trainData = train2016[13:26], responseCol = 2, 
                               nSamples = 1000, model = "rf", kfold = 5, tuneLength = 5, 
                               mode = "classification", verbose = T,  filename = NULL
                                )

UAV_ndvi_2016_knn <- superClass(AU3_2016_ndvi, trainData = train2016[13:26], responseCol = 2, 
                               nSamples = 1000, model = "knn", kfold = 5, tuneLength = 5, 
                               mode = "classification", verbose = T,  filename = NULL
                                )

dev.off()

# ===== 6.2) 2019
sink("Endergebnisse/Klassifikation_2019_rf.txt")
UAV_ndvi_2019_rf <- superClass(AU3_2019_ndvi, trainData = train2019[14:28], responseCol = 1, 
                               nSamples = 1000, model = "rf", kfold = 5, tuneLength = 3, 
                               mode = "classification", verbose = T,  filename = NULL
)

UAV_ndvi_2019_knn <- superClass(AU3_2019_ndvi, trainData = train2016[14:28], responseCol = 1, 
                                nSamples = 1000, model = "knn", kfold = 5, tuneLength = 3, 
                                mode = "classification", verbose = T,  filename = NULL
)

dev.off()

#==================Sentinel
trainData2016 <- train2016[2,14:28]

sink("Endergebnisse/Klassifikation_2016_rf_Sentinel.txt")

S_ndvi_2016_rf <- superClass(S_2016_ndvi, trainData = train2016[13:26], responseCol = 2, 
                                nSamples = 1000, model = "rf", kfold = 5, tuneLength = 3, 
                                mode = "classification", verbose = T,  filename = NULL
)

S_ndvi_2016_knn <- superClass(S_2016_ndvi, trainData = train2016[13:26], responseCol = 2, 
                                nSamples = 1000, model = "knn", kfold = 5, tuneLength = 3, 
                                mode = "classification", verbose = T,  filename = NULL
)

S_ndvi_2019_rf <- superClass(S_2019a_ndvi, trainData = train2016[14:28], responseCol = 1, 
                               nSamples = 1000, model = "rf", kfold = 5, tuneLength = 3, 
                               mode = "classification", verbose = T,  filename = NULL
)

S_ndvi_2019_knn <- superClass(S_2019a_ndvi, trainData = train2016[14:28], responseCol = 1, 
                                nSamples = 1000, model = "knn", kfold = 5, tuneLength = 3, 
                                mode = "classification", verbose = T,  filename = NULL
)

dev.off()

# ===== 7 Vergleich UAV Sentinel =====

par(mfrow = c(2,2))
plot(AU3_2016_proj_ngb_10_ndvi, breaks = breakpoints, col = color, main = "UAV 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(AU3_2019_proj_ngb_10_ndvi,  breaks = breakpoints, col = color, main = "UAV 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2016_ndvi, breaks = breakpoints, col = color,  main = "Sentinel 2016")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)
plot(S_2019a_ndvi,  breaks = breakpoints, col = color, main = "Sentinel 2019")
par(new = T)
plot(UG, axes = T,
     border = "black",
     add = T)

par(mfrow = c(2,2))
hist(AU3_2016_proj_ngb_10_ndvi, main = "UAV 2016", breaks = 20, col = "darkgreen")
hist(AU3_2019_proj_ngb_10_ndvi, main = "UAV 2019", breaks = 20, col = "darkgreen")
hist(S_2016_ndvi, main = "Sentinel 2016", breaks = 20, col = "darkgreen")
hist(S_2019a_ndvi, main = "Sentinel 2019", breaks = 20, col = "darkgreen")

sink("Endergebnisse/Klassifikation_2016_projectRaster_rf_.txt")

dev.off()