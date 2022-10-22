
####
####
####


#!/bin/env Rscript

### purpose:
### to convert kg kg-1 to g m-3

data_dir <- "../../s01.data/era5/out"
missout <- -999.99  ## for output
gas_constant <- 8.314  ## universal gas constant  [j mol-1 k-1]
dry_mass_air <- 28.9644/1000     ## mean molar mass of dry air [kg mol-1]

###
###
###

system("./clean.sh")
suppressPackageStartupMessages(library(ncdf4))

###
###
###

flist_ciwc <- dir(data_dir, pattern = "ciwc")
flist_temp <- dir(data_dir, pattern = "temp")

f = 1

cat("\n")

while (f <= length(flist_ciwc)) {

   cat("\n ... processing", flist_ciwc[f])

   fname_ciwc <- paste0(data_dir, "/", flist_ciwc[f])
   fname_temp <- paste0(data_dir, "/", flist_temp[f])

   nc_ciwc <- nc_open(fname_ciwc)
   nc_temp <- nc_open(fname_temp)

   ciwc <- ncvar_get(nc_ciwc, "ciwc")
   temp <- ncvar_get(nc_temp, "t") 
   
   title <- ncatt_get(nc_ciwc, 0, "title")
   institution <- ncatt_get(nc_ciwc, 0, "institution")
   history <- ncatt_get(nc_ciwc, 0, "history")
   conventions <- ncatt_get(nc_ciwc, 0, "Conventions")
   system <- ncatt_get(nc_ciwc, 0, "system")
   lat <- nc_ciwc$dim$latitude$vals
   lon <- nc_ciwc$dim$longitude$vals
   time <- nc_ciwc$dim$time$vals
   tunits <- nc_ciwc$dim$time$units
   mvalue <- ncatt_get(nc_ciwc, "ciwc", "missing_value")
   missing <- mvalue$value 

   plev <- nc_ciwc$dim$level$vals * 100  ## to pa
  
   temp[which(temp == missing)] <- NA
   ciwc[which(ciwc == missing)] <- NA

   iwc_gm3 <- array(0, dim(ciwc))

   for (i in 1: length(plev)) {
      iwc_gm3[,,i,] <- (ciwc[,,i, ] * plev[i] * 100 * dry_mass_air ) / (gas_constant * temp[,,i,]) * 1000
   }

   iwc_gm3[which(is.na(iwc_gm3) == T)] <- missout

   londim <- ncdim_def("longitude", "degrees_east",
                       as.double(lon),
                       longname = "longitude")
   latdim <- ncdim_def("latitude", "degrees_north",
                       as.double(lat),
                       longname = "latitude")
   pdim <- ncdim_def("level", "millibars",
                       as.double(plev),
                       longname = "pressure_level")
   timedim <- ncdim_def("time",
                        tunits,
                        as.double(time),
                        longname = "time",
                        unlim = T)
   fillvalue <- missout
   dlname <- "cloud liquid water content"

   iwc_gm3_def <- ncvar_def("ciwc", "g m-3",
                        list(londim, latdim, pdim, timedim),
                        fillvalue,
                        dlname,
                        prec = "float")

   ncout <- nc_create("outfile.nc",
                      list(iwc_gm3_def),
                      force_v4 = F)

   ncvar_put(ncout, iwc_gm3_def, iwc_gm3)
   ncatt_put(ncout, "longitude", "axis", "X")
   ncatt_put(ncout, "latitude", "axis", "Y")
   ncatt_put(ncout, "level", "axis", "Z")
   ncatt_put(ncout, "time", "axis", "T")
   ncatt_put(ncout,0,"title",title$value)
   ncatt_put(ncout,0,"system",system$value)
   ncatt_put(ncout,0,"history",history$value)
   ncatt_put(ncout,0,"Conventions",conventions$value)
   nc_close(ncout)

   file.rename("outfile.nc", flist_ciwc[f])
   
   f = f+1
}
cat("\n")
###
###
###
