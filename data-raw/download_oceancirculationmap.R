for(i in 0:5)
  for(j in 0:13) {
    url <- paste0("https://wmts.marine.copernicus.eu/teroWmts/?service=WMTS&request=GetTile&version=2.0.0",
      "&layer=GLOBAL_ANALYSISFORECAST_PHY_001_024/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i&tilematrixset=EPSG:4326",
      i, "&tilecol=", j,
      "&STYLE=vectorStyle:solidAndVector,cmap:velocity")
    download.file(url, destfile = paste0("velocity_row",i, "_col", j, ".png"), mode = "wb")
  }

#https://wmts.marine.copernicus.eu/teroWmts/?service=WMTS&request=GetTile&version=2.0.0&layer=GLOBAL_ANALYSISFORECAST_PHY_001_024/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i_202406/sea_water_velocity&tilematrixset=EPSG:4326&tilematrix=2&tilerow=0&tilecol=0&time=2024-09-09T00:00:00.000000000Z&STYLE=vectorStyle:solidAndVector,cmap:velocity

