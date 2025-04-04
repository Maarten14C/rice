# R has to be within the working directory of this file

# extracted from calib marine database on 2 April 2025:
# http://calib.org/marine/query select * from details
# 1979 records
marinedata <- read.csv("marinedatabase.csv") # 247 kB

# references: downloaded results of 'select * from refs' and saved as .csv file
refs <- read.csv("refs.csv")
# taxa: downloaded results of 'select * from taxa' and saved as .csv file
taxa <- read.csv("taxa.csv")

for(i in 1:nrow(refs)) {
  these <- which(marinedata$RefNo == refs[i,1])	
  marinedata[these,21] <- paste0(refs[i,2], ", ", refs[i,5],  ". ", refs[i,3], ". ", refs[i,4], " ", refs[i,6])	
}

# cleaning up
taxa[which(taxa[,4] == "Suspension"),4] <- "suspension"
taxa[which(taxa[,4] == "Deposit"),4] <- "deposit"
taxa[which(taxa[,4] == " deposit"),4] <- "deposit"
taxa[which(taxa[,4] == "Deposit feeder o"),4] <- "deposit"
taxa[which(taxa[,4] == "Carnivore"),4] <- "carnivore"
taxa[which(taxa[,4] == "Algal Grazer"),4] <- "algal grazer"
taxa[which(taxa[,4] == "NULL"),4] <- "unknown"
taxa[which(is.na(taxa[,4])),4] <- "unknown"
taxa[which(taxa[,4] == ""),4] <- "unknown"
taxa[which(taxa[,4] == "Unknown"),4] <- "unknown"

for(i in 1:nrow(taxa)) {
  these <- which(marinedata$RefNo == taxa[i,1])	
  marinedata[these,22] <- paste(taxa[i,2], taxa[i,3])
  marinedata[these,23] <- tolower(taxa[i,4])	
}
marinedata[which(is.na(marinedata[,23])),23] <- "unknown"

shells <- data.frame(
  lon = marinedata$Lon, 
  lat = marinedata$Lat,
  no = marinedata$MapNo,
  taxonN = marinedata$TaxaNo,
  dR = marinedata$DeltaR,
  dSTD = marinedata$DeltaRErr,
  collected = marinedata$CollectionYear,
  res = marinedata$ReservoirAge,
  res.error = marinedata$ReservoirErr,
  C14 = marinedata$C14age,
  er = marinedata$C14err,
  lab = marinedata$LabID,
  ref = marinedata[,21],
  taxon = marinedata[,22],
  feeding = marinedata[,23]
)

save(shells, file = "../data/shells.rda", compress = "bzip2")
