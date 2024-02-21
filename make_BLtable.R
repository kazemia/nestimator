load("C:/Users/kazem/Documents/Data/myDataWithBL.rdata")
load("C:/Users/kazem/Documents/Data/ekstraBL.rdata")
data <- as.data.frame(my_data)
remove(my_data)


data <- data %>% filter(!is.na(bldt.x))
data$Z_value <- NA
data$m <- NA
NDPC_mid_date <- as.Date(c())
for(y in 2010:2013){
  StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d")
  EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d")
  NDPC_mid_date <- 
    c(NDPC_mid_date, 
      StartDate + floor(as.numeric(
        difftime(EndDate, StartDate, units = "days"))/2))
  data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
  data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
    as.numeric(difftime(data$bldt.x[
      (data$bldt.x > StartDate) & (data$bldt.x < EndDate)], 
      StartDate, units = "days"))
}
y <- 2014
StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d")
EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d")
NDPC_mid_date <- 
  c(NDPC_mid_date, 
    StartDate + floor(as.numeric(
      difftime(EndDate, StartDate, units = "days"))/2))
data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
  as.numeric(difftime(data$bldt.x[
    (data$bldt.x > StartDate) & (data$bldt.x < EndDate)],
    StartDate, units = "days"))
y <- 2015
StartDate <- as.Date(paste0(y, "-02-28"), "%Y-%m-%d")
EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d")
NDPC_mid_date <- 
  c(NDPC_mid_date, 
    StartDate + floor(as.numeric(
      difftime(EndDate, StartDate, units = "days"))/2))
data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
  as.numeric(difftime(data$bldt.x[
    (data$bldt.x > StartDate) & (data$bldt.x < EndDate)],
    StartDate, units = "days"))
y <- 2016
StartDate <- as.Date(paste0(y, "-02-29"), "%Y-%m-%d")
EndDate <- as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d")
NDPC_mid_date <- 
  c(NDPC_mid_date, 
    StartDate + floor(as.numeric(
      difftime(EndDate, StartDate, units = "days"))/2))
data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
  as.numeric(difftime(data$bldt.x[
    (data$bldt.x > StartDate) & (data$bldt.x < EndDate)],
    StartDate, units = "days"))
y <- 2017
StartDate <- as.Date(paste0(y, "-02-28"), "%Y-%m-%d")
EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d")
NDPC_mid_date <- 
  c(NDPC_mid_date, 
    StartDate + floor(as.numeric(
      difftime(EndDate, StartDate, units = "days"))/2))
data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
  as.numeric(difftime(data$bldt.x[
    (data$bldt.x > StartDate) & (data$bldt.x < EndDate)],
    StartDate, units = "days"))
for(y in 2018:2019){
  StartDate <- as.Date(paste0(y, "-01-31"), "%Y-%m-%d")
  EndDate <- as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d")
  NDPC_mid_date <- 
    c(NDPC_mid_date, 
      StartDate + floor(as.numeric(
        difftime(EndDate, StartDate, units = "days"))/2))
  data$Z_value[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- y-2009
  data$m[(data$bldt.x > StartDate) & (data$bldt.x < EndDate)] <- 
    as.numeric(difftime(data$bldt.x[
      (data$bldt.x > StartDate) & (data$bldt.x < EndDate)],
      StartDate, units = "days"))
}

data <- data %>% filter(!is.na(Z_value))
dataWithBL <- data
#Read the data used in the paper from make_adeff
dataWithBL <- dataWithBL %>% filter(tcid %in% data$tcid)
BL_table <- BL_table %>% filter(tcid %in% data$tcid)
colnames(BL_table) <- paste(colnames(BL_table), "BL", sep = "_")
BL_table$tcid <- BL_table$tcid_BL
dataWithBL <- merge(dataWithBL, BL_table, all.x = TRUE, by = "tcid")
dataWithBL <- unique(dataWithBL)
dataWithBL <- dataWithBL %>% filter(tcid %in% data$tcid)
summaryTable <- dataWithBL %>% group_by(trtgrp.x) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
