library(tidyverse)
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
#Make an adeff id
data$myid <- 1:nrow(data)
#merge the data sets
Missingreport <- dataWithBL %>% merge(data, all = TRUE, by = "tcid")
#identify rows in the new data set that aren't in adeff
Missingreport$myid[is.na(Missingreport$myid)] <- Missingreport$tcid[is.na(Missingreport$myid)]
#Remove the duplicates from the new data
Missingreport <- Missingreport[!duplicated(Missingreport$myid),]
#Remove the one patient that isn't in adeff without targetmissing = 2
Missingreport <- Missingreport %>% filter(!((tcid==myid) & (targetmissing != 2)))
#Patients that have target missing = 2 but were in adeff should have targetmissing = 1
Missingreport <- Missingreport %>% mutate(
  targetmissing = if_else((tcid!=myid) & (targetmissing == 2), 1, targetmissing))
#Patients that have target missing = 3 but reached remission in adeff should have targetmissing = 0
Missingreport <- Missingreport %>% mutate(
  targetmissing = if_else((das28crprem.y != 0) & (targetmissing == 3), 0, targetmissing))
Missingreport <- Missingreport %>% select(Z_value.y, trtgrp.x.y, das28crprem.y, targetmissing)
colnames(Missingreport) <- c("Z", "T", "Y", "M")
save(Missingreport, file = "C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator/data/Missingtable.Rdata")

#Make BL table
dataWithBL <- dataWithBL %>% filter(tcid %in% data$tcid)
BL_table <- BL_table %>% filter(tcid %in% data$tcid)
colnames(BL_table) <- paste(colnames(BL_table), "BL", sep = "_")
BL_table$tcid <- BL_table$tcid_BL
dataWithBL <- merge(dataWithBL, BL_table, all.x = TRUE, by = "tcid")
dataWithBL <- unique(dataWithBL)
dataWithBL <- dataWithBL %>% filter(tcid %in% data$tcid)
summaryTable <- dataWithBL %>% group_by(Z_value) %>% summarise_if(is.numeric, mean, na.rm = TRUE)
round(summaryTable$das28crp_BL, digits = 2)
round(summaryTable$diagdur.x, digits = 2)
round(summaryTable$age, digits = 2)
round(summaryTable$phga_BL, digits = 2)
round(summaryTable$anticcp, digits = 2)*100
round((dataWithBL %>% group_by(Z_value) %>% summarise(myper = mean(sex == "Female", na.rm = TRUE)))$myper, digits = 2)*100
round((dataWithBL %>% group_by(Z_value) %>% summarise(myper = mean(smoker == 2, na.rm = TRUE)))$myper, digits = 2)*100

CorrData <- dataWithBL %>% select(tcid, das28crp_BL, diagdur.x, phga_BL, age, sex, anticcp, smoker)
CorrData <- data %>% select(tcid, Z_value) %>% merge(CorrData, by = "tcid")

CorrData[paste0(TLevels, "Price")] <- NA
for (z in 1:10) {
  for (t in TLevels) {
    CorrData[CorrData$Z_value == z, paste0(t, "Price")] <- which(Z[[z]] == t)
  }
}
CorrData <- CorrData %>% select(-tcid, -Z_value)
CorrData$sex <- (CorrData$sex == "Female")
CorrData$smoker <- (CorrData$smoker == 2)
colnames(CorrData) <- c("DAS28 CRP", 
                       "Diagnosis duration",
                       "PhGA",
                       "Age",
                       "Sex",
                       "Anti-CCP positive",
                       "Current smoker",
                       "Inf price",
                       "Ada price",
                       "Cer price",
                       "Eta price",
                       "Gol price"
                       )
corrMat <- cor(CorrData, use = "pairwise.complete.obs", method = "spearman")
corrplot::corrplot(corrMat, method = "color", type = "upper", diag = FALSE, addCoef.col = "black")
