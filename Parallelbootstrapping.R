#We do parallell bootstrapping here)
library(haven)
setwd("C:/Users/kazem/Documents/Data")
iv_price <- read_dta("iv_price.dta")
data <- read_dta("addisact.dta")
RA_iv_price <- iv_price %>% filter(diaggrp == 1)
SPA_iv_price <- iv_price %>% filter(diaggrp == 2)
PSA_iv_price <- iv_price %>% filter(diaggrp == 3)

#Look at only RA and rename the columns
iv_data <- RA_iv_price %>% select(lisperiod, pinx, pada, pcerto, petan, pgoli)
colnames(iv_data) <- c("Year", "Infliximab", "Adalimumab", "Certolizumab", "Etanercept", "Golimumab")
iv_data$Year <- iv_data$Year + 2008

#Take away 2009 for the lack of availabality of treatments
iv_data <- iv_data %>% filter(Year != 2009)
#Define number of treatments and number of instruments
nT <- ncol(iv_data) - 1
nZ <- nrow(iv_data)

#Create the list of instruments and the list of possible treatments
Z <- iv_data %>% select(-Year) %>% split(1:nZ) %>% lapply(function(x){names(x)[order(x)]})
names(Z) <- iv_data$Year
TLevels <- colnames(iv_data)[2:ncol(iv_data)]

#data cleaning

#Pick useful columns only
data <- data %>% select(tcid, trtgrp, diaggrp, visit_no, visit_dt, das28rem, timepoint, bldt) %>%
  filter(!is.na(timepoint))

#Convert to character
for(column in c("trtgrp", "diaggrp"))
  data[[column]] <- names(attributes(data[[column]])$labels)[
    match(data[[column]], attributes(data[[column]])$labels)]

#Remove those with missing DAS
data <- data %>% filter(!is.na(das28rem))

#Create a year and month column
data <- data %>% mutate(Year = as.numeric(format(bldt, format="%Y")))

#identify treatments in the study period
treatmentsInPeriod <- data$tcid[(data$Year >= min(iv_data$Year)) & 
                                  (data$Year <= max(iv_data$Year)) & 
                                  (data$timepoint == 1)]

#Remove those outside the period and only look at visits at 3 months
data <- data %>% filter(tcid %in% treatmentsInPeriod) %>% filter(timepoint == 2)

#Only look at RA
data <- data %>% filter(diaggrp == "RA")

#Remove irrelavant treatments and rename them
data <- data %>% mutate(trtgrp = substr(trtgrp, 1,3)) %>% filter(trtgrp %in% substr(TLevels, 1, 3))
data <- data %>% mutate(trtgrp = TLevels[match(trtgrp, substr(TLevels, 1, 3))])

remove(iv_data, iv_price, PSA_iv_price, RA_iv_price, SPA_iv_price)

BSCI <- BSCICalculator(1000, data, "Year", "trtgrp", "das28rem", TLevels, Z)

myCluster <- makeCluster(14)
registerDoParallel(myCluster)
runtime <- system.time({
  boot_b <- foreach(b=idiv(1000, chunks=getDoParWorkers()), .combine = "cbind",
                    .packages = c("dplyr", "tidyverse", "haven")) %dopar% {
                      sapply(1:b, StatisticCalculator, data,
                             "Year", "trtgrp", "das28rem", TLevels)
                    }
})
stopCluster(myCluster)
bootstrap_data <- t(boot_b)
bootstrap_data <- as.data.frame(bootstrap_data)

for(c in colnames(bootstrap_data)){
  print(c)
  print(quantile(bootstrap_data[[c]], c(0.1,0.9), na.rm = T))
}


