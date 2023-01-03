#Read and organize the data
setwd("C:/Users/kazem/Documents/Data")
TLevels <- c("Inf", "Ada", "Cer", "Eta", "Gol")
Z <- list("2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol") 
)
load("myData.rdata")
data <- as.data.frame(my_data)
remove(my_data)
data$Z_value <- NA
for(y in 2010:2013){
  data$Z_value[(data$bldt.x > as.Date(paste0(y, "-01-31"), "%Y-%m-%d")) & (data$bldt.x < as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d"))] <- y-2009
}
data$Z_value[(data$bldt.x > as.Date("2014-01-31", "%Y-%m-%d")) & (data$bldt.x < as.Date("2015-03-01", "%Y-%m-%d"))] <- 5
for(y in 2015:2016){
  data$Z_value[(data$bldt.x >= as.Date(paste0(y, "-03-01"), "%Y-%m-%d")) & (data$bldt.x < as.Date(paste0(y+1, "-03-01"), "%Y-%m-%d"))] <- y-2009
}
data$Z_value[(data$bldt.x >= as.Date("2017-03-01", "%Y-%m-%d")) & (data$bldt.x < as.Date("2018-02-01", "%Y-%m-%d"))] <- 8
for(y in 2018:2019){
  data$Z_value[(data$bldt.x > as.Date(paste0(y, "-01-31"), "%Y-%m-%d")) & (data$bldt.x < as.Date(paste0(y+1, "-02-01"), "%Y-%m-%d"))] <- y-2009
}

data <- data %>% filter(!is.na(Z_value))

nZ <- length(Z)
names(Z) <- 1:nZ
#Generate all the adherence sets
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

BSCI <- BSCICalculator(1000, data, "Z_value", "trtgrp.x", "das28crprem", TLevels, Z)

myCIs <- BSCI$CIs

