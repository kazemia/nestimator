#Final paper analysis
setwd("C:/Users/kazem/Dropbox/NESTIMATOR PhD/R Scripts/nestimator")
source("Functions.R")
source("stata.R")
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
data$das28crprem <- as.numeric(data$das28crprem)
data$Year = as.integer(2009 + data$Z_value)
data$cal_months <- as.numeric(difftime(data$bldt.x,
                                       as.Date("2010-01-31", "%Y-%m-%d"), 
                                       units = "days"))/30.417
data <- data %>% select(-bldt.x) %>% as.data.frame()

data <- data %>% select(-tcid, -patid.x)
data$trtgrp.x <- factor(data$trtgrp.x, levels = TLevels)
data$Z_value <- factor(data$Z_value, levels = 1:10)

imputationModel <- mice::mice(data, m = 1)
data <- mice::complete(imputationModel)
remove(imputationModel)

adeff <- data
remove(data)

adeff <- adeff %>% mutate(trtgrp = trtgrp.x) %>% select(-trtgrp.x, -visit_dt_BL)

pz_margins <- function(data = adeff) {
  x <- glue::glue(
    "
tempfile tmp

quietly mlogit  trtgrp i.Z_value das28crp_BL
quietly margins  i.Z_value, saving(`tmp')
use `tmp', clear

"
  )
  
  margins <-  stata(
    src = x,
    data.in = data,
    data.out = TRUE,
    stata.path = "\"C:\\Program Files\\Stata17\\StataSE-64\"",
    stata.version = 17,
    stata.echo = TRUE
  )
  
  res <- margins %>%
    rename_all(~ str_replace(., "_", "")) %>%
    mutate(predict = parse_number(as.character(predict)),
           trtgrp = case_when(
             predict == 1 ~ "Inf",
             predict == 2 ~ "Ada",
             predict == 3 ~ "Cer",
             predict == 4 ~ "Eta",
             predict == 5 ~ "Gol",
           )) %>%
    select(trtgrp , estimate = margin, Z_value = m1) %>%
    pivot_wider(names_from = trtgrp, values_from = estimate) %>%
    mutate(Z_value = as.numeric(Z_value))
  
  
  return(res)
  
}

ptzm_margins <- function(data = adeff) {
  x <- glue::glue(
    "
tempfile tmp

quietly logistic das28crprem i.Z_value#i.trtgrp das28crp_BL
quietly margins Z_value#trtgrp,  saving(`tmp')
use `tmp', clear

"
  )
  
  margins <-  stata(
    src = x,
    data.in = data,
    data.out = TRUE,
    stata.path = "\"C:\\Program Files\\Stata17\\StataSE-64\"",
    stata.version = 17,
    stata.echo = TRUE
  )
  
  res <- margins %>%
    rename_all(~ str_replace(., "_", "")) %>%
    select(estimate = margin, Z_value = m1, trtgrp = m2 ) %>%
    pivot_wider(names_from = trtgrp, values_from = estimate) %>% 
    mutate(Z_value = as.numeric(Z_value))
  
  
  return(res)
  
}


P_Z <- pz_margins(adeff)
ptzm_covar <- ptzm_margins(adeff)
ptzm_covar <- ptzm_covar[colnames(P_Z)]
missingCoords <- which(is.na(ptzm_covar))
zs <- as.character(missingCoords%%10)
zs[zs == "0"] <- "10"
t_i <- as.integer(missingCoords/10)
t_i[zs == "10"] <- t_i[zs == "10"] - 1
ts <- TLevels[t_i]
for (counter in 1:length(missingCoords)) {
  ptzm_covar[as.character(ptzm_covar$Z_value) == zs[counter], t_i[counter]+1] <- 
               mean(adeff$das28crprem[(adeff$Z_value == zs[counter]) & 
                                        (adeff$trtgrp == ts[counter])])
  
}
Q_Z <- P_Z * ptzm_covar
Q_Z$Z_value <- 1:10
Q_Z[is.na(Q_Z)] <- 0

nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

data <- adeff
P_Z_list <- list()
Q_Z_list <- list()
P_Sigma_list <- list()
LATEs_list <- list()
options(show.error.messages = FALSE)
set.seed(123)
for (counter2 in 1:2500) {
  adeff <- data[sample(nrow(data), nrow(data), replace=T),]
  P_Z <- try(pz_margins(adeff))
  if(is.character(P_Z)) next
  ptzm_covar <- try(ptzm_margins(adeff))
  if(is.character(ptzm_covar)) next
  ptzm_covar <- try(ptzm_covar[colnames(P_Z)])
  if(is.character(ptzm_covar)) next
  missingCoords <- which(is.na(ptzm_covar))
  zs <- as.character(missingCoords%%10)
  zs[zs == "0"] <- "10"
  t_i <- as.integer(missingCoords/10)
  t_i[zs == "10"] <- t_i[zs == "10"] - 1
  ts <- TLevels[t_i]
  for (counter in 1:length(missingCoords)) {
    ptzm_covar[as.character(ptzm_covar$Z_value) == zs[counter], t_i[counter]+1] <- 
      mean(adeff$das28crprem[(adeff$Z_value == zs[counter]) & 
                               (adeff$trtgrp == ts[counter])])
    
  }
  Q_Z <- try(P_Z * ptzm_covar)
  if(is.character(Q_Z)) next
  Q_Z$Z_value <- 1:10
  Q_Z[is.na(Q_Z)] <- 0
  P_Sigma <- try(P_SigmaIdentifier(P_Z, KB, b))
  if(is.character(P_Sigma)) next
  LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)
  P_Z_list[[counter2]] <- P_Z
  Q_Z_list[[counter2]] <- Q_Z
  P_Sigma_list[[counter2]] <- P_Sigma
  LATEs_list[[counter2]] <- LATEs
}
options(show.error.messages = TRUE)



LATES_BS <- t(matrix(unlist(LATEs_list), nrow = 9))
round(apply(LATES_BS,2,median), digits = 2)
LATES_BS[(LATES_BS > 1)|(LATES_BS < -1)] <- NA
2500-colSums(is.na(LATES_BS))
round(apply(LATES_BS,2,quantile, probs = c(0.025,0.975), na.rm = TRUE), digits = 2)

set.seed(123)
Crude_BS <- BSCICalculator(2500, data %>% mutate(Z_value = as.numeric(Z_value)), "Z_value", "trtgrp", "das28crprem", TLevels, Z)
LATES_BS <- Crude_BS$BSData
round(apply(LATES_BS,2,median), digits = 2)
LATES_BS[(LATES_BS > 1)|(LATES_BS < -1)] <- NA
2500-colSums(is.na(LATES_BS))
round(apply(LATES_BS,2,quantile, probs = c(0.025,0.975), na.rm = TRUE), digits = 2)


P_Z_crude <- MakeP_Z(data,Z_column = "Z_value", T_column = "trtgrp")
P_Z_crude <- P_Z_crude %>% mutate(Z_value = as.numeric(Z_value))
P_Sigma_crude <- P_SigmaIdentifier(P_Z_crude, KB, b)

data %>% group_by(Z_value) %>% summarise(n())

library(ggplot2)
library(ggpubr)
list_plots <- list()
for(counter in 1:length(Z)){
  plot_df <- data.frame(x = 1:5, y = as.vector(t(P_Z_crude[counter,Z[[counter]]])))
  list_plots[[counter]] <- (ggplot(plot_df, aes(x = x,y = y, group = 1)) + 
                              geom_col() +
                              scale_x_continuous(breaks=1:5,
                                                 labels = Z[[counter]]) + 
                              scale_y_continuous(breaks = (1:9)/10,
                                                 limits = c(0, 0.9)) +
                              ylab("Observed probability") + 
                              xlab("Medications ordered by price") +
                              ggtitle(paste("NDPC period", 2009 + counter)))
}
for(counter in c(2:5,7:10)){
  list_plots[[counter]] <- list_plots[[counter]] + 
    rremove("ylab") + rremove("y.text")
}
for(counter in 1:5)
  list_plots[[counter]] <- list_plots[[counter]] + rremove("xlab")
ggarrange(plotlist = list_plots, nrow = 2, ncol = 5, common.legend = TRUE)








model <- glm(das28crprem ~ Z_value:trtgrp + das28crp_BL-1, family = "binomial", data = adeff)
s <- margins::margins_summary(model, variables = "Z_value:trtgrp")
b <- summary(s)
R_Q_Z <- model$coefficients
R_Q_Z <- R_Q_Z[!(names(R_Q_Z) %in% c("das28crp_BL", "(Intercept)"))]
R_Q_Z <- matrix(boot::inv.logit(R_Q_Z), ncol = 5)
R_Q_Z <- as.data.frame(R_Q_Z)
colnames(R_Q_Z) <- TLevels
R_Q_Z$Z_value <- 1:10
R_Q_Z <- R_Q_Z[colnames(P_Z)]

R_Q_Z <- P_Z * R_Q_Z
R_Q_Z$Z_value <- 1:10
Q_Z[is.na(Q_Z)] <- R_Q_Z[is.na(Q_Z)]
Q_Z[is.na(Q_Z)] <- 0


R_Q_Z[is.na(R_Q_Z)] <- 0

Q_Z[is.na(Q_Z)] <- 0
plot(1:10, rowSums(Q_Z %>% select(-Z_value)))
#plot(1:10, rowSums(R_Q_Z %>% select(-Z_value)))
