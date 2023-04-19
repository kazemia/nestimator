#Final paper analysis
source("Functions.R")
source("stata.R")
data <- read.csv("data/adanon.csv")
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
data$Year = as.integer(2009 + data$Z_value)
data$trtgrp <- factor(data$trtgrp, levels = TLevels)
data$Z_value <- factor(data$Z_value, levels = 1:10)

adeff <- data
remove(data)

pz_margins <- function(data = adeff) {
  x <- glue::glue(
    "
tempfile tmp

quietly mlogit  trtgrp i.Z_value das28crp_bl
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

quietly logistic das28crprem i.Z_value#i.trtgrp das28crp_bl
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

#Adjusted analysis
P_Z <- pz_margins(adeff)
ptzm_covar <- ptzm_margins(adeff)
ptzm_covar <- ptzm_covar[colnames(P_Z)]
#Fix the missing values of Q_Z
#Identify missing values
missingCoords <- which(is.na(ptzm_covar))
zs <- as.character(missingCoords%%10)
zs[zs == "0"] <- "10"
t_i <- as.integer(missingCoords/10)
t_i[zs == "10"] <- t_i[zs == "10"] - 1
ts <- TLevels[t_i]
#Impute with the stratified average remission rate whereever available
for (counter in 1:length(missingCoords)) {
  ptzm_covar[as.character(ptzm_covar$Z_value) == zs[counter], t_i[counter]+1] <- 
               mean(adeff$das28crprem[(adeff$Z_value == zs[counter]) & 
                                        (adeff$trtgrp == ts[counter])])
  
}
#Make Q_Z
Q_Z <- P_Z * ptzm_covar
Q_Z$Z_value <- 1:10
#Impute missing values by 0 because P_Z is 0
Q_Z[is.na(Q_Z)] <- 0
#Find the identifiables
nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

#Calculate adjusted P_Sigma and LATES
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

#Run the bootstrap
data <- adeff
P_Z_list <- list()
Q_Z_list <- list()
P_Sigma_list <- list()
LATEs_list <- list()
options(show.error.messages = FALSE)
counter2 <- 1
while (counter2 <= 2500) {
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
  counter2 <- counter2 + 1
}
options(show.error.messages = TRUE)

save.image(file = "data/ClinicalResults.Rdata")

#Print the adjusted confidence intervals and estimates
LATES_BS <- t(matrix(unlist(LATEs_list), nrow = 9))
round(apply(LATES_BS,2,median), digits = 2)
LATES_BS[(LATES_BS > 1)|(LATES_BS < -1)] <- NA
2500-colSums(is.na(LATES_BS))
round(apply(LATES_BS,2,quantile, probs = c(0.025,0.975), na.rm = TRUE), digits = 2)

#Print the crude confidence intervals and estimates
set.seed(123)
Crude_BS <- BSCICalculator(2500, data %>% mutate(Z_value = as.numeric(Z_value)), "Z_value", "trtgrp", "das28crprem", TLevels, Z)
LATES_BS <- Crude_BS$BSData
round(apply(LATES_BS,2,median), digits = 2)
LATES_BS[(LATES_BS > 1)|(LATES_BS < -1)] <- NA
2500-colSums(is.na(LATES_BS))
round(apply(LATES_BS,2,quantile, probs = c(0.025,0.975), na.rm = TRUE), digits = 2)

#Calculate P_Sigma crude
P_Z_crude <- MakeP_Z(data,Z_column = "Z_value", T_column = "trtgrp")
P_Z_crude <- P_Z_crude %>% mutate(Z_value = as.numeric(Z_value))
P_Sigma_crude <- P_SigmaIdentifier(P_Z_crude, KB, b)

#Print table 3
data %>% group_by(Z_value) %>% summarise(n())

#Make figure 5
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

#Make figure 6
SC_summary <- data %>% group_by(Z_value) %>%
  summarise(mean_y = mean(das28crprem),
            mean_month = mean(month))

ggplot() + geom_smooth(data = data, aes(x = month, y = das28crprem), method = "glm", formula = y ~ x, method.args = list(family = "binomial")) + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y)) + ylab("Average Remission rate") + xlab("Months since Jan 2010")

