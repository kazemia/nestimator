library(dplyr)
library(tidyverse)
library(doParallel)
library(haven)
library(combinat)
library(boot)

#Read and organize the data
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



P_Z <- MakeP_Z(data, "Year", "trtgrp")

plot(1:5, P_Z[y,match(Z[[y]], colnames(P_Z))], type = "l")
lines(1:5, exp(5:1)/sum(exp(1:5)))

Q_Z <- MakeQ_Z(data, "Year", "trtgrp", "das28rem")
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 5)
b <- KbSolver(KB,4)
Pis <- PiIdentifier(b)
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
Q_Sigma <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

CIs <- BSCICalculator(1000, data, "Year", "trtgrp", "das28rem", TLevels, Z)

myCluster <- makeCluster(10)
bootstrap <- boot(data, StatisticCalculator, 1000,
                  Z_column = "Year", T_column = "trtgrp", Y_column = "das28rem",
                  parallel = "multicore",
                  cl = myCluster)
stopCluster(myCluster)
bootstrap_data <- bootstrap$t
bootstrap_data <- as.data.frame(bootstrap_data)
colnames(bootstrap_data) <- names(Contrasts)[sapply(Contrasts, function(x) !is_empty(x))]

for(c in colnames(bootstrap_data)){
  print(c)
  print(quantile(bootstrap_data[[c]], c(0.025,0.975), na.rm = T))
}

#---------------------------------------------------------------------------------------------------#

#Binary simulation


BinarySimulator <- function(ZT, VT, VY, TY, n){
  V <- rbinom(n,1,0.5)
  Z <- rbinom(n,1,0.5)
  T <- rbinom(rep(n,n),1,(VT*V) + (ZT*Z))
  Y <- rbinom(rep(n,n),1,(VY*V) + (TY*T))
  return(data.frame(Z = Z, T = T, Y = Y))
}
sim_data <- BinarySimulator(0.3,0.2,0.6,0.4,1000)
mean(sim_data$Y[sim_data$T == 1])- mean(sim_data$Y[sim_data$T == 0])
(mean(sim_data$Y[sim_data$Z == 1]) - mean(sim_data$Y[sim_data$Z == 0]))/
  (mean(sim_data$T[sim_data$Z == 1]) - mean(sim_data$T[sim_data$Z == 0]))




P_Z <- MakeP_Z(sim_data, "Z", "T")
Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
A <- GenerateA(unique(sim_data$T))
Z <- list(c(0,1), c(1,0))
names(Z) <- c(0,1)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, unique(sim_data$T), 0.0001)
b <- KbSolver(KB)
Sigma <- SigmaIdentifier(b)
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
Q_Sigma <- Q_SigmaIdentifier(Q_Z, KB, b, P_Sigma)
commonSigma <- CommonSigmaIdentifier(Sigma)
Contrasts <- ContrastFinder(Sigma, commonSigma, Q_Sigma)
ContrastStrength <- ContrastStrengthFinder(Sigma, commonSigma, P_Sigma)

#---------------------------------------------------------------------------------------------------#

#Multivariate simulation

#n total number of patients
#Z list of ordered treatments
#A list of adherence sets
#TLevels vector of treatments
#VLevels vector of levels of confounders
#VP vector of probability of being in each confounder level (sums to 1)
#VA list of vecotr of probability of A for each level of V (each element sums to 1)
#VY vector of effect of each level of V on Y
#TY vector of effects of T on Y
#max(VY) + max(TY) <= 1
#T_decider the choice function
CatSimulator <- function(n, Z, A, TLevels, VLevels, VP, VA, VY, TY, T_decider){
  #Make a confounder
  V_obs <- sample(VLevels, n, replace = TRUE, prob = VP)
  A_obs <- A[unlist(lapply(VA[as.character(V_obs)],
                             function(x) sample.int(length(A), size = 1, prob = x)))]
  names(A_obs) <- 1:n
  Z_obs_ind <- sample.int(length(Z), size = n, replace = TRUE)
  Z_obs <- Z[Z_obs_ind]
  names(Z_obs) <- 1:n
  T_obs <- sapply(1:n, function(x) T_decider(Z_obs[[x]], A_obs[[x]]))
  
  #Assign patients remission
  YP <- TY[match(T_obs, TLevels)] + VY[match(V_obs,VLevels)]
  Y_obs <- rbinom(rep(n,n),1,YP)
  return(data.frame(Z = Z_obs_ind, T = T_obs, Y = Y_obs))
}



TLevels <- c("a", "b", "c")
all_instruments <- permn(TLevels)
nZ <- 5
Z <- all_instruments[sample.int(length(all_instruments), nZ)]
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 0.0001)
b <- KbSolver(KB)
Sigma <- SigmaIdentifier(b)
commonSigma <- CommonSigmaIdentifier(Sigma)

VLevels <- c(0,1)
VP <- c(0.7,0.3)
VA <- list(c(0,0,0,0.25,0.25,0.25,0.25), c(0.1,0.3,0.3,0,0,0.3,0))
VY <- c(0.4,0.2)
TY <- c(0.2,0.4,0.5)
names(VA) <- VLevels

sim_data <- CatSimulator(n = 10000,
                         Z = Z,
                         A = A,
                         TLevels = TLevels,
                         VLevels = c(0,1),
                         VP = VP,
                         VA = VA,
                         VY = VY,
                         TY = TY,
                         T_decider = T_decider)

P_Z <- MakeP_Z(sim_data, "Z", "T")
Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
Q_Sigma <- Q_SigmaIdentifier(Q_Z, KB, b, P_Sigma)
Contrasts <- ContrastFinder(Sigma, commonSigma, Q_Sigma)
ContrastStrength <- ContrastStrengthFinder(Sigma, commonSigma, P_Sigma)

mean(sim_data$Y[sim_data$T == "a"])- mean(sim_data$Y[sim_data$T == "c"])

sim_data2 <- sim_data %>% filter(T != "c")
(mean(sim_data2$Y[sim_data2$Z == 1]) - mean(sim_data2$Y[sim_data2$Z == 2]))/
  (mean((sim_data2$T[sim_data2$Z == 1]) == "a") - mean((sim_data2$T[sim_data2$Z == 2]) == "a"))



