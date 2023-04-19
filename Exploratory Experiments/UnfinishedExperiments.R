#Read and organize the data
library(ggplot2)
library(ggpubr)
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

P_Z <- MakeP_Z(data, "Z_value", "trtgrp.x")



#Plot the adherence to LIS-anbudet
list_plots <- list()
for(counter in 1:length(Z)){
  plot_df <- data.frame(x = 1:5, y = as.vector(t(P_Z[counter,Z[[counter]]])))
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


#Plot a standard of care model
SC_summary <- data %>% group_by(Z_value) %>%
  summarise(mean_y = mean(das28crprem),
            mean_month = mean(cal_months))
SC_summary$NDPC_month <- as.numeric(difftime(NDPC_mid_date,
                                             as.Date("2010-01-31", "%Y-%m-%d"), 
                                             units = "days"))/30.417
SC_summary$Z_group <- c("Not repeated","Not repeated","Inf Eta","Not repeated",
                        "Inf Cer","Inf Cer","Inf Cer", "Inf Eta", "Inf Eta",
                        "Not repeated")

ggplot() + geom_smooth(data = data, aes(x = cal_months, y = das28crprem), method = "glm", formula = y ~ x, method.args = list(family = "binomial")) + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y)) + ylab("Average Remission rate") + xlab("Months since Jan 2010")



ggplot() + geom_smooth(data = data, aes(x = cal_months, y = das28crprem), method = "glm", formula = y ~ x, method.args = list(family = "binomial")) + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y, colour = as.factor(Z_group))) + 
  scale_colour_manual(values = c("green", "red", "black"), name = "Paired Years") + ylab("Average Remission rate") + xlab("Days since the first NDPC recommendation")


ggplot() + geom_smooth(data = data, aes(x = cal_months, y = das28crprem, colour = as.factor(Z_value + 2009)), method = "glm", formula = y ~ x, method.args = list(family = "binomial")) + 
  scale_colour_hue(name = "NDPC Period") + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y)) + ylab("Average Remission rate") + xlab("Months since the first NDPC recommendation")

Adj_data <- data %>% filter(Z_value %in% c(3,5,6,7,8,9)) %>% 
  mutate(Z_group = ifelse(Z_value %in% c(3,8,9), "InfEta", "InfCer")) %>% 
  filter(!((Z_value %in% c(3,8,9)) & (trtgrp.x %in% c("Gol", "Ada", "Cer")))) %>% 
  filter(!((Z_value %in% c(5,6,7)) & (trtgrp.x %in% c("Eta", "Ada", "Gol")))) %>%
  mutate(S_F = as.factor(paste(Z_group, trtgrp.x, sep = "_"))) %>%
  merge(SC_summary, all.x = TRUE, by = "Z_value")

v1 <- c()
v2 <- summary(Adj_data$S_F)/sum(summary(Adj_data$S_F))
v3 <- c()
Adj_models <- list()
for(s in levels(Adj_data$S_F)){
  Adj_models[[s]] <- glm(das28crprem ~ NDPC_month, data = Adj_data %>% filter(S_F == s), family = "binomial")
  v1 <- c(v1, Adj_models[[s]]$coefficients["NDPC_month"])
  v3 <- c(v3, summary(Adj_models[[s]])$cov.scaled["NDPC_month", "NDPC_month"])
}

slope <- sum(v1*v2)
slopeVar <- sum(v2*v2*v3)
intercept <- -0.27

data <- data %>% mutate(Adjustment = boot::inv.logit(slope*cal_months))

myFun <- function(months)
  boot::inv.logit(slope*months) + intercept

ggplot() + geom_function(fun = myFun, data = data, aes(x = cal_months)) + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y, colour = as.factor(Z_group))) +  
  scale_colour_manual(values = c("green", "red", "black"), name = "Similar periods") + ylab("Average Remission rate") + xlab("Months since the first NDPC recommendation")



data <- data %>% mutate(Adjusted_y = das28crprem - Adjustment) 

#Calculate the observables:
Z_weights <- data %>% group_by(Z_value) %>% summarise(N = n())
Z_weights$Prob <- Z_weights$N/sum(Z_weights$N)
P_Z <- MakeP_Z(data, "Z_value", "trtgrp.x")
Q_Z <- MakeQ_Z(data, "Z_value", "trtgrp.x", "Adjusted_y")

#Find the identifiables:
nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)


#Identify the effects
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

#Make confidence intervals
BSCI <- BSCICalculator(5000, data, "Z_value", "trtgrp.x", "Adjusted_y", TLevels, Z, Cap = TRUE)
myCIs <- BSCI$CIs

#Identify crude effects
Q_Z_crude <- MakeQ_Z(data, "Z_value", "trtgrp.x", "das28crprem")
LATEs_crude <- LATEIdentifier(Q_Z_crude, KB, b, P_Sigma)

#Make confidence intervals for crude
BSCI_crude <- BSCICalculator(5000, data, "Z_value", "trtgrp.x", "das28crprem", TLevels, Z, Cap = TRUE)
myCIs_crude <- BSCI_crude$CIs



#Test out adjustment with z weights
P2_Z <- P_Z
Q2_Z <- Q_Z
for(t in TLevels){
  P2_Z[[t]] <- P_Z[[t]] * Z_weights$Prob
  Q2_Z[[t]] <- Q_Z[[t]] * Z_weights$Prob
}
KB2 <- MakeKB(R, TLevels, 4, balanced = TRUE, weights = Z_weights$Prob)
b2 <- KbSolver(KB2, 3)
P_Sigma2 <- P_SigmaIdentifier(P2_Z, KB2, b2)
LATEs2 <- LATEIdentifier(Q2_Z, KB2, b2, P_Sigma2)

#Make confidence intervals
BSCI2 <- BSCICalculator(1000, data, "Z_value", "trtgrp.x", "Adjusted_y", TLevels, Z, Cap = TRUE, balanced = TRUE)
myCIs2 <- BSCI2$CIs




#Assumption violation?
print(b$Inf_Eta %*% KB[["Inf"]]$B_t_i)
print(b$Inf_Eta %*% KB[["Eta"]]$B_t_i)

print(b$Ada_Eta %*% KB[["Ada"]]$B_t_i)
print(b$Ada_Eta %*% KB[["Eta"]]$B_t_i)

#CODE AFTER THIS NOT IN USE
Z <- list("2010" = c("a", "b"),
          "2011" = c("b", "a"),
          "2012" = c("a", "b")
)
TLevels <- c("a", "b")
nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4, balanced = TRUE, weights = c(1/5,3/5, 2/5))
b <- KbSolver(KB, 3)
b$a_b %*% KB$a$B_t_i

#Some playing
data$X <- rbinom(930,1,0.4)
L_Z <- MakeQ_Z(data, "Z_value", "trtgrp.x", "X")
pp <- PseudoPopulator("Ada_Gol",data, KB, b, P_Sigma,
                      "Z_value", "trtgrp.x", "Adjusted_y",4)

mean((pp$X*pp$w)[pp$trtgrp.x == "Ada"])
mean((pp$X*pp$w)[pp$trtgrp.x == "Gol"])
sum(pp$w[pp$trtgrp.x == "Ada"])/sum(pp$w)
sum(pp$w[pp$trtgrp.x == "Gol"])/sum(pp$w)
b$Ada_Gol %*% KB$Ada$B_t_i %*% L_Z$Ada/P_Sigma$E1$Ada_Gol
b$Ada_Gol %*% KB$Gol$B_t_i %*% L_Z$Gol/P_Sigma$E2$Ada_Gol



Adj_data$Z_group.x <- as.factor(Adj_data$Z_group.x)
v1 <- c()
v2 <- summary(Adj_data$Z_group.x)/sum(summary(Adj_data$Z_group.x))
Adj_models <- list()
for(s in levels(Adj_data$Z_group.x)){
  Adj_models[[s]] <- glm(das28crprem ~ mean_calc_days, data = Adj_data %>% filter(Z_group.x == s), family = "binomial")
  v1 <- c(v1, Adj_models[[s]]$coefficients["mean_calc_days"])
}

slope <- sum(v1*v2)

Adj_data_sum <- Adj_data %>% group_by(S_F, as.factor(Z_value)) %>% summarise(y_mean = mean(das28crprem),
                                                                             cal_days = mean(mean_calc_days))
ggplot() + geom_smooth(data = Adj_data, aes(x = mean_calc_days, y = das28crprem, colour = S_F), method = "glm", formula = y ~ x, se = FALSE, method.args = list(family = "binomial")) + 
  scale_colour_hue(name = "NDPC Period") + 
  geom_point(data = Adj_data_sum, aes(x = cal_days, y = y_mean, colour = S_F)) + ylab("Average Remission rate") + xlab("Days since the first NDPC recommendation")


#Experiments from the simulation

sim_data <- CatSimulator2(1000000, Z, TLevels, VT,
                          VP, VY, TY, Vadj)
P_Z <- MakeP_Z(sim_data, "Z", "T")
Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = FALSE)

SC_summary <- sim_data %>% group_by(Z) %>% summarise(mean_y = mean(Y))
ggplot() + geom_smooth(data = sim_data, aes(x = Z, y = Y), method = lm) + 
  geom_point(data = SC_summary, aes(x = Z, y = mean_y)) + ylab("Average Remission rate") +
  scale_x_continuous(breaks = 2010:2019)

sim_data2 <- sim_data
sim_data2$m <- sample.int(12, size = 1000000, replace = TRUE)
sim_data2 <- sim_data2 %>% mutate(Y = Y + ((Z*12) + m))
P_Z2 <- MakeP_Z(sim_data2, "Z", "T")
Q_Z2 <- MakeQ_Z(sim_data2, "Z", "T", "Y")
P_Sigma2 <- P_SigmaIdentifier(P_Z2, KB, b)
temp2 <- LATEIdentifier(Q_Z2, KB, b, P_Sigma2, RR = FALSE, AverageProb = FALSE)

SC_summary <- sim_data2 %>% group_by(Z) %>% summarise(mean_y = mean(Y))
ggplot() + geom_smooth(data = sim_data2, aes(x = Z, y = Y), method = lm) + 
  geom_point(data = SC_summary, aes(x = Z, y = mean_y)) + ylab("Average Remission rate") +
  scale_x_continuous(breaks = 2010:2019)

SC_model <- lm(Y ~ Z, data = sim_data2)
sim_data2$Adjusted_y <- round(SC_model$residuals)
P_Z3 <- MakeP_Z(sim_data2, "Z", "T")
Q_Z3 <- MakeQ_Z(sim_data2, "Z", "T", "Adjusted_y")
P_Sigma3 <- P_SigmaIdentifier(P_Z3, KB, b)
temp3 <- LATEIdentifier(Q_Z3, KB, b, P_Sigma3, RR = FALSE, AverageProb = FALSE)

SC_model <- lm(Y ~ m + as.factor(Z), data = sim_data2)
sim_data2 <- sim_data2 %>% mutate(Adjustment = (Z*12+m)*SC_model$coefficients["m"])
sim_data2$Adjusted_yZ <- round(sim_data2$Y - sim_data2$Adjustment)
P_Z4 <- MakeP_Z(sim_data2, "Z", "T")
Q_Z4 <- MakeQ_Z(sim_data2, "Z", "T", "Adjusted_yZ")
P_Sigma4 <- P_SigmaIdentifier(P_Z4, KB, b)
temp4 <- LATEIdentifier(Q_Z4, KB, b, P_Sigma4, RR = FALSE, AverageProb = FALSE)


#Create a confounder variable and specify how it affects the probabilities of A and Y
VLevels <- c(0,1)
VP <- c(0.2,0.8)
VA2prT <- list("Inf" = 16,"Ada" = 8,"Cer" = 4,"Eta" = 2, "Gol" = 1)
VA2 <- unlist(lapply(A, function(x){
  out <- 0
  for(t in x){
    out <- out + VA2prT[[t]]
  }
  return(out/length(x))
}))
VA1 <- max(unlist(VA2prT)) + 1 - VA2
VA1 <- c(rep(0, length(A)-1),1)
VA <- list(VA1/sum(VA1),
           VA2/sum(VA2))
VY <- c(0,0.25)
names(VA) <- VLevels


VT = as.matrix(data.frame(V1T = c(1,0,0,0,0), V2T = c(0,0,1,0,0)))
VP = c(0.3,0.3)
VtoZstrength = 0.5
VY = c(0.2,-0.1)
TY = c(0.1,0.25,0.4,0.55,0.7)
sim_data <- CatSimulator2(50000, Z, TLevels, VT, VP, VtoZstrength,VY, TY, ZT = "exp")
P_Z <- MakeP_Z(sim_data, "Z", "T")
Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = FALSE)





box_data <- data.frame(TT = names(unlist(sim_res[[1]])),
                       LATE = unlist(sim_res[[1]]),
                       nP = 1000)
for(nP in 2:10){
  box_data <- rbind(box_data,
                    data.frame(TT = names(unlist(sim_res[[nP]])),
                               LATE = unlist(sim_res[[nP]]),
                               nP = nP*1000))
}

box_data$TT <- as.factor(box_data$TT)
levels(box_data$TT) <- c("Ada-Eta", "Cer-Eta", "Cer-Gol", "Eta-Gol", "Inf-Cer", "Inf-Eta")
box_data$TT <- as.character(box_data$TT)
library(ggplot2)
ggplot(box_data %>% filter(abs(LATE) < 1), aes(x = TT, y = LATE, fill = as.factor(nP))) +
  geom_boxplot(outlier.shape = NA)  +
  facet_wrap(~as.factor(nP), scale="free_y") + 
  xlab("Treatment-Treatment") + scale_fill_discrete(name = "Number of patients")

hist(unlist(sim_res[[4]])[names(unlist(sim_res[[4]])) == "a_b"])

mean_var_data <- box_data %>% group_by(nP, TT) %>% summarise(Mean = mean(LATE),
                                                             STDE = sd(LATE)) %>% ungroup()

ggplot(mean_var_data %>% filter(TT != "Ada-Eta"), aes(x = nP/1000, y = STDE, color = as.factor(TT))) + geom_line() +
  scale_x_continuous(breaks = 1:40)


data <- CatSimulator(n = 10000,
                     Z = Z,
                     A = A,
                     TLevels = TLevels,
                     VLevels = c(0,1),
                     VP = VP,
                     VA = VA,
                     VY = VY,
                     TY = TY,
                     T_decider = T_decider)

P_Z <- MakeP_Z(data, "Z", "T")
Q_Z <- MakeQ_Z(data, "Z", "T", "Y")
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

pData <-  PseudoPopulator(tt = "Inf_Cer", KB = KB, b = b, data = sim_data,
                          P_Z = P_Z, Z_column = "Z", T_column = "T", tolerance = 4)

pData$new_Y <- pData$Y*pData$w

den <- (sum(pData$w[(pData$T == "Inf") & (pData$new_Z == "Inf")])/
          sum(pData$w[pData$new_Z == "Inf"])) - 
  (sum(pData$w[(pData$T == "Inf") & (pData$new_Z == "Cer")])/
     sum(pData$w[pData$new_Z == "Cer"]))

late <- (mean(pData$new_Y[pData$new_Z == "Inf"]) - 
           mean(pData$new_Y[pData$new_Z == "Cer"]))/den


cl <- makeCluster(13)
registerDoParallel(cl)
idPrNz <- foreach(nZ = 2:20,
                  .combine = cbind,
                  .packages = c("dplyr", "tidyverse", "haven", "combinat")) %dopar% {
                    idNrVec <- c()
                    for(counter in 1:50){
                      Z <- all_instruments[sample.int(length(all_instruments), nZ)]
                      names(Z) <- 1:nZ
                      A <- GenerateA(TLevels)
                      R <- MakeR(A, Z, T_decider)
                      KB <- MakeKB(R, TLevels, 3)
                      b <- KbSolver(KB, 2)
                      commonSigma <- SigmaIdentifier(b)
                      idNrVec <- c(idNrVec, length(commonSigma)-sum(unlist(lapply(commonSigma, is_empty))))
                    }
                    idNrVec
                  }
stopCluster(cl)

colnames(idPrNz) <- 2:20
idPrNz_w <- reshape2::melt(idPrNz)
lm(value ~ Var2, data = idPrNz_w) %>% summary()

v <- idPrNz_w %>% group_by(Var2) %>% summarise(v = var(value))

plot(v$Var2, v$v)



idPrNz <- list()
for(nZ in 2:30){
  idNrVec <- c()
  print(nZ)
  for(counter in 1:20){
    print(counter)
    Z <- all_instruments[sample.int(length(all_instruments), nZ)]
    names(Z) <- 1:nZ
    A <- GenerateA(TLevels)
    R <- MakeR(A, Z, T_decider)
    KB <- MakeKB(R, TLevels, 3)
    save.image("myim.Rdata")
    b <- KbSolver(KB, 2)
    Sigma <- SigmaIdentifier(b)
    commonSigma <- CommonSigmaIdentifier(Sigma)
    idNrVec <- c(idNrVec, length(commonSigma)-sum(unlist(lapply(commonSigma, is_empty))))
  }
  idPrNz[[nZ]] <- idNrVec
}


mean(sim_data$Y[sim_data$T == "a"])- mean(sim_data$Y[sim_data$T == "c"])

sim_data2 <- sim_data %>% filter(T != "c")
(mean(sim_data2$Y[sim_data2$Z == 1]) - mean(sim_data2$Y[sim_data2$Z == 2]))/
  (mean((sim_data2$T[sim_data2$Z == 1]) == "a") - mean((sim_data2$T[sim_data2$Z == 2]) == "a"))


#Specify the treatments and the instrument values
TLevels <- c("a", "b", "c")
Z <- list("1" = c("a", "c", "b"),
          "2" = c("b", "c", "a"),
          "3" = c("b", "a", "c"),
          "4" = c("c", "a", "b")
)
names(Z) <- 1:nZ
#Generate all the adherence sets
A <- GenerateA(TLevels)
#Find the identifiable effects
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)



all_instruments <- permn(TLevels)
for(extrainstrument in all_instruments){
  Z[["3"]] <- extrainstrument
  names(Z) <- 1:nZ
  #Generate all the adherence sets
  A <- GenerateA(TLevels)
  #Find the identifiable effects
  R <- MakeR(A, Z, T_decider)
  KB <- MakeKB(R, TLevels, 4)
  b <- KbSolver(KB, 3)
  identified <- FALSE
  used <- FALSE
  identified <- (length(b$a_b) > 0)
  if(identified){
    w <- abs(b$a_b %*% KB$a$B_t_i) + abs(b$a_b %*% KB$b$B_t_i)
    used <- ((w[1] > 0.0001) & (w[2] > 0.0001))
  }
  if(used){
    print(extrainstrument)
  }
}

R1 <- MakeR(A1, Z, T_decider)
KB1 <- MakeKB(R1, TLevels, 4)
b1 <- KbSolver(KB1, 3)
Pis1 <- PiIdentifier(b1)

