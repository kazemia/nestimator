#Simulator with time confounding
#Specify the treatments and the instrument values
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
nZ <- length(Z)
names(Z) <- 1:nZ
#Generate all the adherence sets
A <- GenerateA(TLevels)
#A <- A[sort(c(27, 30, 3, 1, 6, 18, 4, 17, 14, 31))]
#Find the identifiable effects
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)


#Specify arguments again for CatSimulator2
VT = as.matrix(data.frame(V1T = c(0,0,0,exp(4),0), V2T = c(0,0,exp(3),0,0)))
VP = c(0.5,0.3)
VY = c(0.2, 0.2)
TY = c(0.05,0.15,0.25,0.35,0.45)
Vadj = c(TRUE, FALSE)

sim_data <- CatSimulator2(10000, Z, TLevels, VT,
                          VP, VY, TY, Vadj, YPReturn = TRUE)
sim_data <- sim_data %>% arrange(Z)
slope <- 0.06
sim_data$month <- (1:nrow(sim_data))*12*length(Z)/nrow(sim_data)
sim_data <- sim_data %>% mutate(
  TimeEffect = slope*month - (Z-1)*0.6,
  new_YP = boot::inv.logit(boot::logit(YP) + TimeEffect)
)
sim_data$new_Y <- rbinom(rep(10000,10000),1,matrix(sim_data$new_YP, ncol = 1))

data <- sim_data
colnames(data) <- c("Z_value", "trtgrp.x", "Y", "V1", "YP", "cal_months", "TimeEffect", "new_YP", "das28crprem")

NDPC_mid_date <- as.Date(paste0(2009, "-07-31"), "%Y-%m-%d") + (1:length(Z))*365


P_Z <- MakeP_Z(data, "Z_value", "trtgrp.x")

#Plot the adherence to LIS-anbudet
library(ggpubr)
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
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y, colour = as.factor(Z_group))) + 
  scale_colour_manual(values = c("green", "red", "black"), name = "Paired Years") + ylab("Average Remission rate") + xlab("Days since the first NDPC recommendation")


ggplot() + geom_smooth(data = data, aes(x = cal_months, y = das28crprem, colour = as.factor(Z_value + 2009)), method = "glm", formula = y ~ x, method.args = list(family = "binomial")) + 
  scale_colour_hue(name = "NDPC Period") + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y)) + ylab("Average Remission rate") + xlab("Months since the first NDPC recommendation")

#Adjust for standard of care by comparing similar periods

#Create a data set with only similar periods
Adj_data <- data %>% filter(Z_value %in% c(3,5,6,7,8,9)) %>% 
  mutate(Z_group = ifelse(Z_value %in% c(3,8,9), "InfEta", "InfCer")) %>% 
  filter(!((Z_value %in% c(3,8,9)) & (trtgrp.x %in% c("Gol", "Ada", "Cer")))) %>% 
  filter(!((Z_value %in% c(5,6,7)) & (trtgrp.x %in% c("Eta", "Ada", "Gol")))) %>%
  mutate(S_F = as.factor(paste(Z_group, trtgrp.x, sep = "_"))) %>%
  merge(SC_summary, all.x = TRUE, by = "Z_value")

#Within each group and treatment create a log reg

#Merge the log regs and add an intercept for plotting purposes
intercept <- -0.04

adj_m <- glm(das28crprem ~ S_F + NDPC_month:S_F, data = Adj_data, family = "binomial")
marginalEffects <- summary(margins::margins(adj_m, type = "link"))
slope <- marginalEffects$AME["NDPC_month"]
slopeVar <- marginalEffects$SE["Var_dydx_NDPC_month"]^2

#Calculate the adjustment needed
data <- data %>% mutate(Adjustment = boot::inv.logit(slope*cal_months))

#Plot the data and the model
myFun <- function(months)
  boot::inv.logit(slope*months) + intercept

ggplot() + geom_function(fun = myFun, data = data, aes(x = cal_months)) + 
  geom_point(data = SC_summary, aes(x = mean_month, y = mean_y, colour = as.factor(Z_group))) +  
  scale_colour_manual(values = c("green", "red", "black"), name = "Similar periods") + ylab("Average Remission rate") + xlab("Months since the first NDPC recommendation")

#Adjust the outcome
data <- data %>% mutate(Adjusted_y = das28crprem - Adjustment)



#Calculate the basic observables:
Z_weights <- data %>% group_by(Z_value) %>% summarise(N = n())
Z_weights$Prob <- Z_weights$N/sum(Z_weights$N)
P_Z <- MakeP_Z(data, "Z_value", "trtgrp.x")

#Find the identifiables:
nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

#Estimate probability of the subpopulations
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)

#Make the different Q_Zs:

#Q_Z adjusted with similar years:
Q_Z_adjusted <- MakeQ_Z(data, P_Z, "Z_value", "trtgrp.x", "Adjusted_y")

#Q_Z without adjustment
Q_Z_crude <- MakeQ_Z(data, P_Z, "Z_value", "trtgrp.x", "das28crprem")

#Q_Z adjusted with modelling

P_Z_adjusted <- MakeP_Z(data,  "Z_value", "trtgrp.x", Adj_columns = c("cal_months"), adj_method = "binomial")
Q_Z_adjusted1 <- MakeQ_Z(data, P_Z, "Z_value", "trtgrp.x", "das28crprem", Adj_columns = c("cal_months"), adj_method = "binomial")
Q_Z_adjusted2 <- MakeQ_Z(data, P_Z_adjusted, "Z_value", "trtgrp.x", "das28crprem", Adj_columns = c("cal_months"), adj_method = "binomial")

Q_Z_adjusted3 <- MakeQ_Z(data, P_Z, "Z_value", "trtgrp.x", "das28crprem", offset_columns = c("cal_months"), offset_constraints = c(0.0068) , adj_method = "binomial")

Q_Z_true <- MakeQ_Z(data, P_Z, "Z_value", "trtgrp.x", "Y")


sum((Q_Z_crude-Q_Z_true)^2)
sum((Q_Z_adjusted-Q_Z_true)^2)
sum((Q_Z_adjusted1-Q_Z_true)^2)
sum((Q_Z_adjusted2-Q_Z_true)^2)
sum((Q_Z_adjusted3-Q_Z_true)^2)

plot(1:10, rowSums(Q_Z_true)-1:10)
plot(1:10, rowSums(Q_Z_crude)-1:10)
plot(1:10, rowSums(Q_Z_adjusted)-1:10)
plot(1:10, rowSums(Q_Z_adjusted1)-1:10)
plot(1:10, rowSums(Q_Z_adjusted3)-1:10)
P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)

