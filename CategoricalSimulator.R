#Specify the treatments and the instrument values
source("Functions.R")
library(doRNG)
TLevels <- c("t1", "t2", "t3")
Z <- list("2010" = c("t1", "t2", "t3"),
          "2011" = c("t2", "t1", "t3"),
          "2012" = c("t2", "t3", "t1"),
          "2013" = c("t3", "t2", "t1")
)

nZ <- length(Z)
names(Z) <- 1:nZ
#Generate all the adherence sets
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)


#Specify arguments again for CatSimulator2
VT = as.matrix(data.frame(V1T = c(exp(3),0,0), V2T = c(0,exp(3),0)))
VP = c(0.5,0.5)
VY = c(0.3, 0.3)
TY = c(0.05,0.2,0.35)
Vadj = c(TRUE, FALSE)

#Specify maximum number of observations in simulation, the step 
#and number of simulations pr step
sim_max <- 5000
sim_start <- 500
sim_step <- 500
n_sim <- 1000
TT_vec <- c("t1_t2", "t1_t3", "t2_t3")

#Simulate
myCluster <- makeCluster(15)
registerDoParallel(myCluster)
registerDoRNG(123)
sim_res <- list()
for(nP in floor(sim_start/sim_step):floor(sim_max/sim_step)){
  sim_res[[nP]] <- foreach(counter=idiv(n_sim, chunks=getDoParWorkers()), .combine = rbind,
                           .packages = c("dplyr", "tidyverse", "stringr","lpSolve")
  ) %dopar% {
    out <- data.frame(TT = NA, n = NA, method = NA, estimate = NA)
    for(counter2 in 1:counter){
      #Simulate the data
      sim_data <- CatSimulator2(nP*sim_step, Z, TLevels, VT,
                                VP, VY, TY, Vadj)
      #Calculate the CIV estimator
      P_Z <- MakeP_Z(sim_data, "Z", "T")
      Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
      P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
      temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma)
      out <- rbind(out, data.frame(
        TT = TT_vec,
        n = rep(nP*sim_step, length(TT_vec)),
        method = rep("LATE_IV_unadj", length(TT_vec)),
        estimate = unlist(temp, use.names = FALSE)
      ))
      #Calculate the DIV and the NAIVE estimators
      for(tt in names(temp)){
        out <- rbind(out, data.frame(
          TT = rep(tt,2),
          n = rep(nP*sim_step, 2),
          method = c("ATE_naive", "IV_naive"),
          estimate = c(mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][1]]) - 
                         mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][2]]),
                       NaiveIV(tt, sim_data, Z, "T", "Z", "Y"))
        ))
      }
    }
    out[2:nrow(out),]
  }
}
stopCluster(myCluster)

#Put the results into a dataframe
sim_results <- sim_res[[1]]
for (counter in 2:length(sim_res)) {
  sim_results <- rbind(sim_results, sim_res[[counter]])
}
save.image(file = "data/FinalSimulationResults.Rdata")
sim_results_backup <-sim_results
remove(sim_res)
#Cap the effects
#sim_results$estimate[sim_results$estimate < (-1)] <- (-1)
#sim_results$estimate[sim_results$estimate > 1] <- 1
sim_results <- sim_results %>% filter(abs(estimate) <= 1)
#Summarise the simulation results
sim_results <- sim_results %>% group_by(TT) %>%
  mutate(bias = estimate - TY[TLevels == str_split(TT, "_")[[1]][1]] + 
           TY[TLevels == substr(str_split(TT, "_")[[1]][2], 1,3)])

my_sum <- sim_results %>%
  group_by(TT, n, method) %>% summarise(est_bias = median(bias, na.rm = TRUE),
                                        est_sd = log(sd(bias, na.rm = TRUE)),
                                        est_lq = quantile(bias, 0.025 , na.rm = TRUE),
                                        est_hq = quantile(bias, 0.975 , na.rm = TRUE),
                                        est_mean = median(estimate))



#my_sum <- my_sum %>% filter(n>999)
TTs <- unique(my_sum$TT)
meths <- c("ATE_naive", "IV_naive", "LATE_IV_unadj")
#Plot the results
library(ggplot2)
library(gridExtra)
library(latex2exp)
list_plots <- list()
counter <- 1
for(t in TTs){
  t1 <- which(TLevels == str_split(t, "_")[[1]][1])
  t2 <- which(TLevels == substr(str_split(t, "_")[[1]][2],1,3))
  title <- TeX(paste0(r'($\hat{E}(Y(t_)', t1, r'()-Y(t_)', t2, r'())$)' ))
  list_plots[[counter]] <- (ggplot(my_sum %>%
                                     filter((method %in% meths) & (TT == t)),
                                   aes(x = n,
                                       y = est_mean,
                                       color = method, 
                                       linetype = method)) + 
                              geom_line() + 
                              geom_hline(yintercept = (TY[t1] - TY[t2])) +
                              scale_x_continuous(n.breaks = 5) + 
                              scale_y_continuous(n.breaks = 5) +
                              ylab("Median") + 
                              xlab(TeX(r'($n_s$)')) +
                              scale_color_discrete(name = "Estimator:", labels = c("Naive", "DIV", "CIV")) + 
                              scale_linetype_manual(name = "Estimator:", labels = c("Naive", "DIV", "CIV"), values = c("dotdash", "88", "f1")) +
                              ggtitle(title))
  
  list_plots[[counter+3]] <- (ggplot(my_sum %>%
                                       filter((method %in% meths) & (TT == t)),
                                     aes(x = n,
                                         y = exp(est_sd),
                                         color = method, 
                                         linetype = method)) + 
                                geom_line() + 
                                scale_x_continuous(n.breaks = 5) + 
                                scale_y_continuous(n.breaks = 5) +
                                ylab("Standard Deviation") + 
                                xlab(TeX(r'($n_s$)')) +
                                scale_color_discrete(name = "Estimator:", labels = c("Naive", "DIV", "CIV")) + 
                                scale_linetype_manual(name = "Estimator:", labels = c("Naive", "DIV", "CIV"), values = c("dotdash", "88", "f1")))
  counter <- counter + 1
}
library(ggpubr)
for(counter in c(2,3,5,6))
  list_plots[[counter]] <- list_plots[[counter]] + rremove("ylab")
for(counter in 1:3)
  list_plots[[counter]] <- list_plots[[counter]] + rremove("xlab")


ggarrange(plotlist = list_plots, nrow = 2, ncol = 3, common.legend = TRUE)


