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
VT = as.matrix(data.frame(V1T = c(0,0,0,exp(3),0), V2T = c(0,0,exp(3),0,0)))
VP = c(0.3,0.3)
VY = c(0.25, 0.25)
TY = c(0.05,0.15,0.25,0.35,0.45)
Vadj = c(TRUE, FALSE)

#Specify maximum number of observations in simulation, the step 
#and number of simulations pr step
sim_max <- 10000
sim_start <- 750
sim_step <- 250
n_sim <- 5000
TT_vec <- c("Inf_Cer", "Inf_Eta", "Ada_Eta", "Ada_Gol", "Cer_Eta", "Cer_Gol1", "Cer_Gol2", "Cer_Gol3", "Eta_Gol")

#Simulate
myCluster <- makeCluster(15)
registerDoParallel(myCluster)
sim_res <- list()
for(nP in floor(sim_start/sim_step):floor(sim_max/sim_step)){
  sim_res[[nP]] <- foreach(counter=idiv(n_sim, chunks=getDoParWorkers()), .combine = rbind,
                     .packages = c("dplyr", "tidyverse", "stringr","lpSolve")
                     ) %dopar% {
                       out <- data.frame(TT = NA, n = NA, method = NA, estimate = NA)
                       for(counter2 in 1:counter){
                         sim_data <- CatSimulator2(nP*sim_step, Z, TLevels, VT,
                                                   VP, VY, TY, Vadj)
                         P_Z <- MakeP_Z(sim_data, "Z", "T")
                         Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
                         P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
                         temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = FALSE)
                         out <- rbind(out, data.frame(
                           TT = TT_vec,
                           n = rep(nP*sim_step, length(TT_vec)),
                           method = rep("LATE_IV_unadj", length(TT_vec)),
                           estimate = unlist(temp, use.names = FALSE)
                         ))
                         temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = TRUE)
                         out <- rbind(out, data.frame(
                           TT = TT_vec,
                           n = rep(nP*sim_step, length(TT_vec)),
                           method = rep("LATE_IV_wa_unadj", length(TT_vec)),
                           estimate = unlist(temp, use.names = FALSE)
                         ))
                         for(tt in names(temp)[sapply(temp, length)==1]){
                           out <- rbind(out, data.frame(
                             TT = rep(tt,2),
                             n = rep(nP*sim_step, 2),
                             method = c("ATE_naive", "IV_naive"),
                             estimate = c(mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][1]]) - 
                                            mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][2]]),
                                          NaiveIV(tt, sim_data, Z, "T", "Z", "Y"))
                           ))
                         }
                         tt <- "Cer_Gol"
                         out <- rbind(out, data.frame(
                           TT = c("Cer_Gol1", "Cer_Gol2", "Cer_Gol3"),
                           n = rep(nP*sim_step, 3),
                           method = rep("ATE_naive", 3),
                           estimate = rep(mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][1]]) - 
                                          mean(sim_data$Y[sim_data$T == str_split(tt, "_")[[1]][2]]), 3)
                         ))
                         out <- rbind(out, data.frame(
                           TT = c("Cer_Gol1", "Cer_Gol2", "Cer_Gol3"),
                           n = rep(nP*sim_step, 3),
                           method = rep("IV_naive", 3),
                           estimate = rep(NaiveIV(tt, sim_data, Z, "T", "Z", "Y"), 3)
                         ))
                         P_v1 <- mean(sim_data$V1)
                         sim_data_back <- sim_data
                         sim_data <- sim_data_back %>% filter(V1 == 0)
                         P_Z <- MakeP_Z(sim_data, "Z", "T")
                         Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
                         P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
                         temp0 <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE,
                                                 AverageProb = FALSE) %>% 
                           unlist(use.names = FALSE)
                         sim_data <- sim_data_back %>% filter(V1 == 1)
                         P_Z <- MakeP_Z(sim_data, "Z", "T")
                         Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
                         P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
                         temp1 <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE,
                                                 AverageProb = FALSE) %>% 
                           unlist(use.names = FALSE)
                         
                         out <- rbind(out, data.frame(
                           TT = TT_vec,
                           n = rep(nP*sim_step, length(TT_vec)),
                           method = rep("LATE_IV_adj", length(TT_vec)),
                           estimate = (temp1*P_v1) + (temp0*(1-P_v1))
                         ))
                         sim_data <- sim_data_back %>% filter(V1 == 0)
                         P_Z <- MakeP_Z(sim_data, "Z", "T")
                         Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
                         P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
                         temp0 <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE,
                                                 AverageProb = TRUE) %>% 
                           unlist(use.names = FALSE)
                         sim_data <- sim_data_back %>% filter(V1 == 1)
                         P_Z <- MakeP_Z(sim_data, "Z", "T")
                         Q_Z <- MakeQ_Z(sim_data, "Z", "T", "Y")
                         P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
                         temp1 <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE,
                                                 AverageProb = TRUE) %>% 
                           unlist(use.names = FALSE)
                         
                         out <- rbind(out, data.frame(
                           TT = TT_vec,
                           n = rep(nP*sim_step, length(TT_vec)),
                           method = rep("LATE_IV_wa_adj", length(TT_vec)),
                           estimate = (temp1*P_v1) + (temp0*(1-P_v1))
                         ))
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

sim_results$estimate[sim_results$estimate < (-1)] <- (-1)
sim_results$estimate[sim_results$estimate > 1] <- 1

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



my_sum <- my_sum %>% filter(n>999)
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
  if(nchar(t) == 8)
    title <- TeX(paste0(r'($\hat{E})', r'(_)', substring(t,8), r'((Y(t_)', t1, r'()-Y(t_)', t2, r'())$)' ))
  list_plots[[counter]] <- (ggplot(my_sum %>%
                               filter((method %in% meths) & (TT == t)),
                             aes(x = n,
                                 y = est_mean,
                                 color = method)) + 
                        geom_line() + 
                        geom_hline(yintercept = (TY[t1] - TY[t2])) +
                        scale_x_continuous(n.breaks = 5) + 
                        scale_y_continuous(n.breaks = 5) +
                        ylab("Median") + 
                        xlab(TeX(r'($n_s$)')) +
                        scale_color_discrete(name = "Estimator:", labels = c("Naive", "DIV", "CIV")) + 
                        ggtitle(title))
  
  list_plots[[counter+9]] <- (ggplot(my_sum %>%
                               filter((method %in% meths) & (TT == t)),
                             aes(x = n,
                                 y = est_sd,
                                 color = method)) + 
                        geom_line() + 
                        scale_x_continuous(n.breaks = 5) + 
                        scale_y_continuous(n.breaks = 5) +
                        ylab("Log standard deviation") + 
                        xlab(TeX(r'($n_s$)')) +
                        scale_color_discrete(name = "Estimator:", labels = c("Naive", "DIV", "CIV")))
  counter <- counter + 1
}
ggarrange(plotlist = list_plots, nrow = 3, ncol = 6, common.legend = TRUE)
for(counter in (1:18)[-c(1,10,4,13,7,16)])
  list_plots[[counter]] <- list_plots[[counter]] + rremove("ylab")
for(counter in 1:9)
  list_plots[[counter]] <- list_plots[[counter]] + rremove("xlab")

library(ggpubr)
ggarrange(plotlist = list_plots[c(1:3,10:12)], nrow = 2, ncol = 3, common.legend = TRUE)
ggarrange(plotlist = list_plots[c(4:6,13:15)], nrow = 2, ncol = 3, common.legend = TRUE)
ggarrange(plotlist = list_plots[c(7:9,16:18)], nrow = 2, ncol = 3, common.legend = TRUE)

grid.arrange(grobs = list_plots[1:4], ncol = 2)
grid.arrange(grobs = list_plots[5:10], ncol = 2)
grid.arrange(grobs = list_plots[11:16], ncol = 2)
grid.arrange()

for(t in TTs){
  (ggplot(my_sum %>%
           filter((method %in% meths) & (TT == t)),
         aes(x = n, y = est_bias, ymin=est_lq, ymax=est_hq,
             fill = method, linetype = method)) +
    geom_line() + geom_ribbon(alpha=0.5)+ ggtitle(t)) %>% print()
}

for(t in TTs){
  (ggplot(my_sum %>%
            filter((method %in% meths) & (TT == t)),
          aes(x = n, y = est_bias,
              fill = method, linetype = method)) +
     geom_line() + ggtitle(t)) %>% print()
}



#THE CODE AFTER THIS IS NOT IN USE


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

