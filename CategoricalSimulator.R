
TLevels <- c("Inf", "Ada", "Cer", "Eta", "Gol")
all_instruments <- permn(TLevels)
Z <- list("2010" = c("Inf", "Gol", "Cer", "Eta", "Ada"),
          "2011" = c("Eta", "Inf", "Cer", "Gol", "Ada"),
          "2012" = c("Inf", "Eta", "Cer", "Gol", "Ada"),
          "2013" = c("Cer", "Inf", "Gol", "Eta", "Ada"),
          "2014" = c("Inf", "Cer", "Gol", "Ada", "Eta"),
          "2015" = c("Inf", "Cer", "Gol", "Eta", "Ada"),
          "2016" = c("Inf", "Cer", "Eta", "Gol", "Ada"),
          "2017" = c("Inf", "Eta", "Gol", "Cer", "Ada"),
          "2018" = c("Inf", "Eta", "Cer", "Ada", "Gol"),
          "2019" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2020" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2021" = c("Ada", "Inf", "Eta", "Cer", "Gol"),
          "2022" = c("Ada", "Inf", "Eta", "Cer", "Gol")    
)
nZ <- length(Z)
names(Z) <- 1:nZ



A <- GenerateA(TLevels)
#A <- A[c(2,5,7)]
#names(A) <- paste0("A", 1:length(A))

R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

VLevels <- c(0,1)
VP <- c(0.75,0.25)
#VA1 <- sample.int(length(A), replace = T)
#VA2 <- sample.int(length(A), replace = T)
VA2prT <- list("Inf" = 5,"Ada" = 4,"Cer" = 3,"Eta" = 2, "Gol" = 1)
VA2 <- unlist(lapply(A, function(x){
  out <- 0
  for(t in x){
    out <- out + VA2prT[[t]]
  }
  return(out/length(x))
}))
VA1 <- max(unlist(VA2prT)) + 1 - VA2


VA <- list(VA1/sum(VA1),
           VA2/sum(VA2))
VY <- c(0,0.25)
TY <- c(0.1,0.25,0.4,0.55,0.7)
names(VA) <- VLevels
sim_max <- 10000
sim_step <- 500
n_sim <- 1000

myCluster <- makeCluster(14)
registerDoParallel(myCluster)
sim_res <- list()
for(nP in 1:floor(sim_max/sim_step)){
  sim_res[[nP]] <- foreach(counter=idiv(n_sim, chunks=getDoParWorkers()), .combine = rbind,
                     .packages = c("dplyr", "tidyverse", "stringr","lpSolve")
                     ) %dopar% {
                       out <- data.frame(TT = NA, n = NA, method = NA, estimate = NA)
                       for(counter2 in 1:counter){
                         sim_data <- CatSimulator(n = nP*sim_step,
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
                         temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = FALSE)
                         temp$Cer_Gol <- temp$Cer_Gol[2]
                         out <- rbind(out, data.frame(
                           TT = names(temp)[sapply(temp, length)>0],
                           n = rep(nP*sim_step, sum(sapply(temp, length))),
                           method = rep("LATE_IV", sum(sapply(temp, length))),
                           estimate = unlist(temp, use.names = FALSE)
                         ))
                         temp <- LATEIdentifier(Q_Z, KB, b, P_Sigma, RR = FALSE, AverageProb = TRUE)
                         temp$Cer_Gol <- temp$Cer_Gol[2]
                         out <- rbind(out, data.frame(
                           TT = names(temp)[sapply(temp, length)>0],
                           n = rep(nP*sim_step, sum(sapply(temp, length))),
                           method = rep("LATE_IV_wa", sum(sapply(temp, length))),
                           estimate = unlist(temp, use.names = FALSE)
                         ))
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

sim_results <- sim_res[[1]]
for (counter in 2:length(sim_res)) {
  sim_results <- rbind(sim_results, sim_res[[counter]])
}
remove(sim_res)

my_sum <- sim_results %>% group_by(TT, n, method) %>% summarise(est_mean = mean(estimate),
                                                                est_sd = sd(estimate))

library(ggplot2)
for(t in unique(my_sum$TT[my_sum$method == "LATE_IV"])){
  (ggplot(my_sum %>% filter((method %in% c("LATE_IV", "ATE_naive", "IV_naive")) & 
                             (TT == t)), aes(x = n, y = est_mean, color = method)) +
    geom_line() + scale_x_continuous(n.breaks = 15) +
     geom_hline(yintercept = TY[TLevels == str_split(t, "_")[[1]][1]] - 
                  TY[TLevels == str_split(t, "_")[[1]][2]]) + ggtitle(t)) %>% 
    print()
  
  (ggplot(my_sum %>% filter((method %in% c("LATE_IV", "ATE_naive", "IV_naive")) & 
                             (TT == t)), aes(x = n, y = est_sd, color = method)) +
    geom_line() + scale_x_continuous(n.breaks = 15)+ ggtitle(t)) %>% print()
}


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

