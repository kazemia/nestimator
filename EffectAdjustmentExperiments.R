library(doParallel)
library(dplyr)
library(ggplot2)

IV_num <- function(data)
  (mean(data$Y[data$Z == 1]) - mean(data$Y[data$Z == 0]))
IV_den <- function(data)
  (mean(data$T[data$Z == 1]) - mean(data$T[data$Z == 0]))

V1T <- 0.1
V2T <- 0.6
ZT <- 0.3

V1Y <- 0.6
V2Y <- 0.1
TY <- 0.3

sim_list <- list()
sim_max <- 5000
sim_start <- 250
sim_step <- 250
n_sim <- 1000

n_cores <- detectCores()

myCluster <- makeCluster(n_cores-1)
registerDoParallel(myCluster)
for(nP in floor(sim_start/sim_step):floor(sim_max/sim_step)){
  sim_list[[nP]] <-foreach(counter=idiv(n_sim, chunks=getDoParWorkers()),
                           .combine = rbind) %dopar% {
    n <- nP*sim_step
    V1 <- rbinom(n,1,0.3)
    V2 <- rbinom(n,1,0.3)
    Z <- rbinom(n,1,0.5)
    T <- rbinom(rep(n,n),1,(V1T*V1) + (V2T*V2) + (ZT*Z))
    Y <- rbinom(rep(n,n),1,(V1Y*V1) + (V2Y*V2) + (TY*T))
    d <- data.frame(Z = Z, T = T, V1 = V1, V2 = V2, Y = Y)
    P_V1 <- mean(V1)
    P_V2 <- mean(V2)
    
    crude_num <- IV_num(d)
    crude_den <- IV_den(d)
    crude_IV <- crude_num/crude_den
    
    adj_num <- (P_V1*IV_num(d[d$V1 == 1,])) + ((1-P_V1)*IV_num(d[d$V1 == 0,]))
    adj_den <- (P_V1*IV_den(d[d$V1 == 1,])) + ((1-P_V1)*IV_den(d[d$V1 == 0,]))
    V1adj_IV <- (P_V1*IV_num(d[d$V1 == 1,])/IV_den(d[d$V1 == 1,])) + 
      ((1-P_V1)*IV_num(d[d$V1 == 0,])/IV_den(d[d$V1 == 0,]))
    num_V1adj_IV <- adj_num/crude_den
    den_V1adj_IV <- crude_num/adj_den
    
    adj_num <- (P_V2*IV_num(d[d$V2 == 1,])) + ((1-P_V2)*IV_num(d[d$V2 == 0,]))
    adj_den <- (P_V2*IV_den(d[d$V2 == 1,])) + ((1-P_V2)*IV_den(d[d$V2 == 0,]))
    V2adj_IV <- (P_V2*IV_num(d[d$V2 == 1,])/IV_den(d[d$V2 == 1,])) + 
      ((1-P_V2)*IV_num(d[d$V2 == 0,])/IV_den(d[d$V2 == 0,]))
    num_V2adj_IV <- adj_num/crude_den
    den_V2adj_IV <- crude_num/adj_den
    
    data.frame(n = rep(n , 7), 
               method = c("crude_IV", "V1adj_IV", "V2adj_IV",
                          "num_V1adj_IV", "den_V1adj_IV", 
                          "num_V2adj_IV", "den_V2adj_IV"),
               Estimate = c(crude_IV, V1adj_IV, V2adj_IV,
                            num_V1adj_IV, den_V1adj_IV, 
                            num_V2adj_IV, den_V2adj_IV))
    
                           }
}

sim_results <- sim_list[[1]]
for (counter in 2:length(sim_list)) {
  sim_results <- rbind(sim_results, sim_list[[counter]])
}

my_sum <- sim_results %>%
  group_by(n, method) %>% summarise(est_bias = abs(mean(Estimate, na.rm = TRUE) - 0.3) ,
                                    est_sd = sd(Estimate, na.rm = TRUE))

ggplot(my_sum %>% filter(n>999), aes(x = n, y = est_bias, color = method)) +
  geom_line() + scale_x_continuous(n.breaks = 15) %>% 
  print()

ggplot(my_sum %>% filter(n>999), aes(x = n, y = est_sd, color = method)) +
  geom_line() + scale_x_continuous(n.breaks = 15) %>% 
  print()

my_sum2 <- my_sum %>% as.data.frame() %>% filter(n>999) %>% group_by(method) %>% 
  summarise(mean_bias = mean(est_bias),
            mean_sd = mean(est_sd))
