##############################
# Progarm who makes IV estimates in the NOR-DMARD trial, adjusting for time
############################



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

adeff <- readr::read_rds("data/adeff.rds")

source("Functions.R")
source("stata.R")

pz_margins <- function(data = adeff) {
  x <- glue::glue(
    "
tempfile tmp

mlogit  trtgrp i.Z_value
margins  i.Z_value, saving(`tmp')
use `tmp', clear

"
  )
  
  margins <-  stata(
    src = x,
    data.in = data,
    data.out = TRUE,
    stata.path = "/usr/local/bin/stata-se",
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
    select(trtgrp , estimate = margin, z_value = m1) %>%
    pivot_wider(names_from = trtgrp, values_from = estimate) %>%
    mutate(z_value = as.numeric(z_value))


  return(res)
  
}

pz_covar <- pz_margins(adeff) %>% 
  relocate(Z_value = z_value, Ada, Cer, Eta, Gol, `Inf`)
pz_covar


ptzm_margins <- function(data = adeff) {
  x <- glue::glue(
    "
tempfile tmp

constraint 1 cal_month = 0.01249
logistic das28crprem i.Z_value##i.trtgrp c.cal_month, constraint(1)
quietly margins Z_value#trtgrp,  saving(`tmp')
use `tmp', clear

"
  )
  
  margins <-  stata(
    src = x,
    data.in = data,
    data.out = TRUE,
    stata.path = "/usr/local/bin/stata-se",
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


ptzm_covar <- ptzm_margins(adeff) %>% 
  relocate(Z_value, Ada, Cer, Eta, Gol, `Inf`)
ptzm_covar

qz_covar <- pz_covar * ptzm_covar
qz_covar[is.na(qz_covar)] <- 0
qz_covar

plot(1:10, rowSums(qz_covar, na.rm = TRUE)-qz_covar$Z_value)

# 
# #Calculate the observables:
# Z_weights <- data %>% group_by(Z_value) %>% summarise(N = n())
# Z_weights$Prob <- Z_weights$N/sum(Z_weights$N)
# P_Z <- MakeP_Z(data, "Z_value", "trtgrp.x")
# Q_Z <- MakeQ_Z(data, "Z_value", "trtgrp.x", "das28crprem")

#Find the identifiables:
nZ <- length(Z)
names(Z) <- 1:nZ
A <- GenerateA(TLevels)
R <- MakeR(A, Z, T_decider)
KB <- MakeKB(R, TLevels, 4)
b <- KbSolver(KB, 3)
Pis <- PiIdentifier(b)

#P_Sigma <- P_SigmaIdentifier(P_Z, KB, b)
#LATEs <- LATEIdentifier(Q_Z, KB, b, P_Sigma)

P_Sigma_covar <- P_SigmaIdentifier(pz_covar, KB, b)
LATEs_covar <- LATEIdentifier(qz_covar, KB, b, P_Sigma_covar)

readr::write_rds(LATEs_covar, "data/LATEs_covar.rds")
readr::write_rds(qz_covar, "data/qz_covar.rds")

