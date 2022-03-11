################################################################################################
#                                       Simulation study                                       #
#                                           Data                                               #
################################################################################################
#Last file: NA
#Author: Alma Revers
#Libraries needed
library(MASS)
library(Rlab)
library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(broom)


tree_pt_real <- readRDS(file = "example_tree.rds")

set.seed(3483)

data_sim <- list()
input_sim <- list()
outcome_sim <- list()

N_simulations <- 500

sim_N_patients <- 500

l_RR_mean = c(0, 1)
sd_input <- sqrt(0.1)
sim_N_PT <- 1000
filter_value = 4
filter_tree = 2
controle_incidence <- log(0.01)

for(s in 1:N_simulations){
  tree_int <- sample.int(size = sim_N_PT, n = dim(tree_pt_real)[2], replace=TRUE)
  sim_tree <- tree_pt_real[,tree_int]
  
  sim_N_SOC <- length(unique(sim_tree[1,]))
  sim_N_HLGT <- length(unique(sim_tree[2,]))
  sim_N_HLT <- length(unique(sim_tree[3,]))
  
  PT_HLT <- as.numeric(as.factor(sim_tree[3,]))
  sim_tree_hlt <- sim_tree[1:3,]
  sim_tree_hlt <- sim_tree_hlt[,!duplicated(sim_tree_hlt[3,])]
  HLT_HLGT <- as.numeric(as.factor(sim_tree_hlt[2,]))
  sim_tree_hlgt <- sim_tree_hlt[1:2,]
  sim_tree_hlgt <- sim_tree_hlgt[,!duplicated(sim_tree_hlgt[2,])]
  HLGT_SOC <- as.numeric(as.factor(sim_tree_hlgt[1,]))
  
  x <- rbern(n = sim_N_patients, prob = 0.5)
  
  mu_HLGT_t <- rep(0,  sim_N_HLGT)
  mu_HLGT_c <- rep(0,  sim_N_HLGT)
  for(i in 1:sim_N_SOC){
    mu_HLGT_c[HLGT_SOC==i] <- rnorm(n = sum(HLGT_SOC==i), mean = controle_incidence, sd=sd_input)
    mu_HLGT_t[HLGT_SOC==i] <- rnorm(n = sum(HLGT_SOC==i), mean = sample(size=1, x=l_RR_mean), sd=sd_input)
  }
  
  mu_HLT_c <- rep(0,  sim_N_HLT)
  mu_HLT_t <- rep(0,  sim_N_HLT)
  for(i in 1:sim_N_HLGT){
    mu_HLT_c[HLT_HLGT==i] <- rnorm(n = sum(HLT_HLGT==i), mean = mu_HLGT_c[i], sd=sd_input)
    mu_HLT_t[HLT_HLGT==i] <- rnorm(n = sum(HLT_HLGT==i), mean = mu_HLGT_t[i], sd=sd_input)
  }
  mu_PT_c <- rep(0,  sim_N_PT)
  mu_PT_t <- rep(0,  sim_N_PT)
  for(i in 1:sim_N_HLT){
    mu_PT_c[PT_HLT==i] <- rnorm(n = sum(PT_HLT==i), mean = mu_HLT_c[i], sd=sd_input)
    mu_PT_t[PT_HLT==i] <- rnorm(n = sum(PT_HLT==i), mean = mu_HLT_t[i], sd=sd_input)
  }
  
  y <- matrix(0, nrow = sim_N_patients, ncol = sim_N_PT)
  for(i in 1:sim_N_PT){
    y[,i] <- rpois(n = sim_N_patients, lambda = exp(mu_PT_c[i] + mu_PT_t[i]*x))
  }
  
  #y <- matrix(0, nrow = sim_N_patients, ncol = sim_N_PT)
  #totalnumber <- rnorm(n = sim_N_patients, mean = 50, sd = 10)
  #logmu=(cbind(rep(1,sim_N_patients),x) %*% rbind(mu_PT_c, mu_PT_t))
  #probs=exp(logmu)/sum(exp(logmu))
  #for(i in 1:sim_N_patients){
  #  y[i,] = rmultinom(n =1, size = totalnumber[i],prob = probs[i,])
  #}
  
  #Filtering
  y_pt_trt <- apply(y * x, 2, sum)
  y_pt_placebo <- apply(y * (as.numeric(factor(x, levels = c("1", "0")))-1), 2, sum)
  y_pt_def <- c()
  for(i in 1:length(y_pt_trt)){
    y_pt_def <- c(y_pt_def, (y_pt_trt[i] + y_pt_placebo[i])>filter_value)
  }
  t_hlt_pt <- table(as.factor(sim_tree[3,y_pt_def]))
  pt_to_few_in_HLT <- names(t_hlt_pt)[which(t_hlt_pt<filter_tree)]
  y_pt_def <- y_pt_def & !(sim_tree[3,] %in% pt_to_few_in_HLT)
  tree_pt_noPT <- sim_tree[,!y_pt_def]
  tree_pt <- sim_tree[,y_pt_def]
  
  y_pt_trt <- y_pt_trt[y_pt_def]
  y_pt_placebo <- y_pt_placebo[y_pt_def]
  N_PT_perHLT <- table(tree_pt_noPT[3,])
  N_PT_perHLGT <- table(tree_pt_noPT[2,])
  
  y_hlt_trt <- rep(0, length(unique(tree_pt_noPT[3,])))
  y_hlt_c <- rep(0, length(unique(tree_pt_noPT[3,])))
  for ( j in 1:length(unique(tree_pt_noPT[3,]))) {
    y_hlt_trt[j] = sum(y[x==1,sim_tree[3,] == (unique(tree_pt_noPT[3,])[j])])
    y_hlt_c[j] = sum(y[x==0,sim_tree[3,] == (unique(tree_pt_noPT[3,])[j])])
  }
  y_hlt_def <- c()
  for(i in 1:length(y_hlt_trt)){
    y_hlt_def <- c(y_hlt_def, (y_hlt_trt[i] + y_hlt_c[i])>filter_value)
  }
  t_hlgt_hlt <- table(as.factor(sim_tree_hlt[2,y_hlt_def]))
  hlt_to_few_in_HLGT <- names(t_hlgt_hlt)[which(t_hlgt_hlt<filter_tree)]
  y_hlt_def <- y_hlt_def & !(sim_tree_hlt[2,(sim_tree_hlt[3,] %in% tree_pt_noPT[3,])] %in% hlt_to_few_in_HLGT)
  tree_hlt_noHLT <- tree_pt_noPT[,!duplicated(tree_pt_noPT[3,])][,!y_hlt_def]
  tree_hlt <- tree_pt_noPT[,!duplicated(tree_pt_noPT[3,])][,y_hlt_def]
  N_PT_perHLT <- N_PT_perHLT[order(match(names(N_PT_perHLT), unique(tree_pt_noPT[3,])))][y_hlt_def]
  y_hlt_trt <- y_hlt_trt[y_hlt_def]
  y_hlt_c <- y_hlt_c[y_hlt_def]
  
  y_hlgt_trt <- rep(0, length(unique(tree_hlt_noHLT[2,])))
  y_hlgt_c <- rep(0, length(unique(tree_hlt_noHLT[2,])))
  for ( j in 1:length(unique(tree_hlt_noHLT[2,]))) {
    y_hlgt_trt[j] = sum(y[x==1,sim_tree[2,] ==(unique(tree_hlt_noHLT[2,])[j])])
    y_hlgt_c[j] = sum(y[x==0,sim_tree[2,] ==(unique(tree_hlt_noHLT[2,])[j])])
  }
  tree_hlgt <- sim_tree_hlgt[1:2, (sim_tree_hlgt[2,] %in% unique(tree_hlt_noHLT[2,]))]
  N_PT_perHLGT <- N_PT_perHLGT[names(N_PT_perHLGT) %in% unique(tree_hlt_noHLT[2,])]
  N_PT_perHLGT <- N_PT_perHLGT[order(match(names(N_PT_perHLGT), unique(tree_hlt_noHLT[2,])))]
  
  SOC_HLGT <- factor(c(tree_pt[1,!duplicated(tree_pt[2,])], tree_hlt[1,!duplicated(tree_hlt[2,])], tree_hlgt[1,]), 
                     levels = unique(c(tree_pt[1,], tree_hlt[1,], tree_hlgt[1,])))
  
  HLGT_HLT <- factor(c(tree_pt[2,!duplicated(tree_pt[3,])], tree_hlt[2,]),
                     levels = unique(c(tree_pt[2,], tree_hlt[2,])))
  
  data_stan_pt <- list(
    N_subj_trt = sum(x == 1),
    N_subj_cont = sum(x == 0),
    N_ae = length(y_pt_trt) + length(y_hlt_trt) + length(y_hlgt_trt),
    N_pt = length(y_pt_trt),
    N_hlt = length(unique(tree_pt[3,])) + length(y_hlt_trt),
    N_hlt_noPT = length(y_hlt_trt),
    N_hlt_PT = length(unique(tree_pt[3,])),
    N_hlgt = length(SOC_HLGT),
    N_hlgt_noHLT = length(y_hlgt_trt),
    N_hlgt_hlt = length(c(tree_pt[1,!duplicated(tree_pt[2,])], tree_hlt[1,!duplicated(tree_hlt[2,])])),
    N_soc = length(levels(SOC_HLGT)),
    y_hlt_trt = y_hlt_trt,
    y_hlt_placebo = y_hlt_c,
    y_hlgt_trt= y_hlgt_trt,
    y_hlgt_placebo= y_hlgt_c,
    y_pt_trt= y_pt_trt,
    y_pt_placebo= y_pt_placebo,
    SOC_HLGT = as.numeric(SOC_HLGT),
    HLGT_HLT= as.numeric(HLGT_HLT),
    HLT_PT= as.numeric(as.factor(tree_pt[3,])),
    N_PT_perHLT = N_PT_perHLT,
    N_PT_perHLGT = N_PT_perHLGT
  )
  data_sim[[s]] <- data_stan_pt
  
  input_data_list <- list(
    sim_N_patients = sim_N_patients,
    sim_N_PT = sim_N_PT,
    controle_incidence = controle_incidence,
    l_RR_mean = l_RR_mean,
    sim_tree = sim_tree,
    sim_N_SOC = sim_N_SOC,
    sim_N_HLGT = sim_N_HLGT,
    sim_N_HLT = sim_N_HLT,
    x = x,
    mu_HLGT_t = mu_HLGT_t,
    mu_HLGT_c = mu_HLGT_c,
    mu_HLT_c = mu_HLT_c,
    mu_HLT_t = mu_HLT_t,
    mu_PT_c = mu_PT_c,
    mu_PT_t = mu_PT_t,
    y = y,
    y_pt_def = y_pt_def,
    y_hlt_def = y_hlt_def
  )
  input_sim[[s]] <- input_data_list
  
  model_obj <- stan_model(file="AE_simulations.stan")
  stan_pt_poisson <- sampling(model_obj, data=data_stan_pt, chains=4)
  outcome_sim[[s]] <- stan_pt_poisson
  print(s)
}

save.image(file="sim_1000_500_0_1.RData")






