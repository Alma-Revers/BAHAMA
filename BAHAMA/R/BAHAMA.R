#
#'
#' @param BAHAMADataSet A BAHAMADataSet object
BAHAMA <- function(BAHAMADataSet){
  if(BAHAMADataSet@N_hlgt_noHLT > 0){
    model_code = "data {
  int N_subj_trt; //Number of subjects in trt group
  int N_subj_cont; //Number of subjects in control group
  int N_ae;
  int N_pt;
  int N_hlt;
  int N_hlt_noPT;
  int N_hlt_PT;
  int N_hlgt; //Number of hlgt adverse events included
  int N_hlgt_noHLT;
  int N_hlgt_hlt;
  int N_soc; //Number of SOCs included
  int y_hlt_trt[N_hlt_noPT]; //Number of events in trt group
  int y_hlt_cont[N_hlt_noPT]; //number of events in controle group
  int y_hlgt_trt[N_hlgt_noHLT];
  int y_hlgt_cont[N_hlgt_noHLT];
  int y_pt_trt[N_pt];
  int y_pt_cont[N_pt];
  matrix[N_pt, N_hlt_PT] w_pt_hlt; //
  matrix[N_hlt, N_hlgt] w_hlt_hlgt;
  matrix[N_hlgt, N_soc] w_hlgt_soc;
  int N_PT_perHLT[N_hlt_noPT];
  int N_PT_perHLGT[N_hlgt_noHLT];
}
parameters {
  vector[2] mu; //Overall mean effect
  vector<lower=0>[2] sigma; //Overal variance

  vector[N_soc] mu_c_soc; //Mean effect per SOC
  vector<lower=0>[N_soc] sigma_c_soc; //Variance effect per SOC
  vector[N_soc] mu_t_soc; //Mean effect per SOC
  vector<lower=0>[N_soc] sigma_t_soc; //Variance effect per SOC

  vector[N_hlgt] mu_c_hlgt; //Mean effect per SOC
  vector<lower=0>[N_hlgt] sigma_c_hlgt; //Variance effect per SOC
  vector[N_hlgt] mu_t_hlgt; //Mean effect per SOC
  vector<lower=0>[N_hlgt] sigma_t_hlgt; //Variance effect per SOC

   vector[N_hlt] mu_c_hlt; //Mean effect per SOC
  vector<lower=0>[N_hlt] sigma_c_hlt; //Variance effect per SOC
  vector[N_hlt] mu_t_hlt; //Mean effect per SOC
  vector<lower=0>[N_hlt] sigma_t_hlt; //Variance effect per SOC

  vector[N_pt] lambda_pt;
  vector[N_pt] theta_pt;
}
transformed parameters {
  vector[N_hlgt_noHLT] tau_hlgt;
  vector[N_hlt_noPT] tau_hlt;
  vector[N_pt] tau_pt;

  vector[N_pt] mu_c_hlt1;
  vector[N_pt] sigma_c_hlt1;
  vector[N_pt] mu_t_hlt1;
  vector[N_pt] sigma_t_hlt1;

  vector[N_hlt] mu_c_hlgt1;
  vector[N_hlt] sigma_c_hlgt1;
  vector[N_hlt] mu_t_hlgt1;
  vector[N_hlt] sigma_t_hlgt1;

  vector[N_hlgt] mu_c_soc1;
  vector[N_hlgt] sigma_c_soc1;
  vector[N_hlgt] mu_t_soc1;
  vector[N_hlgt] sigma_t_soc1;

  mu_c_hlt1 = w_pt_hlt * mu_c_hlt[1:N_hlt_PT];
  sigma_c_hlt1 = w_pt_hlt * sigma_c_hlt[1:N_hlt_PT];
  mu_t_hlt1 = w_pt_hlt * mu_t_hlt[1:N_hlt_PT];
  sigma_t_hlt1 = w_pt_hlt * sigma_t_hlt[1:N_hlt_PT];

  mu_c_hlgt1 = w_hlt_hlgt * mu_c_hlgt;
  sigma_c_hlgt1 =  w_hlt_hlgt * sigma_c_hlgt;
  mu_t_hlgt1 = w_hlt_hlgt * mu_t_hlgt;
  sigma_t_hlgt1 =  w_hlt_hlgt * sigma_t_hlgt;

  mu_c_soc1 = w_hlgt_soc * mu_c_soc;
  sigma_c_soc1 =  w_hlgt_soc * sigma_c_soc;
  mu_t_soc1 = w_hlgt_soc * mu_t_soc;
  sigma_t_soc1 =  w_hlgt_soc * sigma_t_soc;

  for(j in 1:N_hlgt_noHLT){
    tau_hlgt[j] = mu_t_hlgt[(j+N_hlgt_hlt)] + mu_c_hlgt[(j+N_hlgt_hlt)];
  }
  for(i in 1:N_hlt_noPT){
    tau_hlt[i] = mu_t_hlt[(i+N_hlt_PT)] + mu_c_hlt[(i+N_hlt_PT)];
  }
  tau_pt = theta_pt + lambda_pt;
}
model {
  mu ~ std_normal();
  sigma ~ exponential(1);

  mu_c_soc ~ normal(mu[1],sigma[1]);
  mu_t_soc ~ normal(mu[2],sigma[2]);
  sigma_c_soc ~ exponential(1);
  sigma_t_soc ~ exponential(1);
  sigma_c_hlgt ~ exponential(1);
  sigma_t_hlgt ~ exponential(1);
  sigma_c_hlt ~ exponential(1);
  sigma_t_hlt ~ exponential(1);

  mu_c_hlgt ~  normal( mu_c_soc1, sigma_c_soc1 );
  mu_t_hlgt ~  normal( mu_t_soc1, sigma_t_soc1 );

  for(k in (N_hlgt_hlt+1):N_hlgt){
    y_hlgt_trt[(k-N_hlgt_hlt)] ~ poisson(exp(tau_hlgt[(k-N_hlgt_hlt)] + log(N_subj_trt)) * N_PT_perHLGT[(k-N_hlgt_hlt)]);
    y_hlgt_cont[(k-N_hlgt_hlt)] ~ poisson(exp(mu_c_hlgt[k] + log(N_subj_cont)) * N_PT_perHLGT[(k-N_hlgt_hlt)]);
  }

  mu_c_hlt ~ normal( mu_c_hlgt1, sigma_c_hlgt1 );
  mu_t_hlt ~ normal( mu_t_hlgt1, sigma_t_hlgt1  );
  for ( i in (N_hlt_PT+1):N_hlt){
    y_hlt_trt[(i-N_hlt_PT)] ~ poisson(exp(tau_hlt[(i-N_hlt_PT)] + log(N_subj_trt)) * N_PT_perHLT[(i-N_hlt_PT)]);
    y_hlt_cont[(i-N_hlt_PT)] ~ poisson(exp(mu_c_hlt[i] + log(N_subj_cont)) * N_PT_perHLT[(i-N_hlt_PT)]);
  }

  lambda_pt ~ normal(mu_c_hlt1, sigma_c_hlt1);
  theta_pt ~  normal(mu_t_hlt1, sigma_t_hlt1);
  y_pt_trt ~ poisson_log(tau_pt + log(N_subj_trt));
  y_pt_cont ~ poisson_log(lambda_pt + log(N_subj_cont));
 }"

    data_list = list(N_subj_trt = BAHAMADataSet@N_trt,
                 N_subj_cont = BAHAMADataSet@N_cont,
                 N_ae = BAHAMADataSet@N_ae,
                 N_pt = BAHAMADataSet@N_pt,
                 N_hlt = BAHAMADataSet@N_hlt,
                 N_hlt_noPT = BAHAMADataSet@N_hlt_noPT,
                 N_hlt_PT = BAHAMADataSet@N_hlt_PT,
                 N_hlgt = BAHAMADataSet@N_hlgt,
                 N_hlgt_noHLT = BAHAMADataSet@N_hlgt_noHLT,
                 N_hlgt_hlt = BAHAMADataSet@N_hlgt_hlt,
                 N_soc = BAHAMADataSet@N_soc,
                 y_hlt_trt = BAHAMADataSet@y_hlt_trt,
                 y_hlt_cont = BAHAMADataSet@y_hlt_cont,
                 y_hlgt_trt = BAHAMADataSet@y_hlgt_trt,
                 y_hlgt_cont = BAHAMADataSet@y_hlgt_cont,
                 y_pt_trt = BAHAMADataSet@y_pt_trt,
                 y_pt_cont = BAHAMADataSet@y_pt_cont,
                 w_pt_hlt = BAHAMADataSet@w_pt_hlt,
                 w_hlt_hlgt = BAHAMADataSet@w_hlt_hlgt,
                 w_hlgt_soc = BAHAMADataSet@w_hlgt_soc,
                 N_PT_perHLT = BAHAMADataSet@N_PT_perHLT,
                 N_PT_perHLGT = BAHAMADataSet@N_PT_perHLGT
    )
  }
else{
  model_code = "data {
  int N_subj_trt; //Number of subjects in trt group
  int N_subj_cont; //Number of subjects in control group
  int N_ae;
  int N_pt;
  int N_hlt;
  int N_hlt_noPT;
  int N_hlt_PT;
  int N_hlgt; //Number of hlgt adverse events included
  int N_soc; //Number of SOCs included
  int y_hlt_trt[N_hlt_noPT]; //Number of events in trt group
  int y_hlt_cont[N_hlt_noPT]; //number of events in controle group
  int y_pt_trt[N_pt];
  int y_pt_cont[N_pt];
  matrix[N_pt, N_hlt_PT] w_pt_hlt; //
  matrix[N_hlt, N_hlgt] w_hlt_hlgt;
  matrix[N_hlgt, N_soc] w_hlgt_soc;
  int N_PT_perHLT[N_hlt_noPT];
}
parameters {
  vector[2] mu; //Overall mean effect
  vector<lower=0>[2] sigma; //Overal variance

  vector[N_soc] mu_c_soc; //Mean effect per SOC
  vector<lower=0>[N_soc] sigma_c_soc; //Variance effect per SOC
  vector[N_soc] mu_t_soc; //Mean effect per SOC
  vector<lower=0>[N_soc] sigma_t_soc; //Variance effect per SOC

  vector[N_hlgt] mu_c_hlgt; //Mean effect per SOC
  vector<lower=0>[N_hlgt] sigma_c_hlgt; //Variance effect per SOC
  vector[N_hlgt] mu_t_hlgt; //Mean effect per SOC
  vector<lower=0>[N_hlgt] sigma_t_hlgt; //Variance effect per SOC

   vector[N_hlt] mu_c_hlt; //Mean effect per SOC
  vector<lower=0>[N_hlt] sigma_c_hlt; //Variance effect per SOC
  vector[N_hlt] mu_t_hlt; //Mean effect per SOC
  vector<lower=0>[N_hlt] sigma_t_hlt; //Variance effect per SOC

  vector[N_pt] lambda_pt;
  vector[N_pt] theta_pt;
}
transformed parameters {
  vector[N_hlt_noPT] tau_hlt;
  vector[N_pt] tau_pt;

  vector[N_pt] mu_c_hlt1;
  vector[N_pt] sigma_c_hlt1;
  vector[N_pt] mu_t_hlt1;
  vector[N_pt] sigma_t_hlt1;

  vector[N_hlt] mu_c_hlgt1;
  vector[N_hlt] sigma_c_hlgt1;
  vector[N_hlt] mu_t_hlgt1;
  vector[N_hlt] sigma_t_hlgt1;

  vector[N_hlgt] mu_c_soc1;
  vector[N_hlgt] sigma_c_soc1;
  vector[N_hlgt] mu_t_soc1;
  vector[N_hlgt] sigma_t_soc1;

  mu_c_hlt1 = w_pt_hlt * mu_c_hlt[1:N_hlt_PT];
  sigma_c_hlt1 = w_pt_hlt * sigma_c_hlt[1:N_hlt_PT];
  mu_t_hlt1 = w_pt_hlt * mu_t_hlt[1:N_hlt_PT];
  sigma_t_hlt1 = w_pt_hlt * sigma_t_hlt[1:N_hlt_PT];

  mu_c_hlgt1 = w_hlt_hlgt * mu_c_hlgt;
  sigma_c_hlgt1 =  w_hlt_hlgt * sigma_c_hlgt;
  mu_t_hlgt1 = w_hlt_hlgt * mu_t_hlgt;
  sigma_t_hlgt1 =  w_hlt_hlgt * sigma_t_hlgt;

  mu_c_soc1 = w_hlgt_soc * mu_c_soc;
  sigma_c_soc1 =  w_hlgt_soc * sigma_c_soc;
  mu_t_soc1 = w_hlgt_soc * mu_t_soc;
  sigma_t_soc1 =  w_hlgt_soc * sigma_t_soc;

  for(i in 1:N_hlt_noPT){
    tau_hlt[i] = mu_t_hlt[(i+N_hlt_PT)] + mu_c_hlt[(i+N_hlt_PT)];
  }
  tau_pt = theta_pt + lambda_pt;
}
model {
  mu ~ std_normal();
  sigma ~ exponential(1);

  mu_c_soc ~ normal(mu[1],sigma[1]);
  mu_t_soc ~ normal(mu[2],sigma[2]);
  sigma_c_soc ~ exponential(1);
  sigma_t_soc ~ exponential(1);
  sigma_c_hlgt ~ exponential(1);
  sigma_t_hlgt ~ exponential(1);
  sigma_c_hlt ~ exponential(1);
  sigma_t_hlt ~ exponential(1);

  mu_c_hlgt ~  normal( mu_c_soc1, sigma_c_soc1 );
  mu_t_hlgt ~  normal( mu_t_soc1, sigma_t_soc1 );

  mu_c_hlt ~ normal( mu_c_hlgt1, sigma_c_hlgt1 );
  mu_t_hlt ~ normal( mu_t_hlgt1, sigma_t_hlgt1  );
  for ( i in (N_hlt_PT+1):N_hlt){
    y_hlt_trt[(i-N_hlt_PT)] ~ poisson(exp(tau_hlt[(i-N_hlt_PT)] + log(N_subj_trt)) * N_PT_perHLT[(i-N_hlt_PT)]);
    y_hlt_cont[(i-N_hlt_PT)] ~ poisson(exp(mu_c_hlt[i] + log(N_subj_cont)) * N_PT_perHLT[(i-N_hlt_PT)]);
  }

  lambda_pt ~ normal(mu_c_hlt1, sigma_c_hlt1);
  theta_pt ~  normal(mu_t_hlt1, sigma_t_hlt1);
  y_pt_trt ~ poisson_log(tau_pt + log(N_subj_trt));
  y_pt_cont ~ poisson_log(lambda_pt + log(N_subj_cont));
 }"

data_list = list(N_subj_trt = BAHAMADataSet@N_trt,
                 N_subj_cont = BAHAMADataSet@N_cont,
                 N_ae = BAHAMADataSet@N_ae,
                 N_pt = BAHAMADataSet@N_pt,
                 N_hlt = BAHAMADataSet@N_hlt,
                 N_hlt_noPT = BAHAMADataSet@N_hlt_noPT,
                 N_hlt_PT = BAHAMADataSet@N_hlt_PT,
                 N_hlgt = BAHAMADataSet@N_hlgt,
                 N_soc = BAHAMADataSet@N_soc,
                 y_hlt_trt = BAHAMADataSet@y_hlt_trt,
                 y_hlt_cont = BAHAMADataSet@y_hlt_cont,
                 y_pt_trt = BAHAMADataSet@y_pt_trt,
                 y_pt_cont = BAHAMADataSet@y_pt_cont,
                 w_pt_hlt = BAHAMADataSet@w_pt_hlt,
                 w_hlt_hlgt = BAHAMADataSet@w_hlt_hlgt,
                 w_hlgt_soc = BAHAMADataSet@w_hlgt_soc,
                 N_PT_perHLT = BAHAMADataSet@N_PT_perHLT
)
}

options(mc.cores = parallel::detectCores())
model_obj <- rstan::stan_model(model_code = model_code)
stan_fit <- rstan::sampling(model_obj, data=data_list, chains=4)

tidy_rr_pt <- broom.mixed::tidy(stan_fit, pars="theta_pt", estimate.method = "mean", conf.int=TRUE, conf.level = 0.95,
                   rhat=TRUE,
                   ess = TRUE)
tidy_rr_hlt <- broom.mixed::tidy(stan_fit, pars="mu_t_hlt", estimate.method = "mean", conf.int=TRUE, conf.level = 0.95, rhat=TRUE,
                    ess = TRUE)
tidy_rr_hlgt <- broom.mixed::tidy(stan_fit, pars="mu_t_hlgt", estimate.method = "mean", conf.int=TRUE, conf.level = 0.95, rhat=TRUE,
                     ess = TRUE)
tidy_rr_soc <- broom.mixed::tidy(stan_fit, pars="mu_t_soc", estimate.method = "mean", conf.int=TRUE, conf.level = 0.95, rhat=TRUE,
                    ess = TRUE)

tidy_rr_pt$label <- names(BAHAMADataSet@y_pt_trt)
tidy_rr_hlt$label <- rownames(BAHAMADataSet@w_hlt_hlgt)
tidy_rr_hlgt$label <- rownames(BAHAMADataSet@w_hlgt_soc)
tidy_rr_soc$label <- colnames(BAHAMADataSet@w_hlgt_soc)

PT_abs_diff <- abs(BAHAMADataSet@y_pt_trt - BAHAMADataSet@y_pt_cont)
y_pt_c <- BAHAMADataSet@y_pt_cont
y_pt_t <- BAHAMADataSet@y_pt_trt

HLT_abs_diff <- c()
y_hlt_c <- c()
y_hlt_t <- c()
for(i in 1:BAHAMADataSet@N_hlt_PT){
  y_hlt_c <- c(y_hlt_c, sum(y_pt_c[BAHAMADataSet@w_pt_hlt[,i]>0]))
  y_hlt_t <- c(y_hlt_t, sum(y_pt_t[BAHAMADataSet@w_pt_hlt[,i]>0]))
  HLT_abs_diff <- c(HLT_abs_diff, sum(PT_abs_diff[BAHAMADataSet@w_pt_hlt[,i]>0]))
}
HLT_abs_diff <- c(HLT_abs_diff, abs(BAHAMADataSet@y_hlt_trt - BAHAMADataSet@y_hlt_cont) )
y_hlt_c <- c(y_hlt_c, BAHAMADataSet@y_hlt_cont)
y_hlt_t <- c(y_hlt_t, BAHAMADataSet@y_hlt_trt)

HLGT_abs_diff <- c()
y_hlgt_c <- c()
y_hlgt_t <- c()
for(i in 1:BAHAMADataSet@N_hlgt_hlt){
  y_hlgt_c <- c(y_hlgt_c, sum(y_hlt_c[BAHAMADataSet@w_hlt_hlgt[,i]>0]))
  y_hlgt_t <- c(y_hlgt_t, sum(y_hlt_t[BAHAMADataSet@w_hlt_hlgt[,i]>0]))
  HLGT_abs_diff <- c(HLGT_abs_diff, sum(HLT_abs_diff[BAHAMADataSet@w_hlt_hlgt[,i]>0]) )
}
HLGT_abs_diff <- c(HLGT_abs_diff, abs(BAHAMADataSet@y_hlgt_trt - BAHAMADataSet@y_hlgt_cont) )
y_hlgt_c <- c(y_hlgt_c, BAHAMADataSet@y_hlgt_cont)
y_hlgt_t <- c(y_hlgt_t, BAHAMADataSet@y_hlgt_trt)

SOC_abs_diff <- c()
y_soc_c <- c()
y_soc_t <- c()
for(i in 1:BAHAMADataSet@N_soc){
  y_soc_c <- c( y_soc_c, sum(y_hlgt_c[BAHAMADataSet@w_hlgt_soc[,i]>0]) )
  y_soc_t <- c( y_soc_t, sum(y_hlgt_t[BAHAMADataSet@w_hlgt_soc[,i]>0]) )
  SOC_abs_diff <- c( SOC_abs_diff, sum(HLGT_abs_diff[BAHAMADataSet@w_hlgt_soc[,i]>0]) )
}

tidy_rr_pt$abs_diff <- PT_abs_diff
tidy_rr_hlt$abs_diff <- HLT_abs_diff
tidy_rr_hlgt$abs_diff <- HLGT_abs_diff
tidy_rr_soc$abs_diff <- SOC_abs_diff

tidy_rr_pt$y_pt_c <- y_pt_c
tidy_rr_hlt$y_hlt_c <- y_hlt_c
tidy_rr_hlgt$y_hlgt_c <- y_hlgt_c
tidy_rr_soc$y_soc_c <- y_soc_c
tidy_rr_pt$y_pt_t <- y_pt_t
tidy_rr_hlt$y_hlt_t <- y_hlt_t
tidy_rr_hlgt$y_hlgt_t <- y_hlgt_t
tidy_rr_soc$y_soc_t <- y_soc_t

soc_hlgt_labels <- c()
for(i in 1:BAHAMADataSet@N_hlgt){
  if(sum(BAHAMADataSet@w_hlgt_soc[i,]>0) == 1){
    soc_hlgt_labels <- c(soc_hlgt_labels, colnames(BAHAMADataSet@w_hlgt_soc)[BAHAMADataSet@w_hlgt_soc[i,]>0])
  }
  else{
    soc_hlgt_labels <- c(soc_hlgt_labels, colnames(BAHAMADataSet@w_hlgt_soc)[which.max(BAHAMADataSet@w_hlgt_soc[i,])[1] ])
  }
}

soc_hlt_labels <- c()
hlgt_hlt_labels <- c()
for(i in 1:BAHAMADataSet@N_hlt){
  if(sum(BAHAMADataSet@w_hlt_hlgt[i,]>0) == 1){
    soc_hlt_labels <- c(soc_hlt_labels, soc_hlgt_labels[BAHAMADataSet@w_hlt_hlgt[i,] == 1])
    hlgt_hlt_labels <- c(hlgt_hlt_labels, colnames(BAHAMADataSet@w_hlt_hlgt)[BAHAMADataSet@w_hlt_hlgt[i,] == 1])
  }
  else{
    soc_hlt_labels <- c(soc_hlt_labels, soc_hlgt_labels[which.max(BAHAMADataSet@w_hlt_hlgt[i,])[1] ])
    hlgt_hlt_labels <- c(hlgt_hlt_labels, colnames(BAHAMADataSet@w_hlt_hlgt)[ which.max(BAHAMADataSet@w_hlt_hlgt[i,])[1] ])
  }
}

soc_pt_labels <- c()
hlgt_pt_labels <- c()
hlt_pt_labels <- c()
for(i in 1:BAHAMADataSet@N_pt){
  if(sum(BAHAMADataSet@w_pt_hlt[i,]>0) == 1){
    soc_pt_labels <- c(soc_pt_labels, soc_hlt_labels[BAHAMADataSet@w_pt_hlt[i,] == 1][1])
    hlgt_pt_labels <- c(hlgt_pt_labels, hlgt_hlt_labels[BAHAMADataSet@w_pt_hlt[i,] == 1][1])
    hlt_pt_labels <- c(hlt_pt_labels, colnames(BAHAMADataSet@w_pt_hlt)[BAHAMADataSet@w_pt_hlt[i,] == 1][1])
  }
  else{
    soc_pt_labels <- c(soc_pt_labels, soc_hlt_labels[which.max(BAHAMADataSet@w_pt_hlt[i,])[1] ])
    hlgt_pt_labels <- c(hlgt_pt_labels, hlgt_hlt_labels[ which.max(BAHAMADataSet@w_pt_hlt[i,])[1] ])
    hlt_pt_labels <- c(hlt_pt_labels, colnames(BAHAMADataSet@w_pt_hlt)[ which.max(BAHAMADataSet@w_pt_hlt[i,])[1] ])
  }
}
tidy_rr_soc$label <- as.factor(tidy_rr_soc$label)
tidy_rr_hlgt$soc <- factor(soc_hlgt_labels, levels= levels(tidy_rr_soc$label))
tidy_rr_hlt$soc <- factor(soc_hlt_labels, levels= levels(tidy_rr_soc$label))
tidy_rr_pt$soc <- factor(soc_pt_labels, levels= levels(tidy_rr_soc$label))
tidy_rr_hlt$hlgt <- hlgt_hlt_labels
tidy_rr_pt$hlt <- hlt_pt_labels
tidy_rr_pt$hlgt <- hlgt_pt_labels

object <- list(tidy_rr_pt = tidy_rr_pt,
               tidy_rr_hlt = tidy_rr_hlt,
               tidy_rr_hlgt = tidy_rr_hlgt,
               tidy_rr_soc = tidy_rr_soc)
return(object)
}
