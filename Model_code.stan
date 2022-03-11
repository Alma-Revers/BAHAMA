data {
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
}