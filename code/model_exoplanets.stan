data {
  
  int<lower=1> N;       // number of observations
  int y[N];             // dependent variable  
  vector[N] period;     
  vector[N] depth;      
  vector[N] srho;     
  vector[N] fittype_LS_MCMC;
  vector[N] fittype_MCMC;
  vector[N] fittype_LS;
  vector[N] prad;     
  vector[N] teq;      
  vector[N] insol;    
  vector[N] count;    
  vector[N] smet;     
  vector[N] sparprov_stellar;
  vector[N] sparprov_q16;
  vector[N] sparprov_q17;
  vector[N] sparprov_solar;
  vector[N] ra;     
  vector[N] dec;    
  vector[N] fwm_srao;   
  vector[N] fwm_sdeco;  
  vector[N] dikco_mra;  
  vector[N] dikco_mdec; 
  vector[N] dikco_msky;   
  
}
transformed data{
  
  vector[N] log_period = log(period);     
  vector[N] log_depth = log(depth);      
  vector[N] log_srho = log(srho);     
  vector[N] log_prad = log(prad) ;     
  vector[N] log_teq_c = log(teq) - mean(log(teq));      
  vector[N] log_insol = log(insol);  
  vector[N] period_c = period-mean(period);  
  vector[N] depth_c = depth-mean(depth);  
  

}
parameters{
  
  vector[24] beta;
}
model {
  // Log-likelihood
  target += bernoulli_logit_lpmf(y | beta[1] + beta[2]*log_period + beta[3]*log_depth + 
                                     beta[4]*log_srho + beta[5]*fittype_LS_MCMC +
                                     beta[6]*fittype_MCMC + beta[7]*fittype_LS +
                                     beta[8]*log_prad +
                                     beta[9]*log_teq_c + beta[10]*log_insol + beta[11]*count +
                                     beta[12]*smet + beta[13]*sparprov_stellar +
                                     beta[14]*sparprov_stellar + beta[15]*sparprov_q16 +
                                     beta[16]*sparprov_q17 + beta[17]*sparprov_solar + beta[18]*ra +
                                     beta[19]*dec + beta[20]*fwm_srao + beta[21]*fwm_sdeco +
                                     beta[22]*dikco_mra + beta[23]*dikco_mdec + beta[24]*dikco_msky);
  // Log-priors
  target += normal_lpdf(beta | 0, 10);
  
}
// generated quantities{
//   vector[N] log_lik;
//   vector[N] yrep;
//   for(i in 1:N){
//     real theta = beta[1] + beta[2]*dist_c[i] + beta[3]*arsenic_c[i] + beta[4]*dist_arsenic_c[i];
//     log_lik[i] = bernoulli_logit_lpmf(y | theta);
//     yrep[i] = bernoulli_logit_rng(theta);
//   }
// }
