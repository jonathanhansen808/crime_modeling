
data {
  int<lower=0> N;
  int<lower = 0> J; 
  int<lower = 1, upper = J> group_id[N];
  int<lower = 1, upper = J> group_id_test[N];
  
  vector[N] y;
  vector[N] x1; 
  vector[N] x2; 
  vector[N] x3; 
  vector[N] x4; 
  
  vector[N] y_test;
  vector[N] x1_test; 
  vector[N] x2_test; 
  vector[N] x3_test; 
  vector[N] x4_test; 

  real mu_alpha;
  real sd_alpha; 
  real mu_beta1;
  real sd_beta1;
  real mu_beta2;
  real sd_beta2;
  real mu_beta3;
  real sd_beta3;
  real mu_beta4;
  real sd_beta4;
  
  real<lower = 0> nu;
  real<lower = 0> A;
  
  int<lower = 0> n_test;
  int n_grid; 
  vector[n_grid] x_grid; 
}

parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real<lower = 0> tau;
  real alpha_bar;
  vector[J] alpha_aux; // auxiliary parameter for group-specific intercepts 
}

transformed parameters { 
    vector[J] alpha;

    for(j in 1:J){
      alpha[j] = alpha_bar + tau * alpha_aux[j];
  }
}

model {
  vector[N] mu; 
  tau ~ student_t(nu, 0.0, A);
  alpha_bar ~ normal(mu_alpha, sd_alpha);
  beta1 ~ normal(mu_beta1, sd_beta1);
  beta2 ~ normal(mu_beta2, sd_beta2);
  beta3 ~ normal(mu_beta3, sd_beta3);
  beta4 ~ normal(mu_beta4, sd_beta4);
  
  for(j in 1:J){
    alpha_aux[j] ~ normal(0,1);
  }

  for(i in 1:N){
    mu[i] = alpha[group_id[i]] + (x1[i] * beta1) + (x2[i] * beta2) + (x3[i] * beta3) + (x4[i] * beta4);
  }

  y ~ normal(mu, tau);
}

generated quantities { 
    vector[n_test] ystar;
    vector[n_grid] post_line;
    matrix[J+1,n_grid] prob_grid; // prob. that y = 1 at each grid point
    real alpha_new;  // alpha for the new group
  
    alpha_new = normal_rng(alpha_bar,tau);

    for(i in 1:n_grid){
    post_line[i] = alpha[i] + x_grid[i] * beta1 +  x_grid[i] * beta2 +  x_grid[i] * beta3 +  x_grid[i] * beta4;
   }
    
    for(i in 1:n_test){
      ystar[i] =  normal_rng(alpha[group_id_test[i]] + x1_test[i] * beta1 +  x2_test[i] * beta2 +  x3_test[i] * beta3 +  x4_test[i] * beta4, tau);
    }
    
    for(i in 1:n_grid){
      for(j in 1:J){
      prob_grid[j,i] = normal_rng(alpha[i] + x1_test[i] * beta1 +  x2_test[i] * beta2 +  x3_test[i] * beta3 +  x4_test[i] * beta4, tau); // prediction for j-th kicker at i-th grid point
    }
    prob_grid[J+1,i] = normal_rng(alpha_new + x1_test[i] * beta1 +  x2_test[i] * beta2 +  x3_test[i] * beta3 +  x4_test[i] * beta4, tau); // prediction for new group at i-th grid point
  }
  

}

// sigma in generated when declared and in likelihood
