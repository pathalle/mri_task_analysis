// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists

/* missing parameters for rlddm

// learning rate
etag_mu  ~  dunif(-5,5) 
etag_sd  ~  dunif(0.001,5)
etag_tau <- pow(etag_sd,-2)

// choice consistency (in expression for drift rate)
ig_mu    ~  dunif(-0.5,0.5)
ig_sd    ~  dunif(0.001,5)
ig_tau   <- pow(ig_sd,-2)

// scaling parameter?
mg_mu    ~  dunif(0,10)
mg_sd    ~  dunif(0.001,5)
mg_tau   <- pow(mg_sd,-2)


for (s in 1:S) {
     //Assign starting values to ev[trial,option]. 
     //first is a two-dimensional-array identifying first trial for each subject in each group.
     //ev: expectancy valence!
          ev[first[s],1] <- 0
          ev[first[s],2] <- 0
     }
     //Run through trials
     //last is a two-dimensional-array identifying last trial for each subject in each group.
for (trial in (first[s]):(last[s]-1)) {  
     //Calculate drift rate parameter as ev-delta multiplited by m
     v[trial] <- (ev[trial,1] - ev[trial,2]) * (m[s]) // Define drift rate as 
     //Estimate likelihood of choices response time with dwiener. save log_likelihood for LOO
     // The JAGS Wiener module differentiates choices toward upper and lower boundary by the sign of the RT. Therefore, the RT of all choices in favor the suboptimal options (B, D and F) are set to be negative by multiplying the nal RT with -1. 
    // dwiener(alpha:a,tau:ter,beta:z,delta:drift)
     RT[trial] ~ dwiener(a[s] * pow(iter[trial]/10,i[s]),t[s] ,0.5,v[trial])
     log_lik[trial] <- logdensity.wiener(RT[trial], a[s] * pow(iter[trial]/10,i[s]),t[s],0.5,v[trial])
  //Update ev-values for next trial.
     ev[trial+1,choice[trial]] <- ev[trial,choice[trial]] + ilogit(eta[s,valence[trial]] * (value[trial]-ev[trial,choice[trial]]))
     ev[trial+1,nonchoice[trial]] <- ev[trial,nonchoice[trial]] 
     ev[trial+1,1] <- ev[trial,1]
     ev[trial+1,2] <- ev[trial,2]
}
          //EV-values are not updated in last trial
for (trial in last[s]) {
     v[trial] <- (ev[trial,1] - ev[trial,2]) * m[s]
     RT[trial] ~ dwiener(a[s] * pow(iter[trial]/10,i[s]),t[s],0.5,v[trial])
     log_lik[trial] <- logdensity.wiener(RT[trial], a[s]  * pow(iter[trial]/10),i[s],t[s],0.5,v[trial])
}


*/



data {
  int<lower=1> N;      // Number of subjects
  int<lower=0> Nu_max; // Max (across subjects) number of upper boundary responses
  int<lower=0> Nl_max; // Max (across subjects) number of lower boundary responses
  int<lower=0> Nu[N];  // Number of upper boundary responses for each subj
  int<lower=0> Nl[N];  // Number of lower boundary responses for each subj
  real RTu[N, Nu_max];  // upper boundary response times
  real RTl[N, Nl_max];  // lower boundary response times
  real minRT[N];       // minimum RT for each subject of the observed data
  real RTbound;        // lower bound of RT across all subjects (e.g., 0.1 second)
  int first[N,1];      // first trial of choice 1 for each subject 
  int first[N,2];      // first trial of choice 2 for each subject 
  int last[N,1];       // last trial of choice 1 for each subject 
  int last[N,2];       // last trial of choice 2 for each subject 
}

parameters {
  // parameters of the DDM (parameter names in Ratcliffs DDM), from https://github.com/gbiele/stan_wiener_test/blob/master/stan_wiener_test.R
  // also see: https://groups.google.com/forum///!searchin/stan-users/wiener%7Csort:relevance/stan-users/-6wJfA-t2cQ/Q8HS-DXgBgAJ
  // alpha (a): Boundary separation or Speed-accuracy trade-off (high alpha means high accuracy). alpha > 0
  // beta (b): Initial bias Bias for either response (beta > 0.5 means bias towards "upper" response 'A'). 0 < beta < 1
  // delta (v): Drift rate Quality of the stimulus (delta close to 0 means ambiguous stimulus or weak ability). 0 < delta
  // tau (ter): Nondecision time + Motor response time + encoding time (high means slow encoding, execution). 0 < ter (in seconds)
  // choice consistency (i)
  // scaling parameter (m)
  ///* upper boundary of tau must be smaller than minimum RT
  //to avoid zero likelihood for fast responses.
  //tau can for physiological reasone not be faster than 0.1 s.*/


  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[7] mu_pr;
  vector<lower=0>[7] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[N] alpha_pr;
  vector[N] beta_pr;
  //vector[N] delta_pr;
  vector[N] tau_pr;
  vector[N] eta_pr;
  vector[N] i_pr;
  vector[N] m_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0>[N]         alpha; // boundary separation
  vector<lower=0, upper=1>[N] beta;  // initial bias
  //vector<lower=0>[N]         delta; // drift rate
  vector<lower=RTbound, upper=max(minRT)>[N] tau; // nondecision time
  vector<lower=-5>[N]        eta; // learning parameter
  vector<lower=-0.5>[N]         cc; // choice consistency
  vector<lower=0>[N]        m; // scaling parameter

  for (i in 1:N) {
    beta[i] = Phi_approx(mu_pr[2] + sigma[2] * beta_pr[i]); // phi_approx: inverse CDF of unit normal
    tau[i]  = Phi_approx(mu_pr[4] + sigma[4] * tau_pr[i]) * (minRT[i] - RTbound) + RTbound;
  }
  // I assume this will assign the same alpha value to each subject value
  alpha = exp(mu_pr[1] + sigma[1] * alpha_pr); // exp (.) return the natural exponential of the specified argument.
  //delta = exp(mu_pr[3] + sigma[3] * delta_pr); 
  eta = exp(mu_pr[5] + sigma[5] * i_pr);
  cc = exp(mu_pr[6] + sigma[6] * i_pr);
  m = exp(mu_pr[7] + sigma[7] * m_pr);
}

model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // Individual parameters for non-centered parameterization
  alpha_pr ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);
  //delta_pr ~ normal(0, 1);
  tau_pr   ~ normal(0, 1);
  eta_pr   ~ normal(0, 1);
  i_pr     ~ normal(0, 1);
  m_pr     ~ normal(0, 1);

  // Begin subject loop
  // until second last 
for(i in 1:N){
  ev[first[i],1] = 0;
  ev[first[i],2] = 0;
}

  for (i in 1:N-1) {
    delta[trial] = (ev[trial,1]-ev[trial,2])*m[s];

    // Response time distributed along wiener first passage time distribution
    RTu[i, :Nu[i]] ~ wiener(alpha[i], tau[i], beta[i], delta[i]);
    RTl[i, :Nl[i]] ~ wiener(alpha[i], tau[i], 1-beta[i], -delta[i]);

  } // end of subject loop
}

generated quantities {
  // For group level parameters
  real<lower=0>         mu_alpha; // boundary separation
  real<lower=0, upper=1> mu_beta;  // initial bias
  real<lower=0>         mu_delta; // drift rate
  real<lower=RTbound, upper=max(minRT)> mu_tau; // nondecision time
  real 

  // For log likelihood calculation
  real log_lik[N];

  // Assign group level parameter values
  mu_alpha = exp(mu_pr[1]);
  mu_beta  = Phi_approx(mu_pr[2]);
  mu_delta = exp(mu_pr[3]);
  mu_tau   = Phi_approx(mu_pr[4]) * (mean(minRT)-RTbound) + RTbound;

  { // local section, this saves time and space
    // Begin subject loop
    for (i in 1:N) {
      log_lik[i] = wiener_lpdf(RTu[i, :Nu[i]] | alpha[i], tau[i], beta[i], delta[i]);
      log_lik[i] += wiener_lpdf(RTl[i, :Nl[i]] | alpha[i], tau[i], 1-beta[i], -delta[i]);
    }
  }
}

