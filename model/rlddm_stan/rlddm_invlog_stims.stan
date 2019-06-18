data {
  int<lower=1> N;      // Number of subjects
  real minRT[N];       // minimum RT for each subject of the observed data
  int first[N];        // first trial of subject
  int last[N];         // last trial of subject
  int<lower=1> T;      // Number of observations
  real RTbound;        // lower bound of RT across all subjects (e.g., 0.1 second)
  real iter[T];        // trial of given observation
  int response[T];      // encodes successful trial [1: lower bound (incorrect), 2: upper bound(correct)]
  int nonresponse[T];      // encodes unsuccessful trial
  real RT[T];          // reaction time
  int value[T];        // value of trial: successful / unsuccessful -> encodes rewards
}

parameters {
  // alpha (a): Boundary separation or Speed-accuracy trade-off (high alpha means high accuracy). alpha > 0
  // delta (v): Drift rate Quality of the stimulus (delta close to 0 means ambiguous stimulus or weak ability). 0 < delta
  // tau (ter): Nondecision time + Motor response time + encoding time (high means slow encoding, execution). 0 < ter (in seconds)
  // a_mod: 
  // v_mod: 
  /// upper boundary of tau must be smaller than minimum RT to avoid zero likelihood for fast responses.
  //tau can for physiological reasone not be faster than 0.1 s.*/

  // Hyper(group)-parameters
  vector[6] mu_pr;
  vector<lower=0>[6] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[N] alpha_pr;
  vector[N] eta_pr_pos;
  vector[N] eta_pr_neg;
  vector[N] a_mod_pr;
  vector[N] v_mod_pr;
  vector[N] tau_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0>[N] alpha;                       // boundary separation
  vector[N] eta_pos; // learning parameter for upper boundary
  vector[N] eta_neg; // learning parameter for lower boundary
  vector[N] a_mod;           // choice consistency
  vector<lower=0, upper=10>[N] v_mod;             // scaling parameter
  vector<lower=RTbound, upper=max(minRT)>[N] tau; // nondecision time

  alpha = exp(mu_pr[1] + sigma[1] * alpha_pr); // exp (.): natural exponential of argument.
  eta_neg = exp(mu_pr[2] + sigma[2] * eta_pr_pos);
  eta_pos = exp(mu_pr[3] + sigma[3] * eta_pr_neg);
  a_mod = exp(mu_pr[4] + sigma[4] * a_mod_pr);
  v_mod = exp(mu_pr[5] + sigma[5] * v_mod_pr);
  for (s in 1:N) {
    tau[s]  = Phi_approx(mu_pr[6] + sigma[6] * tau_pr[s]) * (minRT[s] - RTbound) + RTbound;
    tau[s] = fabs(tau[s]);
  }
}

model {
  real ev[T,2];
  vector[T] delta;
  vector[T] log_lik;
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // Individual parameters for non-centered parameterization
  alpha_pr ~ normal(0, 1);
  eta_pr_pos ~ normal(0,1);
  eta_pr_neg ~ normal(0,1);
  a_mod_pr ~ normal(0, 1);
  v_mod_pr ~ normal(0, 1);
  tau_pr   ~ normal(0, 1);
  // Begin subject loop
  // until second last 
  for (s in 1:N) {
    // ev for pos values
    ev[first[s],1] = 0.5;
    // ev for neg values
    ev[first[s],2] = 0.5;
    for(trial in (first[s]):(last[s]-1)) {
      delta[trial] = (ev[trial,2] - ev[trial,1]) * v_mod[s];
      // if lower bound
      if (response[trial]==1){
        RT[trial] ~  wiener(alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s] ,0.5,-delta[trial]);
        log_lik[trial] = wiener_lpdf(RT[trial] | alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s],0.5,-delta[trial]);
        // if the anwswer is wrong, update ev for lower value with neg lr, ev for positive answer stays the same
        ev[trial+1,1] = ev[trial,1] + inv_logit(eta_neg[s]) * (value[trial]-ev[trial,1]);
        ev[trial+1,2] = ev[trial,2];
      }
      // if upper bound (resp = 2)
      else{
        RT[trial] ~  wiener(alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s] ,0.5,delta[trial]);
        log_lik[trial] = wiener_lpdf(RT[trial] | alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s],0.5,delta[trial]);
        // if the anwswer is correct, update ev for upper value with pos lr, ev for neg answer stays the same
        ev[trial+1,2] = ev[trial,2] + inv_logit(eta_pos[s]) * (value[trial]-ev[trial,2]);
        ev[trial+1,1] = ev[trial,1];
      }
    }
    // in last cycle, don't update anymore
    delta[last[s]] = (ev[last[s]-1,2] - ev[last[s]-1,1]) * v_mod[s];
    if (response[last[s]]==1){
      RT[last[s]] ~  wiener(alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s] ,0.5,-delta[last[s]]);
      log_lik[last[s]] = wiener_lpdf(RT[last[s]] | alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s],0.5,-delta[last[s]]);
    }
    else{
      RT[last[s]] ~  wiener(alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s] ,0.5,delta[last[s]]);
      log_lik[last[s]] = wiener_lpdf(RT[last[s]] | alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s],0.5,delta[last[s]]);
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0> mu_alpha;          // boundary separation
  real<lower=0> mu_eta_neg;    // learning rate
  real<lower=0> mu_eta_pos; 
  real<lower=0> mu_a_mod;                  // biundary separation modification
  real<lower=0> mu_v_mod;                  // drift rate modification
  real<lower=RTbound, upper=max(minRT)> mu_tau; // nondecision time
  
  // Assign group level parameter values
  mu_alpha = exp(mu_pr[1]);
  mu_eta_neg = exp(mu_pr[2]);
  mu_eta_pos = exp(mu_pr[3]);
  mu_a_mod =  exp(mu_pr[4]);
  mu_v_mod =  exp(mu_pr[5]);
  mu_tau = Phi_approx(mu_pr[6]) * (mean(minRT)-RTbound) + RTbound;

}
