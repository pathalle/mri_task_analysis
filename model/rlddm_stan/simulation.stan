data {
  int<lower=1> N;      // Number of subjects
  real minRT[N];       // minimum RT for each subject of the observed data
  int first[N];        // first trial of subject
  int last[N];         // last trial of subject
  int<lower=1> T;      // Number of observations
  real RTbound;        // lower bound of RT across all subjects (e.g., 0.1 second)
  real iter[T];        // trial of given observation
  int response[T];      // encodes successful trial [1: lower bound (incorrect), 2: upper bound(correct)]
  real RT[T];          // reaction time
  int value[T];        // value of trial: successful / unsuccessful -> encodes rewards
  int stim_assoc[T];   // index of associated sound-symbol pair
  int stim_nassoc[T];  // index of presented non-associated symbol
  int n_stims[N];      // number of items learned by each subject (represents # blocks)
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
}

parameters {
  // alpha (a): Boundary separation or Speed-accuracy trade-off 
  // delta (v): Drift rate 
  // tau (ter): Nondecision time + Motor response time + encoding time
  // a_mod: modulator for decision boundary
  // v_mod: modulator for drift diffusion rate

  // Hyper-parameters
  vector[6] mu_pr;
  vector<lower=0>[6] sigma;

  // Subject-level raw parameters
  vector[N] alpha_pr;
  //vector[N] eta_pr_pos;
  //vector[N] eta_pr_neg;
  vector[N] a_mod_pr;
  vector[N] v_mod_pr;
  vector[N] tau_pr;
  
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0>[N] alpha;                       // boundary separation
  //vector<lower=0>[N] eta_pos; // learning parameter for upper boundary
  //vector[N] eta_neg; // learning parameter for lower boundary
  vector[N] a_mod;           // choice consistency
  vector<lower=0, upper=10>[N] v_mod;             // scaling parameter
  vector<lower=RTbound, upper=max(minRT)>[N] tau; // nondecision time

  alpha = exp(mu_pr[1] + sigma[1] * alpha_pr); //
  //eta_neg = exp(mu_pr[2] + sigma[2] * eta_pr_pos);
  //eta_pos = exp(mu_pr[3] + sigma[3] * eta_pr_neg);
  a_mod = exp(mu_pr[4] + sigma[4] * a_mod_pr);
  v_mod = exp(mu_pr[5] + sigma[5] * v_mod_pr);
  for (s in 1:N) {
    tau[s]  = Phi_approx(mu_pr[6] + sigma[6] * tau_pr[s]) * (minRT[s] - RTbound) + RTbound;
    tau[s] = fabs(tau[s]);
  }
}

model {
  real eta_neg = logit(0.07);
  real eta_pos = logit(0.07);
  vector[T] log_lik;
  real ev[T,max(n_stims)];
  vector[T] delta;
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // Individual parameters
  alpha_pr ~ normal(0, 1);
  //eta_pr_pos ~ normal(0,1);
  //eta_pr_neg ~ normal(0,1);
  a_mod_pr ~ normal(0, 1);
  v_mod_pr ~ normal(0, 1);
  tau_pr   ~ normal(0, 1);
  // Begin subject loop
  // until second last 
  for (s in 1:N) {
    for(a in 1:n_stims[s]){
      // ev for pos values
      ev[first[s],a] = 0.5;
    }
    for(trial in (first[s]):(last[s]-1)) {
      for(a in 1:n_stims[s]){
        ev[trial+1,a] = ev[trial,a];
      }
      delta[trial] = (ev[trial,stim_assoc[trial]] + ev[trial,stim_nassoc[trial]])/2 * v_mod[s];
      // if lower bound
      if (response[trial]==1){
        if(run_estimation==1){
          RT[trial] ~  wiener(alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s] ,0.5,-(delta[trial]));
        }
        log_lik[trial] = wiener_lpdf(RT[trial] | alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s],0.5,-(delta[trial]));
        ev[trial+1,stim_nassoc[trial]] = ev[trial,stim_nassoc[trial]] - (inv_logit(eta_neg) * (value[trial]-(1-ev[trial,stim_nassoc[trial]])));
        ev[trial+1,stim_assoc[trial]] = ev[trial,stim_assoc[trial]] - (inv_logit(eta_neg) * (value[trial]-ev[trial,stim_assoc[trial]]));
      }
      // if upper bound (resp = 2)
      else{
        if(run_estimation==1){
          RT[trial] ~  wiener(alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s] ,0.5,delta[trial]);
        }
        log_lik[trial] = wiener_lpdf(RT[trial] | alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s],0.5,delta[trial]);
        ev[trial+1,stim_nassoc[trial]] = ev[trial,stim_nassoc[trial]] + (inv_logit(eta_pos) * (value[trial]-(1-ev[trial,stim_nassoc[trial]])));
        ev[trial+1,stim_assoc[trial]] = ev[trial,stim_assoc[trial]] + (inv_logit(eta_pos) * (value[trial]-ev[trial,stim_assoc[trial]]));
      }
    }
    // in last cycle, don't update anymore
    delta[last[s]] = (ev[last[s]-1,stim_assoc[last[s]]] - ev[last[s]-1,stim_nassoc[last[s]]])/2 * v_mod[s];
    if (response[last[s]]==1){
      if(run_estimation==1){
        RT[last[s]] ~  wiener(alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s] ,0.5,-(delta[last[s]]));
      log_lik[last[s]] = wiener_lpdf(RT[last[s]] | alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s],0.5,-(delta[last[s]]));
      }
    }
    if (response[last[s]]==2){
      if(run_estimation==1){
        RT[last[s]] ~  wiener(alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s] ,0.5,delta[last[s]]);
      }
      log_lik[last[s]] = wiener_lpdf(RT[last[s]] | alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s],0.5,delta[last[s]]);
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0> mu_alpha;                  // boundary separation
  real<lower=0> mu_eta_neg;                 // learning rate lower
  real<lower=0> mu_eta_pos;                // learning rate upper
  real<lower=0> mu_a_mod;                  // boundary separation modification
  real<lower=0> mu_v_mod;                  // drift rate modification
  real<lower=RTbound, upper=max(minRT)> mu_tau; // nondecision time
  
  real ev_hat[T,max(n_stims)];
  real pe_hat[T];
  real assoc_active_pair[T];
  real assoc_inactive_pair[T];
  vector[T] delta_hat;
  
  real eta_neg_gen = logit(0.07);
  real eta_pos_gen = logit(0.07);
  
  // Assign group level parameter values
  mu_alpha = exp(mu_pr[1]);
  mu_eta_neg = exp(mu_pr[2]);
  mu_eta_pos = exp(mu_pr[3]);
  mu_a_mod =  exp(mu_pr[4]);
  mu_v_mod =  exp(mu_pr[5]);
  mu_tau = Phi_approx(mu_pr[6]) * (mean(minRT)-RTbound) + RTbound;
  
  for (s in 1:N){
    for(a in 1:n_stims[s]){
      // ev for pos values
      ev_hat[first[s],a] = 0.5;
    }
    assoc_active_pair[first[s]] = 0.5;
    assoc_inactive_pair[first[s]] = 0.5;
    for(trial in (first[s]):(last[s]-1)) {
      for(a in 1:n_stims[s]){
        ev_hat[trial+1,a] = ev_hat[trial,a];
      }
      delta_hat[trial] = (ev_hat[trial,stim_assoc[trial]] + ev_hat[trial,stim_nassoc[trial]])/2 * v_mod[s];
      assoc_active_pair[trial] = ev_hat[trial,stim_assoc[trial]];
      assoc_inactive_pair[trial] = ev_hat[trial,stim_nassoc[trial]];
      pe_hat[trial] = value[trial]-(ev_hat[trial,stim_assoc[trial]]);
      // if lower bound
      if (response[trial]==1){
        ev_hat[trial+1,stim_nassoc[trial]] = ev_hat[trial,stim_nassoc[trial]] - (inv_logit(eta_neg_gen) * (value[trial]-(1-ev_hat[trial,stim_nassoc[trial]])));
        ev_hat[trial+1,stim_assoc[trial]] = ev_hat[trial,stim_assoc[trial]] - (inv_logit(eta_neg_gen) * (value[trial]-ev_hat[trial,stim_assoc[trial]]));
      }
      // if upper bound (resp = 2)
      else{

        ev_hat[trial+1,stim_nassoc[trial]] = ev_hat[trial,stim_nassoc[trial]] + (inv_logit(eta_pos_gen) * (value[trial]-(1-ev_hat[trial,stim_nassoc[trial]])));
        ev_hat[trial+1,stim_assoc[trial]] = ev_hat[trial,stim_assoc[trial]] + (inv_logit(eta_pos_gen) * (value[trial]-ev_hat[trial,stim_assoc[trial]]));
      }
    pe_hat[last[s]] = value[last[s]]-(ev_hat[last[s],stim_nassoc[last[s]]]);
    assoc_active_pair[last[s]] = ev_hat[last[s],stim_assoc[last[s]]];
    assoc_inactive_pair[last[s]] = ev_hat[last[s],stim_nassoc[last[s]]];
    delta_hat[last[s]] = (ev_hat[last[s]-1,stim_assoc[last[s]]] + ev_hat[last[s]-1,stim_nassoc[last[s]]])/2 * v_mod[s];
    }
  }
}

