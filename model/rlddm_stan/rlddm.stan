data {
  int<lower=1> N;      // Number of subjects
  real minRT[N];       // minimum RT for each subject of the observed data
  int first[N];        // first trial of subject
  int last[N];         // last trial of subject
  int<lower=1> T;      // Number of observations
  real RTbound;        // lower bound of RT across all subjects (e.g., 0.1 second)
  real iter[T];         // trial of given observation
  int correct[T];      // encodes successful trial
  int incorrect[T];    // encodes unsuccessful trial (inverse of correct)
  real RT[T];          // reaction time
  int value[T];               // value of trial: successful / unsuccessful -> encodes rewards
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
  vector[N] eta_pr[2];
  vector[N] a_mod_pr;
  vector[N] v_mod_pr;
  vector[N] tau_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0>[N] alpha;                       // boundary separation
  vector<lower=-5>[N] eta[2];                     // learning parameter
  vector<lower=-0.5, upper=2>[N] a_mod;         // choice consistency
  vector<lower=0, upper=10>[N] v_mod;                       // scaling parameter
  vector<lower=RTbound, upper=max(minRT)>[N] tau; // nondecision time

  alpha = exp(mu_pr[1] + sigma[1] * alpha_pr); // exp (.): natural exponential of argument.
  eta[1] = exp(mu_pr[2] + sigma[2] * eta_pr[1]);
  eta[2] = exp(mu_pr[3] + sigma[3] * eta_pr[2]);
  a_mod = exp(mu_pr[4] + sigma[4] * a_mod_pr);
  v_mod = exp(mu_pr[5] + sigma[5] * v_mod_pr);
  for (i in 1:N) {
    tau[i]  = Phi_approx(mu_pr[6] + sigma[6] * tau_pr[i]) * (minRT[i] - RTbound) + RTbound;
  }
}

model {
  vector[N] ev[2];
  vector[T] delta;
  vector[T] log_lik;
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // Individual parameters for non-centered parameterization
  alpha_pr ~ normal(0, 1);
  for(v in 1:2){
    eta_pr[v]  ~ normal(0, 1);
  }
  a_mod_pr ~ normal(0, 1);
  v_mod_pr ~ normal(0, 1);
  tau_pr   ~ normal(0, 1);
  // Begin subject loop
  // until second last 
  for (s in 1:N-1) {
    ev[first[s],1] = 0;
    ev[first[s],2] = 0;
    for(trial in (first[s]):(last[s]-1)) {
      delta[trial] = (ev[trial,1] - ev[trial,2]) * (v_mod[s]);
      RT[trial] ~  wiener(alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s] ,0.5,delta[trial]);
      log_lik[trial] += wiener_lpdf(RT[trial] | alpha[s] * pow(iter[trial]/10,a_mod[s]),tau[s],0.5,delta[trial]);
      ev[trial+1,correct[trial]] = ev[trial,correct[trial]] + inv_logit(eta[s,correct[trial]] * (value[trial]-ev[trial,correct[trial]]));
      ev[trial+1,incorrect[trial]] = ev[trial,incorrect[trial]] + inv_logit(eta[s,incorrect[trial]] * (value[trial]-ev[trial,incorrect[trial]]));
    }
    delta[last[s]] = (ev[last[s]-1,1] - ev[last[s]-1,2]) * (v_mod[s]);
    RT[last[s]] ~  wiener(alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s] ,0.5,delta[last[s]]);
    log_lik[last[s]] += wiener_lpdf(RT[last[s]] | alpha[s] * pow(iter[last[s]]/10,a_mod[s]),tau[s],0.5,delta[last[s]]);
  }
}
generated quantities {
  vector[N] ev[2];
  // For group level parameters
  real<lower=0> mu_alpha;          // boundary separation
  real<lower=0> mu_eta1;    // learning rate
  real<lower=0> mu_eta2; 
  real mu_a_mod;                  // biundary separation modification
  real mu_v_mod;                  // drift rate modification
  real<lower=RTbound, upper=max(minRT)> mu_tau; // nondecision time

  // Assign group level parameter values
  mu_alpha = exp(mu_pr[1]);
  mu_eta1 = exp(mu_pr[2]);
  mu_eta2 = exp(mu_pr[3]);
  mu_a_mod =  exp(mu_pr[4]);
  mu_v_mod =  exp(mu_pr[5]);
  mu_tau = Phi_approx(mu_pr[6]) * (mean(minRT)-RTbound) + RTbound;
}
