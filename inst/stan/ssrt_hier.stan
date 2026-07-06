// =============================================================================
// ssrt_hier.stan
// Hierarchical (multi-subject) independent horse-race model with ex-Gaussian
// go- and stop-process distributions. Non-centred parameterisation on the log
// scale for all six subject-level parameters (avoids funnel pathologies and
// keeps every quantity positive without needing <lower=0> Jacobian warnings).
//
// See ssrt_single.stan for details on the likelihood / log_p_inhibit().
//
// Requires Stan >= 2.26 (array[] syntax, log1m, exp_mod_normal_lccdf).
// =============================================================================

functions {
  real log_p_inhibit(real ssd,
                      real mu_go,   real sigma_go,   real lambda_go,
                      real mu_stop, real sigma_stop, real lambda_stop,
                      int  n_quad) {
    real tau_stop = 1.0 / lambda_stop;
    real upr      = mu_stop + tau_stop
                    + 5.0 * sqrt(square(sigma_stop) + square(tau_stop));
    real dt = upr / n_quad;
    real lp = negative_infinity();
    for (k in 1:n_quad) {
      real t     = (k - 0.5) * dt;
      real log_f = exp_mod_normal_lpdf( t        | mu_stop, sigma_stop, lambda_stop);
      real log_s = exp_mod_normal_lccdf(ssd + t   | mu_go,   sigma_go,   lambda_go);
      lp = log_sum_exp(lp, log_f + log_s + log(dt));
    }
    return lp;
  }
}

data {
  int<lower=1> J;   // number of subjects

  // Go trials
  int<lower=1>                       N_go;
  array[N_go] int<lower=1, upper=J>  go_subj;
  vector<lower=0>[N_go]              go_rt;

  // Signal-respond trials
  int<lower=0>                       N_sr;
  array[N_sr] int<lower=1, upper=J>  sr_subj;
  vector<lower=0>[N_sr]              sr_rt;
  vector<lower=0>[N_sr]              sr_ssd;

  // Inhibited trials
  int<lower=0>                       N_inh;
  array[N_inh] int<lower=1, upper=J> inh_subj;
  vector<lower=0>[N_inh]             inh_ssd;

  int<lower=10, upper=500> n_quad;
  int<lower=0,  upper=1>   use_tf;
}

parameters {
  // Population-level means on the log scale. HARD BOUNDS keep both Stan's
  // default initialisation and every HMC proposal inside a physically sensible
  // range for reaction times in MILLISECONDS. Without these bounds, Stan's
  // default init (uniform on (-2, 2) in the unconstrained space) starts e.g.
  // mu_go at exp(0) = 1 ms, making every go-trial density log(0) and, once the
  // sd * z term is added, overflowing exp() to Inf. The bounds below make that
  // impossible: the implied natural-scale value is always finite and sane.
  real<lower=3.912, upper=7.313> mu_go_pop_log;      // mu_go    in [50, 1500] ms
  real<lower=1.099, upper=6.215> sigma_go_pop_log;   // sigma_go in [3, 500] ms
  real<lower=1.099, upper=6.685> tau_go_pop_log;     // tau_go   in [3, 800] ms
  real<lower=2.996, upper=6.685> mu_stop_pop_log;    // mu_stop  in [20, 800] ms
  real<lower=0.693, upper=5.991> sigma_stop_pop_log; // sigma_stop in [2, 400] ms
  real<lower=0.693, upper=6.215> tau_stop_pop_log;   // tau_stop in [2, 500] ms

  // Population-level SDs of the log random effects. Bounded above so the
  // non-centred term exp(pop_log + sd * z) cannot overflow: with sd <= 1.5 and
  // pop_log <= 7.313, the exponent stays far below the ~709 overflow threshold.
  real<lower=0, upper=1.5> sd_mu_go;     real<lower=0, upper=1.5> sd_sigma_go;     real<lower=0, upper=1.5> sd_tau_go;
  real<lower=0, upper=1.5> sd_mu_stop;   real<lower=0, upper=1.5> sd_sigma_stop;   real<lower=0, upper=1.5> sd_tau_stop;

  // Non-centred subject-level deviations (standard-normal z-scores). Bounded
  // to +/-8 SDs so an extreme HMC proposal during warmup cannot push
  // exp(pop_log + sd * z) to Inf (which caused "Location/Scale parameter is
  // inf" and "Inv_scale parameter is 0" proposal rejections). No real subject
  // deviates 8 population-SDs, so the standard-normal prior is essentially
  // untouched (P(|z| > 8) ~ 1e-15) and the non-centred parameterisation is
  // preserved. With |z| <= 8, sd <= 1.5 and pop_log <= 7.313, the exponent
  // stays <= 19.3, so every subject-level value is finite and positive.
  vector<lower=-8, upper=8>[J] z_mu_go;     vector<lower=-8, upper=8>[J] z_sigma_go;     vector<lower=-8, upper=8>[J] z_tau_go;
  vector<lower=-8, upper=8>[J] z_mu_stop;   vector<lower=-8, upper=8>[J] z_sigma_stop;   vector<lower=-8, upper=8>[J] z_tau_stop;

  // Per-subject trigger-failure probabilities (zero-length when use_tf == 0)
  vector<lower=0, upper=1>[use_tf ? J : 0] p_tf;
}

transformed parameters {
  vector<lower=0>[J] mu_go     = exp(mu_go_pop_log    + sd_mu_go    * z_mu_go);
  vector<lower=0>[J] sigma_go  = exp(sigma_go_pop_log + sd_sigma_go * z_sigma_go);
  vector<lower=0>[J] tau_go    = exp(tau_go_pop_log   + sd_tau_go   * z_tau_go);

  vector<lower=0>[J] mu_stop    = exp(mu_stop_pop_log    + sd_mu_stop    * z_mu_stop);
  vector<lower=0>[J] sigma_stop = exp(sigma_stop_pop_log + sd_sigma_stop * z_sigma_stop);
  vector<lower=0>[J] tau_stop   = exp(tau_stop_pop_log   + sd_tau_stop   * z_tau_stop);

  vector<lower=0>[J] lambda_go   = inv(tau_go);
  vector<lower=0>[J] lambda_stop = inv(tau_stop);
  vector[J]          mean_ssrt   = mu_stop + tau_stop;
}

model {
  // Hyperpriors (centred on typical adult RT/SSRT values, ms)
  mu_go_pop_log    ~ normal(log(400), 0.3);
  sigma_go_pop_log ~ normal(log(50),  0.5);
  tau_go_pop_log   ~ normal(log(100), 0.5);

  mu_stop_pop_log    ~ normal(log(200), 0.5);
  sigma_stop_pop_log ~ normal(log(30),  0.5);
  tau_stop_pop_log   ~ normal(log(80),  0.5);

  // Half-normal-on-log-scale priors for between-subject SDs
  // (~0.3-0.5 on the log scale corresponds to roughly 30-65% multiplicative
  //  variability between subjects)
  sd_mu_go    ~ normal(0, 0.3);   sd_sigma_go    ~ normal(0, 0.4);   sd_tau_go    ~ normal(0, 0.4);
  sd_mu_stop  ~ normal(0, 0.4);   sd_sigma_stop  ~ normal(0, 0.5);   sd_tau_stop  ~ normal(0, 0.5);

  // Non-centred standard-normal priors
  z_mu_go    ~ std_normal();   z_sigma_go    ~ std_normal();   z_tau_go    ~ std_normal();
  z_mu_stop  ~ std_normal();   z_sigma_stop  ~ std_normal();   z_tau_stop  ~ std_normal();

  if (use_tf) {
    p_tf ~ beta(1.5, 8);
  }

  // Likelihood
  for (i in 1:N_go) {
    int j = go_subj[i];
    target += exp_mod_normal_lpdf(go_rt[i] | mu_go[j], sigma_go[j], lambda_go[j]);
  }

  for (i in 1:N_sr) {
    int j = sr_subj[i];
    real lp_go   = exp_mod_normal_lpdf( sr_rt[i]             | mu_go[j],   sigma_go[j],   lambda_go[j]);
    real lp_surv = exp_mod_normal_lccdf(sr_rt[i] - sr_ssd[i] | mu_stop[j], sigma_stop[j], lambda_stop[j]);
    if (use_tf) {
      target += log_sum_exp(log1m(p_tf[j]) + lp_go + lp_surv,
                             log(p_tf[j])   + lp_go);
    } else {
      target += lp_go + lp_surv;
    }
  }

  for (i in 1:N_inh) {
    int j = inh_subj[i];
    real lp_inh = log_p_inhibit(inh_ssd[i],
                                 mu_go[j],   sigma_go[j],   lambda_go[j],
                                 mu_stop[j], sigma_stop[j], lambda_stop[j],
                                 n_quad);
    if (use_tf) {
      target += log1m(p_tf[j]) + lp_inh;
    } else {
      target += lp_inh;
    }
  }
}

generated quantities {
  // Population-level mean SSRT (natural scale)
  real pop_mean_ssrt = exp(mu_stop_pop_log) + exp(tau_stop_pop_log);

  // Per-subject posterior-predictive draws
  array[J] real ssrt_pred;
  array[J] real go_rt_pred;
  for (j in 1:J) {
    ssrt_pred[j]  = normal_rng(mu_stop[j], sigma_stop[j]) + exponential_rng(lambda_stop[j]);
    go_rt_pred[j] = normal_rng(mu_go[j],   sigma_go[j])   + exponential_rng(lambda_go[j]);
  }
}
