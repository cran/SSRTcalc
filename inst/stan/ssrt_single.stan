// =============================================================================
// ssrt_single.stan
// Single-subject independent horse-race model with ex-Gaussian go- and
// stop-process distributions, following Matzke et al. (2013, JEP:General).
//
// Stan's exp_mod_normal(mu, sigma, lambda) parameterises the ex-Gaussian with
// lambda = 1/tau (RATE, not mean). Mean = mu + tau; Variance = sigma^2 + tau^2.
//
// Trial types and their likelihood contributions:
//   Go              : f_go(rt)
//   Signal-respond  : f_go(rt) * S_stop(rt - SSD)
//   Inhibited       : P(inhibit | SSD) = integral_0^inf S_go(SSD+t) f_stop(t) dt
//                      (computed via midpoint-rule quadrature, log-sum-exp'd)
//
// Optional trigger-failure parameter p_tf (Matzke et al., 2013): on a fraction
// p_tf of stop trials the stop process is never triggered, so:
//   - SR trials become a 2-component mixture (race active vs. TF)
//   - inhibited trials require the trigger to have fired: P(inh) * (1 - p_tf)
//
// Requires Stan >= 2.26 (array[] syntax, log1m, exp_mod_normal_lccdf).
// =============================================================================

functions {
  /**
   * log P(inhibit | SSD) via midpoint-rule numerical quadrature.
   *
   * log integral_0^upper  S_go(SSD + t) * f_stop(t)  dt
   *
   * Integration upper bound: mean_stop + 5 * SD_stop (~99.99% of stop-RT mass).
   * Resolution: n_quad midpoint rectangles (100 is a good default; raise if
   * sigma_stop/tau_stop become large relative to the bound during warmup).
   */
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
  // Go trials
  int<lower=1>           N_go;
  vector<lower=0>[N_go]  go_rt;       // ms

  // Signal-respond trials (responded despite the stop signal)
  int<lower=0>           N_sr;
  vector<lower=0>[N_sr]  sr_rt;       // ms
  vector<lower=0>[N_sr]  sr_ssd;      // stop-signal delay, ms

  // Inhibited (successful stop) trials
  int<lower=0>           N_inh;
  vector<lower=0>[N_inh] inh_ssd;     // ms

  // Settings
  int<lower=10, upper=500> n_quad;    // quadrature resolution
  int<lower=0,  upper=1>   use_tf;    // 1 = include trigger-failure parameter
}

parameters {
  // Go process (ex-Gaussian). Upper bounds prevent exp_mod_normal's location/
  // scale from ever becoming Inf during sampling (all values in ms).
  real<lower=50, upper=2000> mu_go;       // Gaussian mean, ms
  real<lower=1,  upper=800>  sigma_go;    // Gaussian SD, ms
  real<lower=1,  upper=1000> tau_go;      // exponential mean (1/lambda), ms

  // Stop process (ex-Gaussian) -- the SSRT distribution
  real<lower=1,  upper=1000> mu_stop;
  real<lower=1,  upper=500>  sigma_stop;
  real<lower=1,  upper=800>  tau_stop;

  // Trigger-failure probability (zero-length array when use_tf == 0)
  array[use_tf ? 1 : 0] real<lower=0, upper=1> p_tf;
}

transformed parameters {
  real lambda_go   = inv(tau_go);
  real lambda_stop = inv(tau_stop);
  real mean_ssrt   = mu_stop + tau_stop;   // E[SSRT]
  real mean_go_rt  = mu_go   + tau_go;     // E[go RT]
}

model {
  // Weakly informative priors (RT in ms)
  mu_go    ~ normal(400, 150);
  sigma_go ~ normal(50,   40);
  tau_go   ~ normal(100,  80);

  mu_stop    ~ normal(200, 150);
  sigma_stop ~ normal(30,   25);
  tau_stop   ~ normal(80,   60);

  if (use_tf) {
    p_tf[1] ~ beta(1.5, 8);   // trigger-failure rate is typically low
  }

  // 1. Go trials
  for (i in 1:N_go) {
    target += exp_mod_normal_lpdf(go_rt[i] | mu_go, sigma_go, lambda_go);
  }

  // 2. Signal-respond trials: f_go(rt) * S_stop(rt - SSD)
  for (i in 1:N_sr) {
    real lp_go   = exp_mod_normal_lpdf( sr_rt[i]             | mu_go,   sigma_go,   lambda_go);
    real lp_surv = exp_mod_normal_lccdf(sr_rt[i] - sr_ssd[i] | mu_stop, sigma_stop, lambda_stop);

    if (use_tf) {
      // Mixture: race active (1 - p_tf) vs. trigger failed (p_tf)
      target += log_sum_exp(log1m(p_tf[1]) + lp_go + lp_surv,
                             log(p_tf[1])   + lp_go);
    } else {
      target += lp_go + lp_surv;
    }
  }

  // 3. Inhibited trials: P(go loses the race | SSD)
  for (i in 1:N_inh) {
    real lp_inh = log_p_inhibit(inh_ssd[i],
                                 mu_go,   sigma_go,   lambda_go,
                                 mu_stop, sigma_stop, lambda_stop,
                                 n_quad);
    if (use_tf) {
      // Inhibition is only possible if the trigger fired
      target += log1m(p_tf[1]) + lp_inh;
    } else {
      target += lp_inh;
    }
  }
}

generated quantities {
  // Posterior-predictive draws for density/PP-check plots
  real ssrt_pred  = normal_rng(mu_stop, sigma_stop) + exponential_rng(lambda_stop);
  real go_rt_pred = normal_rng(mu_go,   sigma_go)   + exponential_rng(lambda_go);

  // Diagnostic: inhibition probability at the posterior mean SSRT
  real p_inh_at_mean = exp(
    log_p_inhibit(mean_ssrt,
                  mu_go, sigma_go, lambda_go,
                  mu_stop, sigma_stop, lambda_stop, 50)
  );
}
