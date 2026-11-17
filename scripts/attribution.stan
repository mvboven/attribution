////////////////////////////////////////////////////////////////////////////////
//
// Stan code for pathogen-attribution analyses
//
// Copyright (c) 2026 Michiel van Boven
// All rights reserved
//
// This software is distributed under the BSD 3-Clause License.
// Redistribution and use, with or without modification, are permitted subject
// to the conditions specified in the LICENSE file accompanying this software.
//
// This software is provided “AS IS”, without warranty of any kind, express or
// implied, including but not limited to warranties of merchantability or
// fitness for a particular purpose.
//
// Users of this code in scientific work are kindly requested to acknowledge
// Michiel van Boven
//
////////////////////////////////////////////////////////////////////////////////

functions {
  // B-spline implementation adapted from code by Milad Kharratzadeh:
  // https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order == 1) {
      for (i in 1:size(t)) {
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind + 1]);
      }
    } else {
      if (ext_knots[ind] != ext_knots[ind + order - 1]) {
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind + order - 1] - ext_knots[ind]);
      }
      if (ext_knots[ind + 1] != ext_knots[ind + order]) {
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind + 1], size(t))) /
                 (ext_knots[ind + order] - ext_knots[ind + 1]);
      }
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1) +
                 w2 .* build_b_spline(t, ext_knots, ind + 1, order - 1);
    }
    return b_spline;
  }
}

data {
  int<lower=1> num_weeks;                                            // number of weeks
  int<lower=1> num_knots;                                            // number of spline knots
  vector[num_knots] knots;                                           // knot locations
  int<lower=1> spline_degree;                                        // spline degree (typically 3)
  int<lower=1> num_types;                                            // number of age groups
  int<lower=1> num_path;                                             // number of pathogens
  int<lower=3> num_basis;                                            // num_knots + spline_degree - 1
  array[num_weeks] int<lower=0> x_data;                              // week numbers
  array[num_types, num_weeks] int<lower=0> y_data;                   // outcome data for all age groups
  array[num_path, num_weeks] real<lower=0> VIR;                      // covariates for all pathogens
  int min_lag;                                                       // minimum lag (e.g., -3 weeks)
  int max_lag;                                                       // maximum lag (e.g., 3 weeks)
  vector[num_types] prior_sigma;                                     // data-informed prior for sigma of spline weights
  vector[8] prior_season_log_effect;                                 // prior of season effect
  array[num_weeks] int<lower=1, upper=8> season_index;               // season-specific effects
}

transformed data {
  int padded_size = num_weeks + (max_lag - min_lag);                 // size of padded datasets
  int lag_offset = -min_lag;                                         // convert min_lag to a positive value
  array[padded_size] int x_data_padded;                              // padded x_data
  array[num_types, padded_size] int y_data_padded;                   // padded y_data for all age groups
  array[num_path, padded_size] real VIR_padded;                      // padded covariates for all pathogens
  matrix[num_basis, num_weeks] B;                                    // matrix of B-splines (corrected to num_weeks)
  matrix[num_weeks, num_basis] Q;                                    // orthogonal matrix from thin QR decomposition
  matrix[num_basis, num_basis] R;                                    // upper triangular matrix from thin QR decomposition
  matrix[num_basis, num_weeks] C;                                    // orthogonal design matrix
  array[spline_degree + num_knots] real ext_knots_temp;              // extended knots
  array[2 * spline_degree + num_knots] real ext_knots;               // extended knots
  int num_lags = max_lag - min_lag + 1;                              // number of lags

  // lagged design (precomputed): VIR_lag[p][w, i] = VIR at week w with lag i (aligned to pmf indexing)
  array[num_path] matrix[num_weeks, num_lags] VIR_lag;

  // padding of data to include distributed lags
  for (w in 1:lag_offset) {
    x_data_padded[w] = x_data[1];
    
  }
  for (w in (lag_offset + 1):(num_weeks - max_lag)) {
    x_data_padded[w] = x_data[w + min_lag];
  }
  for (w in (num_weeks - max_lag + 1):padded_size) {
    x_data_padded[w] = x_data[num_weeks];
  }

  for (t in 1:num_types) {
    for (w in 1:lag_offset) {
      y_data_padded[t, w] = y_data[t, 1];
    }
    for (w in (lag_offset + 1):(num_weeks - max_lag)) {
      y_data_padded[t, w] = y_data[t, w + min_lag];
    }
    for (w in (num_weeks - max_lag + 1):padded_size) {
      y_data_padded[t, w] = y_data[t, num_weeks];
    }
  }

  for (p in 1:num_path) {
    for (w in 1:lag_offset) {
      VIR_padded[p, w] = VIR[p, 1];
    }
    for (w in (lag_offset + 1):(num_weeks - max_lag)) {
      VIR_padded[p, w] = VIR[p, w + min_lag];
    }
    for (w in (num_weeks - max_lag + 1):padded_size) {
      VIR_padded[p, w] = VIR[p, num_weeks];
    }
  }

  // precompute lagged VIR design
  // note: index in padded arrays is w + (i-1), where i = 1..num_lags corresponds to lags min_lag..max_lag
  for (p in 1:num_path) {
    for (w in 1:num_weeks) {
      for (i in 1:num_lags) {
        int index = w + (i - 1);
        VIR_lag[p][w, i] = VIR_padded[p, index];
      }
    }
  }

  // make B using num_weeks
  ext_knots_temp = append_array(rep_array(knots[1], spline_degree), to_array_1d(knots));
  ext_knots = append_array(ext_knots_temp, rep_array(knots[num_knots], spline_degree));
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(x_data, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  // The recursion uses right-open intervals. At the final knot all basis
  // values would otherwise be zero. Include the final observed week in the
  // final basis function so its background is evaluated like other weeks.
  B[num_basis, num_weeks] = 1;

  // thin QR decomposition to reduce correlations between spline weights
  // notice that design matrix B is already sparse and this yields only a minor (to no) speed increase
  //Q = qr_thin_Q(B') * sqrt(num_weeks - 1);
  //R = qr_thin_R(B') / sqrt(num_weeks - 1);
  //C = R * Q';
  C = B;                                                             // keep the B-spline rather than p-spline, use LOO_IC to select knots  
}

parameters {
  matrix<lower=0>[num_types, num_path] beta;                         // regression coefficients for pathogens
  matrix<lower=min_lag, upper=max_lag>[num_types, num_path] mu;      // location parameters for pathogens
  matrix<lower=0>[num_types, num_path] alpha;                        // alpha parameters for pathogen lags for all age groups
  vector[num_types] background_intercept;                            // intercept background
  matrix[num_types, num_basis] coefficients;                         // spline coefficients (B-spline now)
  // array[num_path, num_types] sum_to_zero_vector[8] season_log_effect;// direct averaging
  array[num_path] sum_to_zero_vector[8] season_log_effect;           // no age dependencies
}

transformed parameters {
  matrix<lower=0>[num_types, num_weeks] lambda;
  array[num_types, num_path] vector[num_lags] pmf;
  // array[num_path, 8, num_types] real<lower=0> seasoneffect;
  array[num_path, 8] real<lower=0> seasoneffect;

  // lag values corresponding to the columns of VIR_lag
  vector[num_lags] lag_values =
    linspaced_vector(num_lags, min_lag, max_lag);

  // compute distributed-lag weights
  for (t in 1:num_types) {
    for (p in 1:num_path) {
      pmf[t, p] = softmax(
        -alpha[t, p] *
        abs(lag_values - mu[t, p])
      );
    }
  }

  // normalise season effects to have mean one across seasons
  /* with age
  for (p in 1:num_path) {
    for (t in 1:num_types) {
      vector[8] tmp = exp(season_log_effect[p, t]);
      // Preserve the interpretation of beta as the
      // arithmetic-average effect across seasons
      tmp /= mean(tmp);
      for (s in 1:8) {
        seasoneffect[p, s, t] = tmp[s];
      }
    }
  }
  */
  
  // no age
  for (p in 1:num_path) {
    vector[8] tmp = exp(season_log_effect[p]);
    tmp /= mean(tmp);
  for (s in 1:8) {
    seasoneffect[p, s] = tmp[s];
    }
  }

  // compute the Poisson rate using an additive model
  {
    // local variable, not dsaved
    matrix[num_types, num_weeks] linear_predictor;

    for (t in 1:num_types) {
      // center the spline coefficients so that the spline represents
      // deviations around the age-specific background intercept.
      row_vector[num_basis] centered_coefficients =
        coefficients[t, ] -
        rep_row_vector(mean(coefficients[t, ]), num_basis);

      linear_predictor[t, ] =
        rep_row_vector(background_intercept[t], num_weeks) +
        centered_coefficients * C;

      // local variable: total pathogen contribution for this age group
      row_vector[num_weeks] total_contribution =
        rep_row_vector(0, num_weeks);

      for (p in 1:num_path) {
        // distributed-lag pathogen activity for all weeks simultaneously
        vector[num_weeks] lagged_vir =
          VIR_lag[p] * pmf[t, p];

        for (w in 1:num_weeks) {
          total_contribution[w] +=
            // seasoneffect[p, season_index[w], t] * beta[t, p] * lagged_vir[w];
            seasoneffect[p, season_index[w]] * beta[t, p] * lagged_vir[w];  // no age in season effect
        }
      }

      lambda[t, ] =
        exp(linear_predictor[t, ]) +
        total_contribution +
        1e-7;
    }
  }
}

model {
  // priors for distributed-lag parameters
  to_vector(alpha) ~ gamma(10, 10);
  to_vector(mu) ~ double_exponential(0, 0.5);

  // weakly informative priors for age-specific background levels
  background_intercept ~ normal(0, 5);

  // priors for B-spline coefficients
  for (t in 1:num_types) {
    coefficients[t, ] ~ normal(0, prior_sigma[t]);
  }

  // priors for pathogen regression coefficients
  //to_vector(beta) ~ exponential(0.5); 
  //to_vector(beta) ~ lognormal(0, 1.175); // approximately 95% mass between 0.1 and 10; 1.5284 if 0.05<beta<20
  //to_vector(beta) ~ gamma(1.1168, 0.3930); //0.1-10 but with leaner tail! or gamma(2,1) -> 0.24–5.57
  to_vector(beta) ~ gamma(2,2); // main analysis
  
  // priors for season-specific effects
  /*
  for (p in 1:num_path) {
    for (t in 1:num_types) {
      season_log_effect[p, t] ~ normal(0, 0.38); // 0.15 is default; 0.38 results in 0.5-2 variation a priori
    }
  }
  */
  // no age in season effect
  for (p in 1:num_path) {
    season_log_effect[p] ~ normal(0, prior_season_log_effect); // 0.38 is default
  }
  
  // vectorised Poisson likelihood, for each age group
  for (t in 1:num_types) {
    y_data[t, ] ~ poisson(lambda[t, ]');
  }
}

generated quantities {
  // expected pathogen-specific contributions
  array[num_path, num_types, num_weeks] real lambda_contributions;

  // pointwise log likelihood for LOO-CV
  array[num_types, num_weeks] real log_lik;

  // posterior predictive total admissions
  array[num_types, num_weeks] int y_pred;

  // log likelihood and posterior predictive distribution
  for (t in 1:num_types) {
    for (w in 1:num_weeks) {
      log_lik[t, w] =
        poisson_lpmf(y_data[t, w] | lambda[t, w]);
      y_pred[t, w] =
        poisson_rng(lambda[t, w] + 1e-7);
    }
  }

  // pathogen-specific contributions
  for (p in 1:num_path) {
    for (t in 1:num_types) {
      for (w in 1:num_weeks) {
        // real se = seasoneffect[p, season_index[w], t];
        real se = seasoneffect[p, season_index[w]];
        lambda_contributions[p, t, w] =
          se *
          beta[t, p] *
          dot_product(
            pmf[t, p],
            (VIR_lag[p][w, ])'
          );
      }
    }
  }
}
