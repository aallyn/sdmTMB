#define TMB_LIB_INIT R_init_sdmTMB
#include <TMB.hpp>
#include "utils.h"
// #include <omp.h>

enum valid_family {
  gaussian_family = 0,
  binomial_family = 1,
  tweedie_family  = 2,
  poisson_family  = 3,
  Gamma_family    = 4,
  nbinom2_family  = 5,
  lognormal_family= 6,
  student_family  = 7,
  Beta_family     = 8,
  truncated_nbinom2_family  = 9,
  nbinom1_family  = 10,
  truncated_nbinom1_family  = 11,
  censored_poisson_family  = 12,
  gamma_mix_family = 13,
  lognormal_mix_family = 14,
  nbinom2_mix_family = 15
};

enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3,
  cloglog_link  = 4
};

template <class Type>
Type InverseLink(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  case cloglog_link:
    out = Type(1) - exp(-exp(eta));
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

// logit transformed inverse_linkfun without losing too much accuracy
template<class Type>
Type LogitInverseLink(Type eta, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = eta;
    break;
  case cloglog_link:
    ans = sdmTMB::logit_invcloglog(eta);
    break;
  default:
    ans = logit(InverseLink(eta, link));
  }
  return ans;
}

template <class Type>
Type Link(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = log(eta);
    break;
  case logit_link:
    out = logit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

// ------------------ Main TMB template ----------------------------------------

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Vectors of real data
  DATA_ARRAY(y_i);      // response
  DATA_STRUCT(X_ij, sdmTMB::LOM_t); // list of model matrices
  DATA_MATRIX(z_i);      // model matrix for spatial covariate effect
  DATA_MATRIX(X_rw_ik);  // model matrix for random walk covariate(s)

  DATA_STRUCT(Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_STRUCT(proj_Zs, sdmTMB::LOM_t); // [L]ist [O]f (basis function matrices) [Matrices]
  DATA_MATRIX(Xs); // smoother linear effect matrix
  DATA_MATRIX(proj_Xs); // smoother linear effect matrix

  // DATA_VECTOR_INDICATOR(keep, y_i); // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  DATA_VECTOR(weights_i); // optional weights
  DATA_VECTOR(offset_i); // optional offset
  DATA_VECTOR(proj_offset_i); // optional offset

  DATA_INTEGER(n_t);  // number of years

  // Random intercepts:
  DATA_IMATRIX(RE_indexes);
  DATA_IMATRIX(proj_RE_indexes);
  DATA_IVECTOR(nobs_RE);
  DATA_IVECTOR(ln_tau_G_index);

  DATA_SPARSE_MATRIX(A_st); // INLA 'A' projection matrix for unique stations
  DATA_IVECTOR(A_spatial_index); // Vector of stations to match up A_st output

  // Indices for factors
  DATA_FACTOR(year_i);

  DATA_INTEGER(normalize_in_r);
  DATA_INTEGER(flag);
  DATA_IVECTOR(share_range);

  // Prediction?
  DATA_INTEGER(do_predict);
  // With standard errors on the full projections?
  DATA_INTEGER(calc_se);
  // Should predictions be population (vs. individual-level) predictions?
  DATA_INTEGER(pop_pred);
  // Calculate total summed by year (e.g. biomass)?
  DATA_INTEGER(calc_index_totals);
  DATA_INTEGER(calc_cog);
  // DATA_INTEGER(calc_quadratic_range); // DELTA TODO
  DATA_VECTOR(area_i); // area per prediction grid cell for index standardization

  DATA_VECTOR(priors_b_mean);
  DATA_MATRIX(priors_b_Sigma); // beta priors matrix
  DATA_INTEGER(priors_b_n);
  DATA_IVECTOR(priors_b_index);
  DATA_VECTOR(priors); // all other priors as a vector
  DATA_IVECTOR(ar1_fields);
  DATA_IVECTOR(rw_fields);
  DATA_INTEGER(include_spatial);
  DATA_INTEGER(omit_spatial_intercept);
  DATA_INTEGER(random_walk);
  DATA_INTEGER(ar1_time);
  DATA_IVECTOR(exclude_RE); // DELTA TODO currently shared...
  DATA_INTEGER(no_spatial); // omit all spatial calculations

  DATA_VECTOR(proj_lon);
  DATA_VECTOR(proj_lat);

  // Distribution
  DATA_IVECTOR(family);
  DATA_IVECTOR(link);
  DATA_SCALAR(df);  // Student-t DF
  DATA_VECTOR(size); // binomial, via glmmTMB

  // SPDE objects from R-INLA
  DATA_STRUCT(spde_aniso, spde_aniso_t);
  DATA_STRUCT(spde, spde_t);
  PARAMETER_ARRAY(ln_H_input);
  DATA_INTEGER(anisotropy);

  // Barrier
  DATA_INTEGER(barrier);
  DATA_STRUCT(spde_barrier, sdmTMB::spde_barrier_t);
  DATA_VECTOR(barrier_scaling); // scaling of range

  // Projections
  DATA_SPARSE_MATRIX(proj_mesh);
  DATA_STRUCT(proj_X_ij, sdmTMB::LOM_t);
  DATA_MATRIX(proj_X_rw_ik);
  DATA_FACTOR(proj_year);
  DATA_MATRIX(proj_z_i);
  DATA_IVECTOR(proj_spatial_index);

  DATA_IVECTOR(spatial_only);
  DATA_INTEGER(spatial_covariate);

  DATA_VECTOR(X_threshold);
  DATA_VECTOR(proj_X_threshold);
  DATA_INTEGER(threshold_func);
  // optional model for nonstationary st variance
  DATA_INTEGER(est_epsilon_model);
  DATA_INTEGER(est_epsilon_slope);
  DATA_INTEGER(est_epsilon_re);
  DATA_VECTOR(epsilon_predictor);

  // optional stuff for penalized regression splines
  DATA_INTEGER(has_smooths);  // whether or not smooths are included
  DATA_IVECTOR(b_smooth_start);

  DATA_IVECTOR(sim_re); // sim random effects? 0,1; order: omega, epsilon, zeta, IID, RW, smoothers
  DATA_IVECTOR(simulate_t); // sim this specific time step? (used for forecasting)

  DATA_VECTOR(lwr); // lower bound for censpois on counts
  DATA_VECTOR(upr); // upper bound for censpois on counts
  DATA_INTEGER(poisson_link_delta); // logical

  DATA_INTEGER(stan_flag); // logical whether to pass the model to Stan
  // ------------------ Parameters ---------------------------------------------

  // Parameters
  // Fixed effects
  PARAMETER_VECTOR(b_j);  // fixed effect parameters
  PARAMETER_VECTOR(b_j2);  // fixed effect parameters delta2 part
  PARAMETER_ARRAY(bs); // smoother linear effects
  PARAMETER_VECTOR(ln_tau_O);    // spatial process
  PARAMETER_ARRAY(ln_tau_Z);    // optional spatially varying covariate process
  PARAMETER_VECTOR(ln_tau_E);    // spatio-temporal process
  PARAMETER_ARRAY(ln_kappa);    // Matern parameter

  PARAMETER(thetaf);           // tweedie only
  PARAMETER(logit_p_mix);           // ECE / positive mixture only
  PARAMETER(log_ratio_mix);           // ECE / positive mixture only

  PARAMETER_VECTOR(ln_phi);           // sigma / dispersion / etc.
  PARAMETER_ARRAY(ln_tau_V);  // random walk sigma
  PARAMETER_ARRAY(rho_time_unscaled); // (k, m) dimension ar1 time correlation rho -Inf to Inf
  PARAMETER_VECTOR(ar1_phi);          // AR1 fields correlation
  PARAMETER_ARRAY(ln_tau_G);  // random intercept sigmas
  PARAMETER_ARRAY(RE);        // random intercept deviations
  // Random effects
  PARAMETER_ARRAY(b_rw_t);  // random walk effects
  PARAMETER_ARRAY(omega_s);    // spatial effects; n_s length
  PARAMETER_ARRAY(zeta_s);    // spatial effects on covariate; n_s length, n_z cols, n_m
  PARAMETER_ARRAY(epsilon_st);  // spatio-temporal effects; n_s by n_t by n_m array
  PARAMETER_ARRAY(b_threshold);  // coefficients for threshold relationship (3) // DELTA TODO
  PARAMETER_VECTOR(b_epsilon); // slope coefficient for log-linear model on epsilon
  PARAMETER_VECTOR(ln_epsilon_re_sigma);
  PARAMETER_ARRAY(epsilon_re);
  PARAMETER_ARRAY(b_smooth);  // P-spline smooth parameters
  PARAMETER_ARRAY(ln_smooth_sigma);  // variances of spline REs if included

  // Joint negative log-likelihood
  Type jnll = 0.;

  // ------------------ End of parameters --------------------------------------

  // DELTA DONE
  int n_i = y_i.rows();   // number of observations
  int n_m = y_i.cols();   // number of models (delta)
  int n_RE = RE_indexes.cols();  // number of random effect intercepts

  // DELTA TODO
  // ------------------ Derived variables -------------------------------------------------
  vector<Type> s_slope(n_m), s_cut(n_m), s50(n_m), s95(n_m), s_max(n_m);
  // these are for linear model
  for (int m = 0; m < n_m; m++) {
    s_slope(m) = b_threshold(0,m);
    s_cut(m) = b_threshold(1,m);
    if (threshold_func == 2) {
      s50(m) = b_threshold(0,m); // threshold at which function is 50% of max
      s95(m) = b_threshold(0,m) + exp(b_threshold(1,m)); // threshold at which function is 95% of max
      s_max(m) = b_threshold(2,m);
    }
  }
  // DELTA DONE
  vector<Type> rho(n_m);
  for (int m = 0; m < n_m; m++) rho(m) = sdmTMB::minus_one_to_one(ar1_phi(m));
  vector<Type> phi = exp(ln_phi);

  // ------------------ Geospatial ---------------------------------------------

  // DELTA DONE
  // Matern:
  array<Type> range(2,n_m);
  for (int m = 0; m < n_m; m++) {
    for (int r = 0; r < 2; r++) {
      range(r,m) = sqrt(Type(8.)) / exp(ln_kappa(r,m));
    }
  }

  // DELTA DONE
  vector<Type> sigma_O(n_m);
  int n_z = ln_tau_Z.rows();
  array<Type> sigma_Z(n_z, n_m);
  if (include_spatial) {
    for (int m = 0; m < n_m; m++) {
      sigma_O(m) = sdmTMB::calc_rf_sigma(ln_tau_O(m), ln_kappa(0,m));
      for (int z = 0; z < n_z; z++) {
        sigma_Z(z,m) = sdmTMB::calc_rf_sigma(ln_tau_Z(z,m), ln_kappa(0,m));
      }
    }
    vector<Type> log_sigma_O = log(sigma_O);
    ADREPORT(log_sigma_O);
    REPORT(sigma_O);
    ADREPORT(sigma_O);
    array<Type> log_sigma_Z(n_z,n_m); // for SE
    for (int z = 0; z < n_z; z++)
      for (int m = 0; m < n_m; m++)
        log_sigma_Z(z,m) = log(sigma_Z(z,m));
    ADREPORT(log_sigma_Z);
    REPORT(sigma_Z);
    ADREPORT(sigma_Z);
  }

  // TODO can we not always run this for speed?
  //vector<Type> sigma_E(n_m);
  //for (int m = 0; m < n_m; m++) {
  //  sigma_E(m) = sdmTMB::calc_rf_sigma(ln_tau_E(m), ln_kappa(1,m));
  //}

  // optional non-stationary model on epsilon
  array<Type> sigma_E(n_t, n_m);
  array<Type> ln_tau_E_vec(n_t, n_m);
  if (!est_epsilon_model) { // constant model
    for (int m = 0; m < n_m; m++) {
      // do calculation once,
      sigma_E(0,m) = sdmTMB::calc_rf_sigma(ln_tau_E(m), ln_kappa(1,m));
      ln_tau_E_vec(0,m) = ln_tau_E(m);
      for (int i = 1; i < n_t; i++) {
        sigma_E(i,m) = sigma_E(0,m);
        ln_tau_E_vec(i,m) = ln_tau_E_vec(0,m);
      }
    }
  }
  if (est_epsilon_model) { // loglinear model
    // epsilon_intcpt is the intercept parameter, derived from ln_tau_E.
    // For models with time as covariate, this is interpreted as sigma when covariate = 0.
    for (int m = 0; m < n_m; m++) {

      Type epsilon_intcpt = sdmTMB::calc_rf_sigma(ln_tau_E(m), ln_kappa(1,m));
      Type log_epsilon_intcpt = log(epsilon_intcpt);
      Type log_epsilon_temp = 0.0;
      Type epsilon_cnst = - log(Type(4.0) * M_PI) / Type(2.0) - ln_kappa(1,m);
      if (est_epsilon_re) {
        Type epsilon_re_sigma = exp(ln_epsilon_re_sigma(m));
        for (int i = 0; i < n_t; i++) {
          jnll -= dnorm(epsilon_re(i,m), Type(0), Type(epsilon_re_sigma), true);
        }
      }

      for(int i = 0; i < n_t; i++) {
        log_epsilon_temp = log_epsilon_intcpt;
        if (est_epsilon_slope) log_epsilon_temp += b_epsilon(m) * epsilon_predictor(i);
        if (est_epsilon_re) log_epsilon_temp += epsilon_re(i,m);
        sigma_E(i,m) = exp(log_epsilon_temp); // log-linear model
        ln_tau_E_vec(i,m) = -log_epsilon_temp + epsilon_cnst;
      }
    }
  }

  Eigen::SparseMatrix<Type> Q_s; // Precision matrix
  Eigen::SparseMatrix<Type> Q_st; // Precision matrix
  Eigen::SparseMatrix<Type> Q_s2; // Precision matrix
  Eigen::SparseMatrix<Type> Q_st2; // Precision matrix

  if (barrier) {
    Q_s = Q_spde(spde_barrier, exp(ln_kappa(0,0)), barrier_scaling);
    if (n_m > 1) Q_s2 = Q_spde(spde_barrier, exp(ln_kappa(0,1)), barrier_scaling);
    if (!share_range(0)) Q_st = Q_spde(spde_barrier, exp(ln_kappa(1,0)), barrier_scaling);
    if (!share_range(1) && n_m > 1) Q_st2 = Q_spde(spde_barrier, exp(ln_kappa(1,1)), barrier_scaling);
  } else {
    if (anisotropy) {
      matrix<Type> H = sdmTMB::MakeH(vector<Type>(ln_H_input.col(0)));
      Q_s = R_inla::Q_spde(spde_aniso, exp(ln_kappa(0,0)), H);
      if (!share_range(0)) Q_st = R_inla::Q_spde(spde_aniso, exp(ln_kappa(1,0)), H);
      REPORT(H);
      if (n_m > 1) {
        matrix<Type> H2 = sdmTMB::MakeH(vector<Type>(ln_H_input.col(1)));
        Q_s2 = R_inla::Q_spde(spde_aniso, exp(ln_kappa(0,1)), H2);
        if (!share_range(1)) Q_st2 = R_inla::Q_spde(spde_aniso, exp(ln_kappa(1,1)), H2);
        REPORT(H2);
      }
    }
    if (!anisotropy) {
      Q_s = R_inla::Q_spde(spde, exp(ln_kappa(0,0)));
      if (!share_range(0)) Q_st = R_inla::Q_spde(spde, exp(ln_kappa(1,0)));
      if (n_m > 1) Q_s2 = R_inla::Q_spde(spde, exp(ln_kappa(0,1)));
      if (!share_range(1) && n_m > 1) Q_st2 = R_inla::Q_spde(spde, exp(ln_kappa(1,1)));
    }
  }
  if (share_range(0)) Q_st = Q_s;
  if (share_range(1)) Q_st2 = Q_s2;

  bool s = true;
  if (normalize_in_r) s = false;

  // Spatial (intercept) random effects:
  for (int m = 0; m < n_m; m++) {
    Eigen::SparseMatrix<Type> Q_temp; // Precision matrix
    if (include_spatial) {
      if (m == 0) {
        Q_temp = Q_s;
      } else {
        Q_temp = Q_s2;
      }
      if (!omit_spatial_intercept) {
        PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1. / exp(ln_tau_O(m)))(omega_s.col(m));
        if (sim_re(0)) {
          vector<Type> omega_s_tmp(omega_s.rows());
          SIMULATE {
            GMRF(Q_temp, s).simulate(omega_s_tmp);
            omega_s.col(m) = omega_s_tmp / exp(ln_tau_O(m));
          }
        }
      }
      if (spatial_covariate) {
        for (int z = 0; z < n_z; z++) {
          PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1. / exp(ln_tau_Z(z,m)))(zeta_s.col(m).col(z));
          if (sim_re(3)) {
            vector<Type> zeta_s_tmp(zeta_s.col(m).rows());
            SIMULATE {
              GMRF(Q_s, s).simulate(zeta_s_tmp);
              zeta_s.col(m).col(z) = zeta_s_tmp / exp(ln_tau_Z(z,m));
            }
          }
        }
      }
    }
  }

  // Spatiotemporal random effects:
  for (int m = 0; m < n_m; m++) {
    Eigen::SparseMatrix<Type> Q_temp; // Precision matrix
    if (m == 0) {
      Q_temp = Q_st;
    } else {
      Q_temp = Q_st2;
    }
    if (!spatial_only(m)) {
      if (!ar1_fields(m) && !rw_fields(m)) {
        for (int t = 0; t < n_t; t++)
          PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1. / exp(ln_tau_E_vec(t,m)))(epsilon_st.col(m).col(t));
        if (sim_re(1)) {
          for (int t = 0; t < n_t; t++) {
            if (simulate_t(t)) {
              vector<Type> epsilon_st_tmp(epsilon_st.col(m).rows());
              SIMULATE {GMRF(Q_temp, s).simulate(epsilon_st_tmp);
                epsilon_st.col(m).col(t) = epsilon_st_tmp / exp(ln_tau_E_vec(t,m));}
            }
          }
        }
      } else {
        if (ar1_fields(m)) { // not using separable(ar1()) so we can simulate by time step
          // PARALLEL_REGION jnll += SCALE(SEPARABLE(AR1(rho(m)), GMRF(Q_temp, s)), 1./exp(ln_tau_E(m)))(epsilon_st.col(m));
          // Split out by year so we can turn on/off simulation by year and model covariates of ln_tau_E:
          PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1. / exp(ln_tau_E_vec(0,m)))(epsilon_st.col(m).col(0));
          for (int t = 1; t < n_t; t++) {
            PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1. /exp(ln_tau_E_vec(t,m)))((epsilon_st.col(m).col(t) -
              rho(m) * epsilon_st.col(m).col(t - 1))/sqrt(1. - rho(m) * rho(m)));
          }
          Type n_cols = epsilon_st.col(m).cols();
          Type n_rows = epsilon_st.col(m).rows();
          // Penalty to match TMB AR1_t() implementation:
          PARALLEL_REGION jnll += Type((n_cols - 1.) * n_rows) * log(sqrt(1. - rho(m) * rho(m)));
          if (sim_re(1)) {
            // array<Type> epsilon_st_tmp(epsilon_st.col(m).rows(),n_t);
            // SIMULATE {SEPARABLE(AR1(rho(m)), GMRF(Q_temp, s)).simulate(epsilon_st_tmp);
            //   epsilon_st.col(m) = epsilon_st_tmp / exp(ln_tau_E(m));}
            for (int t = 0; t < n_t; t++) {
              if (simulate_t(t)) {
                vector<Type> epsilon_st_tmp(epsilon_st.col(m).rows());
                SIMULATE {
                  GMRF(Q_temp, s).simulate(epsilon_st_tmp);
                  epsilon_st_tmp *= 1./exp(ln_tau_E_vec(t,m));
                  // https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html
                  Type ar1_scaler = sqrt(1. - rho(m) * rho(m));
                  if (t == 0) {
                    epsilon_st.col(m).col(0) = epsilon_st_tmp; // no scaling of first step
                  } else {
                    epsilon_st.col(m).col(t) = rho(m) * epsilon_st.col(m).col(t-1) + epsilon_st_tmp * ar1_scaler;
                  }
                }
              }
            }
          }
          ADREPORT(rho);
        } else if (rw_fields(m)) {
          PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s), 1./exp(ln_tau_E_vec(0,m)))(epsilon_st.col(m).col(0));
          for (int t = 1; t < n_t; t++) {
            PARALLEL_REGION jnll += SCALE(GMRF(Q_temp, s),
                1./exp(ln_tau_E_vec(t,m)))(epsilon_st.col(m).col(t) - epsilon_st.col(m).col(t - 1));
          }
          if (sim_re(1)) {
            for (int t = 0; t < n_t; t++) {
              if (simulate_t(t)) {
                vector<Type> epsilon_st_tmp(epsilon_st.col(m).rows());
                SIMULATE {
                  GMRF(Q_temp, s).simulate(epsilon_st_tmp);
                  epsilon_st_tmp *= 1./exp(ln_tau_E_vec(t,m));
                  if (t == 0) {
                    epsilon_st.col(m).col(0) = epsilon_st_tmp;
                  } else {
                    epsilon_st.col(m).col(t) = epsilon_st.col(m).col(t-1) + epsilon_st_tmp;
                  }
                }
              }
            }
          }
        } else {
          error("Field type not implemented.");
        }
      }
    }
  }
  if (flag == 0) return jnll;

  // ------------------ Probability of random effects --------------------------

  // IID random intercepts:
  for (int m = 0; m < n_m; m++) {
    for (int g = 0; g < RE.rows(); g++) {
      PARALLEL_REGION jnll -= dnorm(RE(g,m), Type(0), exp(ln_tau_G(ln_tau_G_index(g), m)), true);
      if (sim_re(3)) SIMULATE{RE(g,m) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g),m)));
      }
    }
  }

  // Time-varying effects (dynamic regression):
  if (random_walk == 1 || ar1_time || random_walk == 2) {
    array<Type> rho_time(X_rw_ik.cols(), n_m);
    rho_time.setZero();
    for (int m = 0; m < n_m; m++) {
      for (int k = 0; k < X_rw_ik.cols(); k++) {
        if (random_walk == 1) { // type = 'rw'
          // flat prior on the initial value... then:
          for (int t = 1; t < n_t; t++) {
            PARALLEL_REGION jnll -=
              dnorm(b_rw_t(t, k, m), b_rw_t(t - 1, k, m), exp(ln_tau_V(k,m)), true);
            if (sim_re(4) && simulate_t(t))
              SIMULATE{b_rw_t(t, k, m) = rnorm(b_rw_t(t - 1, k, m), exp(ln_tau_V(k,m)));}
          }
        } else if (random_walk == 2) { // type = 'rw0'
          // N(0, SD) prior on the initial value... then:
          for (int t = 0; t < n_t; t++) {
            if (t == 0) {
              PARALLEL_REGION jnll -=
                dnorm(b_rw_t(t, k, m), Type(0.), exp(ln_tau_V(k,m)), true);
              if (sim_re(4) && simulate_t(t))
                SIMULATE{b_rw_t(t, k, m) = rnorm(Type(0.), exp(ln_tau_V(k,m)));}
            } else {
              PARALLEL_REGION jnll -=
                dnorm(b_rw_t(t, k, m), b_rw_t(t - 1, k, m), exp(ln_tau_V(k,m)), true);
              if (sim_re(4) && simulate_t(t))
                SIMULATE{b_rw_t(t, k, m) = rnorm(b_rw_t(t - 1, k, m), exp(ln_tau_V(k,m)));}
            }
          }
        } else if (ar1_time) { // type = 'ar1'
          rho_time(k, m) = sdmTMB::minus_one_to_one(rho_time_unscaled(k, m));
          jnll += SCALE(AR1(rho_time(k, m)), exp(ln_tau_V(k, m)))(vector<Type>(b_rw_t.col(m).col(k)));
          if (sim_re(4)) {
            // https://kaskr.github.io/adcomp/classdensity_1_1AR1__t.html
            vector<Type> tmp(n_t);
            Type ar1_sigma = sqrt(Type(1) - rho_time(k, m) * rho_time(k, m));
            Type x0 = rnorm(Type(0), exp(ln_tau_V(k, m)));
            tmp(0) = rho_time(k, m) * x0 + ar1_sigma * rnorm(Type(0), exp(ln_tau_V(k, m)));
            for (int t = 1; t < n_t; t++)
              tmp(t) = rho_time(k, m) * tmp(t-1) + ar1_sigma * rnorm(Type(0), exp(ln_tau_V(k, m)));
            for(int t = 0; t < n_t; t++)
              if (simulate_t(t)) SIMULATE{b_rw_t(t, k, m) = tmp(t);}
          }
        } else {
          error("Time-varying type not found.");
        }
      }
    }
  }
  // ------------------ INLA projections ---------------------------------------

  // Here we are projecting the spatiotemporal and spatial random effects to the
  // locations of the data using the INLA 'A' matrices.
  array<Type> omega_s_A(n_i, n_m);
  array<Type> zeta_s_A(n_i, n_z, n_m);
  array<Type> epsilon_st_A(n_i, n_t, n_m);
  array<Type> epsilon_st_A_vec(n_i, n_m);
  omega_s_A.setZero();
  zeta_s_A.setZero();
  epsilon_st_A.setZero();
  epsilon_st_A_vec.setZero();

  if (!no_spatial) {
    for (int m = 0; m < n_m; m++) {
      for (int t = 0; t < n_t; t++)
        if (!spatial_only(m)) epsilon_st_A.col(m).col(t) =
          A_st * vector<Type>(epsilon_st.col(m).col(t));
      if (!omit_spatial_intercept) omega_s_A.col(m) = A_st * vector<Type>(omega_s.col(m));
      for (int z = 0; z < n_z; z++)
        zeta_s_A.col(m).col(z) = A_st * vector<Type>(zeta_s.col(m).col(z));
    }
  }

  // ------------------ Linear predictor ---------------------------------------

  array<Type> eta_fixed_i(n_i, n_m);
  for (int m = 0; m < n_m; m++) {
    if (m == 0) eta_fixed_i.col(m) = X_ij(m) * b_j;
    if (m == 1) eta_fixed_i.col(m) = X_ij(m) * b_j2;
  }

  // FIXME delta must be same in 2 components:
  // p-splines/smoothers
  array<Type> eta_smooth_i(n_i, n_m);
  eta_smooth_i.setZero();
  if (has_smooths) {
    for (int m = 0; m < n_m; m++) {
      for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
        array<Type> beta_s(Zs(s).cols(),n_m);
        beta_s.setZero();
        for (int j = 0; j < beta_s.rows(); j++) {
          beta_s(j,m) = b_smooth(b_smooth_start(s) + j,m);
          PARALLEL_REGION jnll -= dnorm(beta_s(j,m), Type(0), exp(ln_smooth_sigma(s,m)), true);
          if (sim_re(5)) SIMULATE{beta_s(j,m) = rnorm(Type(0), exp(ln_smooth_sigma(s,m)));}
        }
        eta_smooth_i.col(m) += Zs(s) * vector<Type>(beta_s.col(m));
      }
      eta_smooth_i.col(m) += Xs * vector<Type>(bs.col(m));
    }
  }

  // add threshold effect if specified
  if (threshold_func > 0) {
    if (threshold_func == 1) {
      // linear
      for (int m = 0; m < n_m; m++) {
        for (int i = 0; i < n_i; i++) {
          eta_fixed_i(i,m) += sdmTMB::linear_threshold(X_threshold(i), s_slope(m), s_cut(m));
        }
      }
    } else {
      // logistic
      for (int m = 0; m < n_m; m++) {
        for (int i = 0; i < n_i; i++) {
          eta_fixed_i(i,m) += sdmTMB::logistic_threshold(X_threshold(i), s50(m), s95(m), s_max(m));
        }
      }
    }
  }

  matrix<Type> mu_i(n_i,n_m), eta_i(n_i,n_m), eta_rw_i(n_i,n_m), eta_iid_re_i(n_i,n_m);
  eta_rw_i.setZero();
  eta_iid_re_i.setZero();
  mu_i.setZero();
  eta_i.setZero();

  // combine parts:
  for (int m = 0; m < n_m; m++) {
    for (int i = 0; i < n_i; i++) {
      eta_i(i,m) = eta_fixed_i(i,m) + eta_smooth_i(i,m);
      if ((n_m == 2 && m == 1) || n_m == 1) {
        if (!poisson_link_delta) eta_i(i,m) += offset_i(i);
      }
      if (random_walk == 1 || ar1_time || random_walk == 2) {
        for (int k = 0; k < X_rw_ik.cols(); k++) {
          eta_rw_i(i,m) += X_rw_ik(i, k) * b_rw_t(year_i(i), k, m); // record it
          eta_i(i,m) += eta_rw_i(i,m);
        }
      }

      // Spatially varying effects:
      if (include_spatial) {
        if (!omit_spatial_intercept)
          eta_i(i,m) += omega_s_A(i,m);  // spatial omega
        if (spatial_covariate)
          for (int z = 0; z < n_z; z++)
            eta_i(i,m) += zeta_s_A(i,z,m) * z_i(i,z); // spatially varying covariate DELTA
      }
      if (!no_spatial) epsilon_st_A_vec(i,m) = epsilon_st_A(A_spatial_index(i), year_i(i),m); // record it
      eta_i(i,m) += epsilon_st_A_vec(i,m); // spatiotemporal

      // IID random intercepts:
      int temp = 0;
      for (int k = 0; k < n_RE; k++) {
        if (k == 0) eta_iid_re_i(i,m) += RE(RE_indexes(i, k),m); // record it
        if (k > 0) {
          temp += nobs_RE(k - 1);
          eta_iid_re_i(i,m) += RE(RE_indexes(i, k) + temp,m); // record it
        }
      }
      eta_i(i,m) += eta_iid_re_i(i,m);

      // bool poisson_link_delta = false;
      // if (n_m > 1) if (family(0) == 1 && family(1) == 4 && link(0) == 1 && link(1) == 1)
      //     poisson_link_delta = true;

      if (family(m) == 1 && !poisson_link_delta) { // regular binomial
        // binomial(link = "logit"); don't touch (using robust density function in logit space)
        mu_i(i,m) = LogitInverseLink(eta_i(i,m), link(m));
      } else if (poisson_link_delta) { // clogog, but put in logit space for robust density function:
        Type n = exp(eta_i(i,0));
        // Type p = Type(1) - exp(-exp(offset_i(i)) * n);
        Type p = Type(1) - exp(-exp(offset_i(i) + eta_i(i,0)));
        if (m == 0) mu_i(i,0) = logit(p);
        // if (m == 1) mu_i(i,1) = (n/p) * exp(eta_i(i,1));
        if (m == 1) mu_i(i,1) = exp(log(n/p) + eta_i(i,1));
      } else {
        mu_i(i,m) = InverseLink(eta_i(i,m), link(m));
      }
    }
  }

  // ------------------ Probability of data given random effects ---------------

  // from glmmTMB:
  // close to zero: use for count data (cf binomial()$initialize)
#define zt_lik_nearzero(x,loglik_exp) ((x < Type(0.001)) ? -INFINITY : loglik_exp)

  Type s1, s2, s3, lognzprob, tmp_ll, ll_1, ll_2, p_mix, mix_ratio, ll_max;
  REPORT(phi);
  for (int m = 0; m < n_m; m++) PARALLEL_REGION {
    for (int i = 0; i < n_i; i++) {
      if (!sdmTMB::isNA(y_i(i,m))) {
        switch (family(m)) {
          case gaussian_family: {
            tmp_ll = dnorm(y_i(i,m), mu_i(i,m), phi(m), true);
            SIMULATE{y_i(i,m) = rnorm(mu_i(i,m), phi(m));}
            break;
          }
          case tweedie_family: {
            s1 = invlogit(thetaf) + Type(1.0);
            if (!sdmTMB::isNA(priors(12))) {
              jnll -= dnorm(s1, priors(12), priors(13), true);
              // derivative: https://www.wolframalpha.com/input?i=e%5Ex%2F%281%2Be%5Ex%29+%2B+1
              if (stan_flag) jnll -= thetaf - 2 * log(1 + exp(thetaf)); // Jacobian adjustment
            }
            tmp_ll = dtweedie(y_i(i,m), mu_i(i,m), phi(m), s1, true);
            SIMULATE{y_i(i,m) = rtweedie(mu_i(i,m), phi(m), s1);}
            break;
          }
          case binomial_family: {
            tmp_ll = dbinom_robust(y_i(i,m), size(i), mu_i(i,m), true);
            // SIMULATE{y_i(i,m) = rbinom(size(i), InverseLink(mu_i(i,m), link(m)));}
            SIMULATE{y_i(i,m) = rbinom(size(i), invlogit(mu_i(i,m)));} // hardcoded invlogit b/c mu_i in logit space
            break;
          }
          case poisson_family: {
            tmp_ll = dpois(y_i(i,m), mu_i(i,m), true);
            SIMULATE{y_i(i,m) = rpois(mu_i(i,m));}
            break;
          }
          case censored_poisson_family: {
            tmp_ll = sdmTMB::dcenspois(y_i(i,m), mu_i(i,m), lwr(i), upr(i), true);
            SIMULATE{y_i(i,m) = rpois(mu_i(i,m));}
            break;
          }
          case Gamma_family: {
            s1 = exp(ln_phi(m));        // shape
            s2 = mu_i(i,m) / s1;        // scale
            tmp_ll = dgamma(y_i(i,m), s1, s2, true);
            SIMULATE{y_i(i,m) = rgamma(s1, s2);}
            // s1 = Type(1) / (pow(phi, Type(2)));  // s1=shape, ln_phi=CV,shape=1/CV^2
            // tmp_ll = dgamma(y_i(i,m), s1, mu_i(i,m) / s1, true);
            break;
          }
          case nbinom2_family: {
            s1 = log(mu_i(i,m)); // log(mu_i)
            s2 = 2. * s1 - ln_phi(m); // log(var - mu)
            tmp_ll = dnbinom_robust(y_i(i,m), s1, s2, true);
            SIMULATE { // from glmmTMB
              s1 = mu_i(i,m);
              s2 = mu_i(i,m) * (Type(1) + mu_i(i,m) / phi(m));
              y_i(i,m) = rnbinom2(s1, s2);
            }
            break;
          }
          case truncated_nbinom2_family: {
            s1 = log(mu_i(i,m)); // log(mu_i)
            s2 = 2. * s1 - ln_phi(m); // log(var - mu)
            tmp_ll = dnbinom_robust(y_i(i,m), s1, s2, true);
            s3 = logspace_add(Type(0), s1 - ln_phi(m));
            lognzprob = logspace_sub(Type(0), -phi(m) * s3);
            tmp_ll -= lognzprob;
            tmp_ll = zt_lik_nearzero(y_i(i,m), tmp_ll); // from glmmTMB
            SIMULATE{y_i(i,m) = sdmTMB::rtruncated_nbinom(asDouble(phi(m)), 0, asDouble(mu_i(i,m)));}
            break;
          }
          case nbinom1_family: {
            s1 = log(mu_i(i,m));
            s2 = s1 + ln_phi(m);
            tmp_ll = dnbinom_robust(y_i(i,m), s1, s2, true);
            SIMULATE {y_i(i,m) = rnbinom2(mu_i(i,m), mu_i(i,m) * (Type(1) + phi(m)));}
            break;
          }
          case truncated_nbinom1_family: {
            s1 = log(mu_i(i,m));
            s2 = s1 + ln_phi(m);
            tmp_ll = dnbinom_robust(y_i(i,m), s1, s2, true);
            s3 = logspace_add(Type(0), ln_phi(m));
            lognzprob = logspace_sub(Type(0), -mu_i(i,m) / phi(m) * s3); // 1-prob(0)
            tmp_ll -= lognzprob;
            tmp_ll = zt_lik_nearzero(y_i(i,m), tmp_ll);
            SIMULATE{y_i(i,m) = sdmTMB::rtruncated_nbinom(asDouble(mu_i(i,m)/phi(m)), 0, asDouble(mu_i(i,m)));}
            break;
          }
          case lognormal_family: {
            tmp_ll = sdmTMB::dlnorm(y_i(i,m), log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m), true);
            SIMULATE{y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));}
            break;
          }
          case student_family: {
            tmp_ll = sdmTMB::dstudent(y_i(i,m), mu_i(i,m), exp(ln_phi(m)), df, true);
            SIMULATE{y_i(i,m) = mu_i(i,m) + phi(m) * rt(df);}
            break;
          }
          case Beta_family: { // Ferrari and Cribari-Neto 2004; betareg package
            s1 = mu_i(i,m) * phi(m);
            s2 = (Type(1) - mu_i(i,m)) * phi(m);
            tmp_ll = dbeta(y_i(i,m), s1, s2, true);
            SIMULATE{y_i(i,m) = rbeta(s1, s2);}
            break;
          }
          case gamma_mix_family: {
            p_mix = invlogit(logit_p_mix); // probability of larger event
            s1 = exp(ln_phi(m));        // shape
            s2 = mu_i(i,m) / s1;        // scale
            ll_1 = log(Type(1.0 - p_mix)) + dgamma(y_i(i,m), s1, s2, true);

            mix_ratio = exp(log_ratio_mix) + 1.0; // ratio of large:small values, constrained > 1.0
            s3 = s2 * mix_ratio; // mean of large component = mean of smaller * ratio
            ll_2 = log(p_mix) + dgamma(y_i(i,m), s1, s3, true);

            ll_max = ll_1; // determine larger of ll_1 and ll_2
            if(ll_2 > ll_1) ll_max = ll_2;

            tmp_ll = ll_max + log(exp(ll_1 - ll_max) + exp(ll_2 - ll_max)); // log sum exp function

            SIMULATE{
              if(rbinom(Type(1),p_mix) == 0) {
                y_i(i,m) = rgamma(s1, s2);
              } else {
                y_i(i,m) = rgamma(s1, s3);
              }
            }
            // s1 = Type(1) / (pow(phi, Type(2)));  // s1=shape, ln_phi=CV,shape=1/CV^2
            // tmp_ll = dgamma(y_i(i,m), s1, mu_i(i,m) / s1, true);
            break;
          }
        case lognormal_mix_family: {
          p_mix = invlogit(logit_p_mix); // probability of larger event

          ll_1 = log(Type(1.0 - p_mix)) + sdmTMB::dlnorm(y_i(i,m), log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m), true);
          mix_ratio = exp(log_ratio_mix) + 1.0; // ratio of large:small values, constrained > 1.0
          ll_2 = log(p_mix) + sdmTMB::dlnorm(y_i(i,m), log(mix_ratio) + log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m), true);

          ll_max = ll_1; // determine larger of ll_1 and ll_2
          if(ll_2 > ll_1) ll_max = ll_2;

          tmp_ll = ll_max + log(exp(ll_1 - ll_max) + exp(ll_2 - ll_max)); // log sum exp function

          SIMULATE{
            if(rbinom(Type(1),p_mix) == 0) {
              y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));;
            } else {
              y_i(i,m) = exp(rnorm(log(mix_ratio) + log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));;
            }
          }
          // s1 = Type(1) / (pow(phi, Type(2)));  // s1=shape, ln_phi=CV,shape=1/CV^2
          // tmp_ll = dgamma(y_i(i,m), s1, mu_i(i,m) / s1, true);
          break;
        }
        case nbinom2_mix_family: {
          p_mix = invlogit(logit_p_mix); // probability of larger event
          Type mu_i_large = mu_i(i,m) * mix_ratio;  // mean of large component = mean of smaller * ratio
          s1 = log(mu_i(i,m)); // log(mu_i)
          s2 = 2. * s1 - ln_phi(m); // log(var - mu)
          Type s1_large = log(mu_i_large);
          Type s2_large = 2. * s1_large - ln_phi(m);
          mix_ratio = exp(log_ratio_mix) + 1.0; // ratio of large:small values, constrained > 1.0
          ll_1 = log(Type(1.0 - p_mix)) + dnbinom_robust(y_i(i,m), s1, s2, true);
          ll_2 = log(p_mix) + dnbinom_robust(y_i(i,m), s1_large, s2_large, true);
          tmp_ll = ll_max + log(exp(ll_1 - ll_max) + exp(ll_2 - ll_max)); // log sum exp function
          SIMULATE{
            if (rbinom(Type(1), p_mix) == 0) {
              y_i(i,m) = rnbinom2(s1, s2);
            } else {
              y_i(i,m) = rnbinom2(s1_large, s2_large);
            }
          }
          break;
        }
        default:
          error("Family not implemented.");
        }
        tmp_ll *= weights_i(i);
        jnll -= tmp_ll; // * keep
      }
    }
  }

  // ------------------ Priors -------------------------------------------------

  // Construct special object for MVN distribution; always has mean 0.
  MVNORM_t<Type> neg_log_dmvnorm(priors_b_Sigma);
  // Apply nll on residual. Note that other univariate densities are positive
  // log-likelihoods but the dmvnorm is negative.
  // We're accumulating the neg LL, which is why this is a + sign.

  // DELTA TODO split b_j priors!?
  // otherwise done
  for (int m = 0; m < n_m; m++) {
    if (priors_b_n > 0) {
      vector<Type> b_j_subset(priors_b_n),b_mean_subset(priors_b_n);
      for(int j = 0; j < priors_b_n; j++) {
        if (m == 0) b_j_subset(j) = b_j(priors_b_index(j));
        if (m == 1) b_j_subset(j) = b_j2(priors_b_index(j));
        b_mean_subset(j) = priors_b_mean(j);
      }
      jnll += neg_log_dmvnorm(b_j_subset - b_mean_subset);
    }

    if (!sdmTMB::isNA(priors(0)) && !sdmTMB::isNA(priors(1)) &&
        !sdmTMB::isNA(priors(2)) && !sdmTMB::isNA(priors(3))) {
      // std::cout << "Using spatial PC prior" << "\n";
      jnll -= sdmTMB::pc_prior_matern(
          ln_tau_O(m), ln_kappa(0,m),
          priors(0), priors(1), priors(2), priors(3),
          true, /* log */
          false, /* share range */
          stan_flag);
    }
    if (!sdmTMB::isNA(priors(4)) && !sdmTMB::isNA(priors(5)) &&
        !sdmTMB::isNA(priors(6)) && !sdmTMB::isNA(priors(7))) {
      // std::cout << "Using spatiotemporal PC prior" << "\n";
      jnll -= sdmTMB::pc_prior_matern(
          ln_tau_E(m), ln_kappa(1,m),
          priors(4), priors(5), priors(6), priors(7),
          true, /* log */
          share_range(m), stan_flag);
    }
    if (!sdmTMB::isNA(priors(8))) { // phi
      jnll -= dnorm(phi(m), priors(8), priors(9), true);
      if (stan_flag) jnll -= ln_phi(m); // Jacobian adjustment
    }
    if (!sdmTMB::isNA(priors(10))) { // AR1 random field rho
      jnll -= dnorm(rho(m), priors(10), priors(11), true);
      // Jacobian adjustment:
      // transform = 2 * (e^x/(1+e^x)) - 1
      // https://www.wolframalpha.com/input?i=2+*+%28e%5Ex%2F%281%2Be%5Ex%29%29+-+1
      // log abs derivative = log((2 * exp(x)) / (1 + exp(x))^2)
      if (stan_flag) jnll -= log(2.) + ar1_phi(m) - 2. * log(1. + exp(ar1_phi(m)));
    }
  }

  // ------------------ Predictions on new data --------------------------------

  if (do_predict) {
    int n_p = proj_X_ij(0).rows(); // n 'p'redicted newdata
    int n_p_mesh = proj_mesh.rows(); // n 'p'redicted mesh (less than n_p if duplicate locations)
    array<Type> proj_fe(n_p, n_m);

    for (int m = 0; m < n_m; m++) {
      if (m == 0) proj_fe.col(m) = proj_X_ij(m) * b_j;
      if (n_m == 1) proj_fe.col(m) += proj_offset_i;
      if (m == 1) proj_fe.col(m) = proj_X_ij(m) * b_j2 + proj_offset_i;
    }

    // add threshold effect if specified
    if (threshold_func > 0) {
      if (threshold_func == 1) {
        // linear
        for (int m = 0; m < n_m; m++) {
          for (int i = 0; i < n_p; i++) {
            // TODO: does proj_X_threshold(i) need to be dimensioned by model?
            proj_fe(i,m) += sdmTMB::linear_threshold(proj_X_threshold(i), s_slope(m), s_cut(m));
          }
        }
      } else {
        // logistic
        for (int m = 0; m < n_m; m++) {
          for (int i = 0; i < n_p; i++) {
            // TODO: does proj_X_threshold(i) need to be dimensioned by model?
            proj_fe(i,m) += sdmTMB::logistic_threshold(proj_X_threshold(i), s50(m), s95(m), s_max(m));
          }
        }
      }
    }

    // Smoothers:
    array<Type> proj_smooth_i(n_p, n_m);
    proj_smooth_i.setZero();
    if (has_smooths) {
      for (int m = 0; m < n_m; m++) {
        for (int s = 0; s < b_smooth_start.size(); s++) { // iterate over # of smooth elements
          array<Type> beta_s(proj_Zs(s).cols(),n_m);
          beta_s.setZero();
          for (int j = 0; j < beta_s.rows(); j++) {
            beta_s(j,m) = b_smooth(b_smooth_start(s) + j,m);
          }
          proj_smooth_i.col(m) += proj_Zs(s) * vector<Type>(beta_s.col(m));
        }
        proj_smooth_i.col(m) += proj_Xs * vector<Type>(bs.col(m));
        proj_fe.col(m) += proj_smooth_i.col(m);
      }
    }

    // IID random intercepts:
    array<Type> proj_iid_re_i(n_p,n_m);
    proj_iid_re_i.setZero();
    for (int m = 0; m < n_m; m++) {
      for (int i = 0; i < n_p; i++) {
        int temp = 0;
        for (int k = 0; k < n_RE; k++) {
          if (k == 0 && !exclude_RE(0)) proj_iid_re_i(i,m) += RE(proj_RE_indexes(i, k),m);
          if (k > 0) {
            temp += nobs_RE(k - 1);
            if (!exclude_RE(k)) proj_iid_re_i(i,m) += RE(proj_RE_indexes(i, k) + temp,m);
          }
        }
        proj_fe(i) += proj_iid_re_i(i,m);
      }
    }

    // Random walk covariates:
    array<Type> proj_rw_i(n_p,n_m);
    proj_rw_i.setZero();
    if (random_walk == 1 || ar1_time || random_walk == 2) {
      for (int m = 0; m < n_m; m++) {
        for (int i = 0; i < proj_X_rw_ik.rows(); i++) {
          for (int k = 0; k < proj_X_rw_ik.cols(); k++) {
            proj_rw_i(i,m) += proj_X_rw_ik(i, k) * b_rw_t(proj_year(i), k, m);
            proj_fe(i) += proj_rw_i(i,m);
          }
        }
      }
    }

    // Spatial and spatiotemporal random fields (by unique location):
    array<Type> proj_omega_s_A_unique(n_p_mesh, n_m);
    array<Type> proj_zeta_s_A_unique(n_p_mesh, n_z, n_m);
    array<Type> proj_epsilon_st_A_unique(n_p_mesh, n_t, n_m);
    proj_epsilon_st_A_unique.setZero();

    // Expanded to full length:
    array<Type> proj_omega_s_A(n_p, n_m);
    array<Type> proj_zeta_s_A(n_p, n_z, n_m);
    array<Type> proj_epsilon_st_A_vec(n_p, n_m);
    proj_omega_s_A.setZero(); // may not get filled
    proj_zeta_s_A.setZero(); // may not get filled
    proj_epsilon_st_A_vec.setZero(); // may not get filled

    array<Type> proj_zeta_s_A_cov(n_p, n_z, n_m);
    proj_zeta_s_A_cov.setZero();

    if (!no_spatial) {
      for (int m = 0; m < n_m; m++) {
        for (int t = 0; t < n_t; t++) {
          proj_epsilon_st_A_unique.col(m).col(t) = proj_mesh * vector<Type>(epsilon_st.col(m).col(t));
        }
        if (!omit_spatial_intercept) proj_omega_s_A_unique.col(m) = proj_mesh * vector<Type>(omega_s.col(m));
      }

      // Spatially varying coefficients:
      if (spatial_covariate) {
        for (int m = 0; m < n_m; m++) {
          for (int z = 0; z < n_z; z++) {
            proj_zeta_s_A_unique.col(m).col(z) = proj_mesh * vector<Type>(zeta_s.col(m).col(z));
          }
        }
      }

      // Pick out the appropriate spatial and/or or spatiotemporal values:
      for (int m = 0; m < n_m; m++) {
        for (int i = 0; i < n_p; i++) {
          proj_omega_s_A(i,m) = proj_omega_s_A_unique(proj_spatial_index(i),m);
          proj_epsilon_st_A_vec(i,m) = proj_epsilon_st_A_unique(proj_spatial_index(i), proj_year(i),m);
          for (int z = 0; z < n_z; z++) {
            proj_zeta_s_A(i,z,m) = proj_zeta_s_A_unique(proj_spatial_index(i),z,m);
          }
        }
      }

      if (spatial_covariate) {
        for (int m = 0; m < n_m; m++) {
          for (int z = 0; z < n_z; z++) {
            for (int i = 0; i < n_p; i++) {
              proj_zeta_s_A_cov(i,z,m) = proj_zeta_s_A(i,z,m) * proj_z_i(i,z);
            }
          }
        }
      }
    }

    // for (int m = 0; m < n_m; m++) {
    //   for (int i = 0; i < n_i; i++) {
    //     if ((n_m == 2 && m == 2) || n_m == 1) proj_fe(i,m) += proj_offset_i(i);
    //   }
    // }

    array<Type> proj_rf(n_p, n_m);
    array<Type> proj_eta(n_p, n_m);
    for (int m = 0; m < n_m; m++)
      proj_rf.col(m) = proj_omega_s_A.col(m) + proj_epsilon_st_A_vec.col(m);

    for (int m = 0; m < n_m; m++)
      for (int z = 0; z < n_z; z++)
        proj_rf.col(m) += proj_zeta_s_A_cov.col(m).col(z); // SVC effects

      // proj_fe includes s(), RW, and IID random effects:
    for (int m = 0; m < n_m; m++)
      proj_eta.col(m) = proj_fe.col(m) + proj_rf.col(m);

    if (n_m > 1 && pop_pred) { // grab SE on fixed effects combined if delta model:
      Type t1, t2;
      vector<Type> proj_rf_delta(n_p);
      for (int i = 0; i < n_p; i++) {
        if (poisson_link_delta) {
          proj_rf_delta(i) = proj_fe(i,0) + proj_fe(i,1); // check
        } else {
          t1 = InverseLink(proj_fe(i,0), link(0));
          t2 = InverseLink(proj_fe(i,1), link(1));
          proj_rf_delta(i) = Link(t1 * t2, link(1));
        }
      }
      if (calc_se) ADREPORT(proj_rf_delta);
    }

    // FIXME save memory by not reporting all these or optionally so for MVN/Bayes?
    REPORT(proj_fe);            // fixed effect projections
    REPORT(proj_omega_s_A);     // spatial random effect projections
    REPORT(proj_epsilon_st_A_vec);  // spatiotemporal random effect projections
    REPORT(proj_zeta_s_A);      // spatial slope projections
    // REPORT(proj_zeta_s_A_cov);  // spatial slope * covariate projections
    REPORT(proj_eta);           // combined projections (in link space)
    REPORT(proj_rf);            // combined random field projections
    REPORT(proj_rw_i);          // random walk projections
    REPORT(proj_iid_re_i);      // IID random intercept projections

    if (calc_se) {
      if (pop_pred) {
        ADREPORT(proj_fe);
      } else {
        ADREPORT(proj_eta);
      }
    }

    // Total biomass etc.:
    vector<Type> total(n_t);
    total.setZero();
    vector<Type> mu_combined(n_p);
    mu_combined.setZero();

    if (calc_index_totals || calc_cog) {
      // ------------------ Derived quantities ---------------------------------
      Type t1;
      Type t2;
      int link_tmp;

      if (n_m > 1) { // delta model
        for (int i = 0; i < n_p; i++) {
          if (poisson_link_delta) {
            // Type R1 = Type(1.) - exp(-exp(proj_eta(i,0)));
            // Type R2 = exp(proj_eta(i,0)) / R1 * exp(proj_eta(i,1))
            mu_combined(i) = exp(proj_eta(i,0) + proj_eta(i,1)); // prevent numerical issues
          } else {
            t1 = InverseLink(proj_eta(i,0), link(0));
            t2 = InverseLink(proj_eta(i,1), link(1));
            mu_combined(i) = t1 * t2;
          }
          total(proj_year(i)) += mu_combined(i) * area_i(i);
        }
      } else { // non-delta model
        for (int i = 0; i < n_p; i++) {
          mu_combined(i) = InverseLink(proj_eta(i,0), link(0));
          total(proj_year(i)) += mu_combined(i) * area_i(i);
        }
      }
      vector<Type> link_total(n_t);
      if (n_m > 1) {
        link_tmp = link(1); // 2nd link should always be log/exp in this case
      } else {
        link_tmp = link(0);
      }
      for (int i = 0; i < n_t; i++) {
        link_total(i) = Link(total(i), link_tmp);
      }
      if (calc_index_totals) {
        REPORT(link_total);
        ADREPORT(link_total);
        ADREPORT(total);
      }

      // Low-rank sparse hessian bias-correction
      PARAMETER_VECTOR(eps_index);
      if (eps_index.size() > 0) {
        Type S;
        for (int t=0; t < n_t; t++) {
          S = total(t); // Set lowrank tag on S = sum(exp(x))
          jnll += eps_index(t) * S;
        }
      }
      if (calc_cog) {
        // Centre of gravity:
        vector<Type> cog_x(n_t);
        vector<Type> cog_y(n_t);
        cog_x.setZero();
        cog_y.setZero();
        for (int i = 0; i < n_p; i++) {
          cog_x(proj_year(i)) += proj_lon(i) * mu_combined(i) * area_i(i);
          cog_y(proj_year(i)) += proj_lat(i) * mu_combined(i) * area_i(i);
        }
        for (int t = 0; t < n_t; t++) {
          cog_x(t) /= total(t);
          cog_y(t) /= total(t);
        }
        REPORT(cog_x);
        ADREPORT(cog_x);
        REPORT(cog_y);
        ADREPORT(cog_y);
      }
    }
  }

   if (threshold_func == 1) { // linear breakpoint model
     REPORT(s_slope);
     ADREPORT(s_slope);
     REPORT(s_cut);
     ADREPORT(s_cut);
   }
   if (threshold_func == 2) { // logistic function model
     REPORT(s50);
     ADREPORT(s50);
     REPORT(s95);
     ADREPORT(s95);
     REPORT(s_max);
     ADREPORT(s_max);
   }
//    if (calc_quadratic_range && b_j(1) < Type(0)) {
//      vector<Type> quadratic_roots = sdmTMB::GetQuadraticRoots(b_j(1), b_j(0), Type(0.05));
//      Type quadratic_low = quadratic_roots(0);
//      Type quadratic_hi = quadratic_roots(1);
//      Type quadratic_range = quadratic_roots(1) - quadratic_roots(0);
//      if (quadratic_range < 0) quadratic_range = quadratic_range * -1.;
//      Type quadratic_peak = quadratic_roots(2);
//      Type quadratic_reduction = quadratic_roots(3);
//
//      REPORT(quadratic_low);
//      REPORT(quadratic_hi);
//      REPORT(quadratic_range);
//      REPORT(quadratic_peak);
//      REPORT(quadratic_reduction);
//
//      ADREPORT(quadratic_low);
//      ADREPORT(quadratic_hi);
//      ADREPORT(quadratic_range);
//      ADREPORT(quadratic_peak);
//      ADREPORT(quadratic_reduction);
//    }
   if (est_epsilon_slope) {
     REPORT(b_epsilon);
     ADREPORT(b_epsilon);
   }
  if (est_epsilon_re) {
    REPORT(ln_epsilon_re_sigma);
    ADREPORT(ln_epsilon_re_sigma);
  }

  //  // ------------------ Reporting ----------------------------------------------
  // FIXME save memory by not reporting all these or optionally so for MVN/Bayes?

  array<Type> log_range(range.rows(),range.cols()); // for SE
  for (int i = 0; i < range.rows(); i++) {
    for (int m = 0; m < range.cols(); m++) {
      log_range(i,m) = log(range(i,m));
    }
  }
  array<Type> log_sigma_E(sigma_E.rows(),sigma_E.cols()); // for SE
  for (int i = 0; i < sigma_E.rows(); i++) {
    for (int m = 0; m < sigma_E.cols(); m++) {
      log_sigma_E(i,m) = log(sigma_E(i,m));
    }
  }

  ADREPORT(logit_p_mix);
  ADREPORT(log_ratio_mix);
  REPORT(p_mix);
  REPORT(mix_ratio);
  ADREPORT(phi);
  ADREPORT(log_sigma_E);      // log spatio-temporal SD
  REPORT(sigma_E);      // spatio-temporal SD
  ADREPORT(sigma_E);      // log spatio-temporal SD
  REPORT(epsilon_st_A_vec);   // spatio-temporal effects; vector
  REPORT(b_rw_t);   // time-varying effects
  REPORT(omega_s_A);      // spatial effects; n_s length vector
  REPORT(zeta_s_A);     // spatial covariate effects; n_s length vector
  REPORT(eta_fixed_i);  // fixed effect predictions in the link space
  REPORT(eta_smooth_i); // smooth effect predictions in the link space
  REPORT(eta_i);        // fixed and random effect predictions in link space
  REPORT(eta_rw_i);     // time-varying predictions in link space
  REPORT(eta_iid_re_i); // IID intercept random effect estimates
  REPORT(rho);          // AR1 correlation in -1 to 1 space
  REPORT(range);        // Matern approximate distance at 10% correlation
  ADREPORT(range);        // Matern approximate distance at 10% correlation
  REPORT(log_range);  // log Matern approximate distance at 10% correlation
  ADREPORT(log_range);  // log Matern approximate distance at 10% correlation
  REPORT(b_smooth);     // smooth coefficients for penalized splines
  REPORT(ln_smooth_sigma); // standard deviations of smooth random effects, in log-space
  SIMULATE {
    REPORT(y_i);
    REPORT(omega_s);
    REPORT(omega_s_A);
    REPORT(epsilon_st);
    REPORT(epsilon_st_A_vec);
    REPORT(zeta_s);
    REPORT(zeta_s_A);
  }
  return jnll;
}
