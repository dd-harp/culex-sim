#include "rates.h"


// --------------------------------------------------------------------------------
//   mortality
// --------------------------------------------------------------------------------

// for each life stage, mortality is computed; returns vector of rates, input is temperature
mortality_rate_fn make_egg_mortality_rate_ewing(const Rcpp::List& pars) {
  double nu_0E = Rcpp::as<double>(pars["nu_0E"]); // U3
  double nu_1E = Rcpp::as<double>(pars["nu_1E"]); // U4
  double nu_2E = Rcpp::as<double>(pars["nu_2E"]); // U5
  double death_max = Rcpp::as<double>(pars["death_max"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [nu_0E, nu_1E, nu_2E, death_max, n_patch] (const arma::Col<double>& temp) -> arma::Col<double> {
    arma::Col<double> mort(n_patch, arma::fill::zeros);
    for (auto i = 0u; i < n_patch; ++i) {
      mort(i) = std::min(nu_0E * exp(pow((temp(i) - nu_1E) / nu_2E, 2.0)), death_max);
    }
    return mort;
  };
};

mortality_rate_fn make_larvae_mortality_rate_ewing(const Rcpp::List& pars) {
  double nu_0L = Rcpp::as<double>(pars["nu_0L"]); // U3
  double nu_1L = Rcpp::as<double>(pars["nu_1L"]); // U4
  double nu_2L = Rcpp::as<double>(pars["nu_2L"]); // U5
  double death_max = Rcpp::as<double>(pars["death_max"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [nu_0L, nu_1L, nu_2L, death_max, n_patch] (const arma::Col<double>& temp) -> arma::Col<double> {
    arma::Col<double> mort(n_patch, arma::fill::zeros);
    for (auto i = 0u; i < n_patch; ++i) {
      mort(i) = std::min(nu_0L * exp(pow((temp(i) - nu_1L) / nu_2L, 2.0)), death_max);
    }
    return mort;
  };
};

mortality_rate_fn make_pupae_mortality_rate_ewing(const Rcpp::List& pars) {
  double nu_0P = Rcpp::as<double>(pars["nu_0P"]); // U3
  double nu_1P = Rcpp::as<double>(pars["nu_1P"]); // U4
  double nu_2P = Rcpp::as<double>(pars["nu_2P"]); // U5
  double death_max = Rcpp::as<double>(pars["death_max"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [nu_0P, nu_1P, nu_2P, death_max, n_patch] (const arma::Col<double>& temp) -> arma::Col<double> {
    arma::Col<double> mort(n_patch, arma::fill::zeros);
    for (auto i = 0u; i < n_patch; ++i) {
      mort(i) = std::min(nu_0P * exp(pow((temp(i) - nu_1P) / nu_2P, 2.0)), death_max);
    }
    return mort;
  };
};

mortality_rate_fn make_adult_mortality_rate_ewing(const Rcpp::List& pars) {
  double alpha_A = Rcpp::as<double>(pars["alpha_A"]); // ALPHA
  double beta_A = Rcpp::as<double>(pars["beta_A"]); // # BETA
  double death_min_a = Rcpp::as<double>(pars["death_min_a"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [alpha_A, beta_A, death_min_a, n_patch] (const arma::Col<double>& temp) -> arma::Col<double> {
    arma::Col<double> mort(n_patch, arma::fill::zeros);
    for (auto i = 0u; i < n_patch; ++i) {
      mort(i) = std::max(alpha_A * pow(temp(i), beta_A), death_min_a);
    }
    return mort;
  };
};


// --------------------------------------------------------------------------------
//   gonotrophic rate
// --------------------------------------------------------------------------------

gonotrophic_fn make_gonotrophic_ewing(const Rcpp::List& pars) {
  double q1 = Rcpp::as<double>(pars["q1"]); // KG
  double q2 = Rcpp::as<double>(pars["q2"]); // QG
  double q3 = Rcpp::as<double>(pars["q3"]); // BG
  double gon_zero_temp = Rcpp::as<double>(pars["gon_zero_temp"]);
  double gon_min = Rcpp::as<double>(pars["gon_min"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [q1, q2, q3, gon_min, gon_zero_temp, n_patch] (const arma::Col<double>& temp, const arma::Col<double>& diapause) -> arma::Col<double> {
    arma::Col<double> gon(n_patch, arma::fill::zeros);
    for (auto i = 0u; i < n_patch; ++i) {
      if (temp(i) < 0.0) {
        gon(i) = gon_zero_temp;
      } else {
        gon(i) = q1 / (1.0 + q2*exp(-q3*temp(i)));
      }
      
      if (gon(i) < gon_min) {
        gon(i) = gon_min;
      }
      
      gon(i) *= diapause(i);
    }
    return gon;
  };
};


// --------------------------------------------------------------------------------
//   oviposition
// --------------------------------------------------------------------------------

oviposition_fn make_oviposition_ewing(const Rcpp::List& pars) {
  double max_egg = Rcpp::as<double>(pars["max_egg"]);
  size_t n_patch = Rcpp::as<size_t>(pars["n_patch"]);
  
  return [max_egg, n_patch] (const arma::Col<double>& gonotrophic) -> arma::Col<double> {
    return gonotrophic * max_egg * 0.5;
  };
};