#include <RcppArmadillo.h>

#include "rates.h"

// --------------------------------------------------------------------------------
//   ODEs describing change in delay duration
// --------------------------------------------------------------------------------

//' @title Differential equations describing change in delays
//' @description This function is meant to be passed to [deSolve::ode].
//' @param t time
//' @param y state
//' @param params a list of parameters
//' @export
// [[Rcpp::export]]
Rcpp::List tau_diffeqn(const double t, const Rcpp::NumericVector& y, const Rcpp::List& params) {
  
  // state variables
  double DE = y[0]; // tau_E(t)
  double DL = y[1]; // tau_L(t)
  double DP = y[2]; // tau_P(t)
  double DEIP = y[3]; // tau_EIP(t)
  
  // temperature
  double temp = temperature_ewing(t, params);
  double temp_E = temperature_ewing(t - DE, params);
  double temp_L = temperature_ewing(t - DL, params);
  double temp_P = temperature_ewing(t - DP, params);
  double temp_EIP = temperature_ewing(t - DEIP, params);

  // development
  double larvae_maturation = larvae_maturation_rate_ewing(temp, params);
  double larvae_maturation_L = larvae_maturation_rate_ewing(temp_L, params);
        
  double egg_maturation = egg_maturation_rate_ewing(temp, params);
  double egg_maturation_E = egg_maturation_rate_ewing(temp_E, params);

  double pupae_maturation = pupae_maturation_rate_ewing(temp, params);
  double pupae_maturation_P = pupae_maturation_rate_ewing(temp_P, params);
  
  double incubation = eip_rate_ewing(temp, params);
  double incubation_EIP = eip_rate_ewing(temp_EIP, params);

  // ODEs describing change in state duration
  double dDEdt = 1.0 - egg_maturation/egg_maturation_E;
  double dDLdt = 1.0 - larvae_maturation/larvae_maturation_L;
  double dDPdt = 1.0 - pupae_maturation/pupae_maturation_P;
  double dDEIPdt = 1.0 - incubation/incubation_EIP;

  // ODE system
  Rcpp::NumericVector du(4);
  du[0] = dDEdt;
  du[1] = dDLdt;
  du[2] = dDPdt;
  du[3] = dDEIPdt;
  
  return Rcpp::List::create(du);
}