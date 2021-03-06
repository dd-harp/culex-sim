/*
 * rates.hpp
 * rates, photoperiod, and temperature
 * mathematical model by David Ewing
 * stochastic adaption by Sean L. Wu
 * April 2022
*/

#ifndef RATES_HPP
#define RATES_HPP

#include <functional>
#include <RcppArmadillo.h>

// --------------------------------------------------------------------------------
//  external forcing (temp, photoperiod)
// --------------------------------------------------------------------------------

// temperature as modified cosine function
inline double temperature_ewing(const double t, const Rcpp::List& pars) {
  
  double phi = pars["phi"]; // PHASE
  double lambda = pars["lambda"]; // A
  double mu = pars["mu"]; // M
  double gamma = pars["gamma"]; // POWER
  
  double temp = 0.0;
  
  if (t < 0.0) {
    temp = (mu - lambda) + lambda * 2.0 * pow(0.5 * (1.0 + cos(2.0 * M_PI * (0.0 - phi) / 365.0)), gamma);
  } else {
    temp = (mu - lambda) + lambda * 2.0 * pow(0.5 * (1.0 + cos(2.0 * M_PI * (t - phi) / 365.0)), gamma);
  }
  
  return temp;
}

// photoperiod
inline double photoperiod(const double t, const Rcpp::List& pars) {
  
  double L = pars["L"]; // latitude (51 in thesis)
  
  // define photoperiod values
  double EPS = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (t - 3.5)))));
  double NUM = sin(0.8333 * M_PI/ 180.0) + (sin(L * M_PI / 180.0) * sin(EPS));
  double DEN = cos(L * M_PI / 180.0) * cos(EPS);
  double DAYLIGHT = 24.0 - (24.0 / M_PI) * acos(NUM / DEN);
  
  return DAYLIGHT;
}


// --------------------------------------------------------------------------------
//   diapause
// --------------------------------------------------------------------------------

// pp: photoperiod
inline double diapause_spring_ewing(const double pp){
  return 1.0 / (1.0 + exp(5.0 * (14.0 - pp)));
}

inline double diapause_autumn_ewing(const double pp) {
  return 1.0 / (1.0 + exp(5.0 * (13.0 - pp)));
}


// --------------------------------------------------------------------------------
//   mortality
// --------------------------------------------------------------------------------

// for each life stage, mortality is computed; returns vector of rates, input is temperature
using mortality_rate_fn = std::function<arma::Col<double>(const arma::Col<double>&)>;

mortality_rate_fn make_egg_mortality_rate_ewing(const Rcpp::List& pars);

mortality_rate_fn make_larvae_mortality_rate_ewing(const Rcpp::List& pars);

mortality_rate_fn make_pupae_mortality_rate_ewing(const Rcpp::List& pars);

mortality_rate_fn make_adult_mortality_rate_ewing(const Rcpp::List& pars);

// // egg mortality
// inline double death_egg_rate(const double temp, const Rcpp::List& pars) {
//   
//   double nu_0E = pars["nu_0E"]; // U3
//   double nu_1E = pars["nu_1E"]; // U4
//   double nu_2E = pars["nu_2E"]; // U5
//   double death_max = pars["death_max"];
//   
//   // calculate egg death rate
//   double egg_d = nu_0E * exp(pow((temp - nu_1E) / nu_2E, 2.0));
//   
//   if (egg_d > death_max) {
//     egg_d = death_max;
//   }
//   
//   return egg_d;
// }
// 
// // larvae mortality
// inline double death_larvae_rate(const double temp, const Rcpp::List& pars) {
//   
//   double nu_0L = pars["nu_0L"]; // U3
//   double nu_1L = pars["nu_1L"]; // U4
//   double nu_2L = pars["nu_2L"]; // U5
//   double death_max = pars["death_max"];
//   
//   // calculate egg death rate
//   double larvae_d = nu_0L * exp(pow((temp - nu_1L) / nu_2L, 2.0));
//   
//   if (larvae_d > death_max) {
//     larvae_d = death_max;
//   }
//   
//   return larvae_d;
// }
// 
// // pupal mortality
// inline double death_pupae_rate(const double temp, const Rcpp::List& pars) {
//   
//   double nu_0P = pars["nu_0P"]; // U3
//   double nu_1P = pars["nu_1P"]; // U4
//   double nu_2P = pars["nu_2P"]; // U5
//   double death_max = pars["death_max"];
//   
//   // calculate egg death rate
//   double pupal_d = nu_0P * exp(pow((temp - nu_1P)/nu_2P, 2.0));
//   
//   if (pupal_d > death_max) {
//     pupal_d = death_max;
//   }
//   
//   return pupal_d;
// }
// 
// // adult mortality
// inline double death_adult_rate(const double temp, const Rcpp::List& pars) {
//   
//   double alpha_A = pars["alpha_A"]; // ALPHA
//   double beta_A = pars["beta_A"]; // # BETA
//   double death_min_a = pars["death_min_a"];
//   
//   // calculate adult death rate
//   double adult_d = alpha_A * pow(temp, beta_A);
//   
//   if(adult_d < death_min_a){
//     adult_d = death_min_a;
//   }
//   
//   return adult_d;
// }


// --------------------------------------------------------------------------------
//   gonotrophic rate
// --------------------------------------------------------------------------------

// gonotrophic cycle duration is directly proportional to feeding rate
// args: temperature, diapause, returns rates
using gonotrophic_fn = std::function<arma::Col<double>(const arma::Col<double>&, const arma::Col<double>&)>;

gonotrophic_fn make_gonotrophic_ewing(const Rcpp::List& pars);

// // G: duration of gonotrophic cycle
// inline double gonotrophic(const double temp, const Rcpp::List& pars) {
//   
//   double q1 = pars["q1"]; // KG
//   double q2 = pars["q2"]; // QG
//   double q3 = pars["q3"]; // BG
//   double gon_min = pars["gon_min"];
//   
//   // calculate gonotrophic cycle length
//   double grate;
//   if (temp < 0.0) {
//     grate = 0.0333;
//   } else {
//     grate = q1 / (1.0 + q2*exp(-q3*temp));
//   }
//   
//   if (grate < gon_min) {
//     grate = gon_min;
//   }
//   
//   return 1.0 / grate;
// }


// --------------------------------------------------------------------------------
//   oviposition
// --------------------------------------------------------------------------------

// per-capita oviposition rate (eggs/female/day)
// args: gonotrophic rate
using oviposition_fn = std::function<arma::Col<double>(const arma::Col<double>&)>;

oviposition_fn make_oviposition_ewing(const Rcpp::List& pars);

// // per-capita oviposition rate
// // d: diapause
// // G: duration of gonotrophic cycle
// // pars: parameters
// inline double oviposition(const double d, const double G, const Rcpp::List& pars) {
//   
//   double max_egg = pars["max_egg"];
//   
//   double egg_raft = d * max_egg * 0.5;
//   double ovi = egg_raft / G;
//   
//   return ovi;
// }


// --------------------------------------------------------------------------------
//   lifecycle stage progression rates
// --------------------------------------------------------------------------------

inline double egg_maturation_rate_ewing(const double temp, const Rcpp::List& pars) {
  
  double alpha_E = pars["alpha_E"]; // ALPHA
  double beta_E = pars["beta_E"]; // BETA
  double maturation_min = pars["maturation_min"];
  
  // calculate egg development rate
  double egg_maturation;
  if (temp < 0.0) {
    egg_maturation = 0.016667;
  } else {
    egg_maturation = alpha_E * pow(temp, beta_E);
  }
  
  if (egg_maturation < maturation_min) {
    egg_maturation = maturation_min;
  }
  
  return egg_maturation;
}

inline double larvae_maturation_rate_ewing(const double temp, const Rcpp::List& pars) {
  
  double alpha_L = pars["alpha_L"]; // ALPHA
  double beta_L = pars["beta_L"]; // BETA
  double maturation_min = pars["maturation_min"];
  
  // calculate larvae development rate
  double larvae_maturation;
  if (temp < 0.0) {
    larvae_maturation = 0.016667;
  } else {
    larvae_maturation = alpha_L * pow(temp, beta_L);
  }
  
  if (larvae_maturation < maturation_min) {
    larvae_maturation = maturation_min;
  }
  
  return larvae_maturation;
}

inline double pupae_maturation_rate_ewing(const double temp, const Rcpp::List& pars) {
  
  double alpha_P = pars["alpha_P"]; // ALPHA
  double beta_P = pars["beta_P"]; // BETA
  double maturation_min = pars["maturation_min"];
  
  // calculate larvae development rate
  double pupae_maturation;
  if (temp < 0.0) {
    pupae_maturation = 0.016667;
  } else {
    pupae_maturation = alpha_P * pow(temp, beta_P);
  }
  
  if (pupae_maturation < maturation_min) {
    pupae_maturation = maturation_min;
  }
  
  return pupae_maturation;
}


// --------------------------------------------------------------------------------
//   extrinsic incubation period
// --------------------------------------------------------------------------------

inline double eip_rate_ewing(const double temp, const Rcpp::List& pars) {
  
  double q = pars["eip_q"];
  double Tmax = pars["eip_tmax"];
  double Tmin = pars["eip_tmin"];
  double eip_min = pars["eip_min"];
  
  if (temp > Tmax) {
    return eip_min;
  } else if (temp < Tmin) {
    return eip_min;
  } else {
    double eip = q*temp*(temp - Tmin) * sqrt(Tmax - temp);
    return std::max(eip, eip_min);
  }
  
}

#endif