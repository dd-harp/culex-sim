/*
 * culex.hpp
 * struct to define the model and temperature/photoperiod functions (from https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90)
 * mathematical model by David Ewing
 * stochastic adaption by Sean L. Wu
 * March 2022
 */

#ifndef CULEX_HPP
#define CULEX_HPP

#include <Rmath.h>
#include <RcppArmadillo.h>

#include <algorithm> // for max_element
#include <vector>

#include "rates.h"


// --------------------------------------------------------------------------------
//   model struct
// --------------------------------------------------------------------------------

// fields:
// E: eggs (steps X patches)
// L: larvae (steps X patches)
// P: pupae (steps X patches)
// A: adults (patches)
// psi: a movement matrix; rows should sum to 1 (patches X patches)
// shiftE: sparse matrix when multiplies E on the left to advance by 1 step
// shiftL: sparse matrix when multiplies L on the left to advance by 1 step
// shiftP: sparse matrix when multiplies P on the left to advance by 1 step
// step: current step
// p: number of patches
// dt: size of time step
// tau_E: tau_E[t] gives the number of steps required for a new E at step t to mature
// tau_L: tau_L[t] gives the number of steps required for a new L at step t to mature
// tau_P: tau_P[t] gives the number of steps required for a new P at step t to mature
template <typename T>
struct culex {
  
  arma::Mat<T> E;
  arma::Mat<T> L;
  arma::Mat<T> P;
  arma::Row<T> A;
  
  arma::Mat<double> psi;
  
  arma::SpMat<int> shiftE;
  arma::SpMat<int> shiftL;
  arma::SpMat<int> shiftP;
  
  int step;
  int p;
  double dt;
  
  std::vector<int> tau_E;
  std::vector<int> tau_L;
  std::vector<int> tau_P;
  
  // rate functions
  mortality_rate_fn mort_egg;
  mortality_rate_fn mort_larvae;
  mortality_rate_fn mort_pupae;
  mortality_rate_fn mort_adult;
  
  gonotrophic_fn gonotrophic_cycle;
  oviposition_fn oviposition;
  
  // external forcing
  arma::Mat<double> temperature;
  arma::Mat<double> diapause;
  
  culex(const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_, const Rcpp::List& parameters);
  ~culex() = default;
  
  void update(const Rcpp::List& parameters);
  
};

// constructor
template <typename T>
inline culex<T>::culex(const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_, const Rcpp::List& parameters) :
  psi(Rcpp::as<arma::Mat<double>>(parameters["psi"])), step(0), p(Rcpp::as<int>(parameters["n_patch"])), 
  dt(dt_), tau_E(tau_E_), tau_L(tau_L_), tau_P(tau_P_),
  temperature(Rcpp::as<arma::Mat<double>>(parameters["temperature"])), diapause(Rcpp::as<arma::Mat<double>>(parameters["diapause"]))
{
  int maxE = *max_element(tau_E.begin(), tau_E.end());
  int maxL = *max_element(tau_L.begin(), tau_L.end());
  int maxP = *max_element(tau_P.begin(), tau_P.end());
  
  // state
  E = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P = arma::Mat<T>(maxP, p, arma::fill::zeros);
  A = arma::Row<T>(p, arma::fill::zeros);
  
  // shift matrices (multiply on left)
  arma::umat fillE(2, maxE);
  for (auto i = 0u; i < (maxE - 1); ++i) {
    fillE(0, i) = i;
    fillE(1, i) = i+1;
  }
  arma::Col<int> fillEvals(maxE, arma::fill::ones);
  shiftE = arma::SpMat<int>(fillE, fillEvals, maxE, maxE);
  
  arma::umat fillL(2, maxL);
  for (auto i = 0u; i < (maxL - 1); ++i) {
    fillL(0, i) = i;
    fillL(1, i) = i+1;
  }
  arma::Col<int> fillLvals(maxL, arma::fill::ones);
  shiftL = arma::SpMat<int>(fillL, fillLvals, maxL, maxL);
  
  arma::umat fillP(2, maxP);
  for (auto i = 0u; i < (maxP - 1); ++i) {
    fillP(0, i) = i;
    fillP(1, i) = i+1;
  }
  arma::Col<int> fillPvals(maxP, arma::fill::ones);
  shiftP = arma::SpMat<int>(fillP, fillPvals, maxP, maxP);
  
  // rate functions
  mort_egg = make_egg_mortality_rate_ewing(parameters);
  mort_larvae = make_larvae_mortality_rate_ewing(parameters);
  mort_pupae = make_pupae_mortality_rate_ewing(parameters);
  mort_adult = make_adult_mortality_rate_ewing(parameters);
  
  gonotrophic_cycle = make_gonotrophic_ewing(parameters);
  oviposition = make_oviposition_ewing(parameters);
  
};

// stochastic update
template <>
inline void culex<int>::update(const Rcpp::List& parameters) {
  
  // for stochastic update, better to work with psi transposed
  // due to artifact of how to retrieve pointer to memory needed in multinomial sample
  if (this->step == 0) {
    this->psi = this->psi.t();
  }
  
  double tnow = this->step * this->dt;
  int tau_E = this->tau_E[this->step];
  int tau_L = this->tau_L[this->step];
  int tau_P = this->tau_P[this->step];
  
  double p0 = parameters["p0"];
  double p1 = parameters["p1"];
  
  // temperature
  arma::Col<double> temp = this->temperature.col(step);
  
  // diapause
  arma::Col<double> dia = this->diapause.col(step);

  // gonotrophic cycle
  arma::Col<double> gon = this->gonotrophic_cycle(temp, dia);
  
  // oviposition
  arma::Col<double> ovi = this->oviposition(gon);
  
  // mortality
  arma::Col<double> death_egg = this->mort_egg(temp);
  arma::Col<double> death_larvae = this->mort_larvae(temp);
  arma::Col<double> death_pupae = this->mort_pupae(temp);
  arma::Col<double> death_adult = this->mort_adult(temp);
  
  // larvae mortality with density dependence
  arma::Row<int> larvae_tot = arma::sum(this->L, 0);
  arma::Col<double> death_larvae_dd = death_larvae;
  for (auto i = 0u; i < this->p; ++i) {
    death_larvae_dd[i] += p0 * larvae_tot(i) / (p1 + larvae_tot(i));
  }
  
  // egg laying
  arma::Row<int> lambda(this->p, arma::fill::zeros);
  int i{0};
  lambda.for_each([&i, this, &ovi](arma::Row<int>::elem_type& val) {
    double lambda_mean = ovi(i) * this->A(i) * this->dt;
    val = R::rpois(lambda_mean);
    i++;
  });
  
  // survival
  i = 0;
  this->E.each_col([&i, &death_egg, dt = this->dt](arma::Col<int>& val) {
    double surv = R::pexp(death_egg(i) * dt, 1.0, 0, 0);
    val.for_each([surv](arma::Col<int>::elem_type& val) {
      if (val > 0) {
        val = R::rbinom(val, surv);
      }
    });
    i++;
  });
  
  i = 0;
  this->L.each_col([&i, &death_larvae_dd, dt = this->dt](arma::Col<int>& val) {
    double surv = R::pexp(death_larvae_dd(i) * dt, 1.0, 0, 0);
    val.for_each([surv](arma::Col<int>::elem_type& val) {
      if (val > 0) {
        val = R::rbinom(val, surv);
      }
    });
    i++;
  });
  
  i = 0;
  this->P.each_col([&i, &death_pupae, dt = this->dt](arma::Col<int>& val) {
    double surv = R::pexp(death_pupae(i) * dt, 1.0, 0, 0);
    val.for_each([surv](arma::Col<int>::elem_type& val) {
      if (val > 0) {
        val = R::rbinom(val, surv);
      }
    });
    i++;
  });
  
  i = 0;
  this->A.for_each([&i, &death_adult, dt = this->dt](arma::Row<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_adult(i) * dt, 1.0, 0, 0));
    }
    i++;
  });
  
  // dispersal
  arma::Row<int> A_move(this->p, arma::fill::zeros);
  arma::Row<int> tmp(this->p, arma::fill::zeros);
  for (auto i = 0u; i < this->p; ++i) {
    R::rmultinom(this->A(i), this->psi.colptr(i), this->p, tmp.memptr());
    A_move += tmp;
  }
  this->A = A_move;
  
  // advancement
  arma::Row<int> E2L = this->E.row(0);
  this->E.row(0).zeros();
  this->E = this->shiftE * this->E;
  this->E.row(tau_E-1) += lambda;
  
  arma::Row<int> L2P = this->L.row(0);
  this->L.row(0).zeros();
  this->L = this->shiftL * this->L;
  this->L.row(tau_L-1) += E2L;
  
  arma::Row<int> P2A = this->P.row(0);
  this->P.row(0).zeros();
  this->P = this->shiftP * this->P;
  this->P.row(tau_P-1) += L2P;;
  
  this->A += P2A;
  
  this->step += 1;
  
};

// deterministic update
template <>
inline void culex<double>::update(const Rcpp::List& parameters) {
  
  double tnow = this->step * this->dt;
  int tau_E = this->tau_E[this->step];
  int tau_L = this->tau_L[this->step];
  int tau_P = this->tau_P[this->step];
  
  double p0 = parameters["p0"];
  double p1 = parameters["p1"];

  // temperature
  arma::Col<double> temp = this->temperature.col(step);
  
  // diapause
  arma::Col<double> dia = this->diapause.col(step);
  
  // gonotrophic cycle
  arma::Col<double> gon = this->gonotrophic_cycle(temp, dia);
  
  // oviposition
  arma::Col<double> ovi = this->oviposition(gon);
  
  // mortality
  arma::Col<double> death_egg = this->mort_egg(temp);
  arma::Col<double> death_larvae = this->mort_larvae(temp);
  arma::Col<double> death_pupae = this->mort_pupae(temp);
  arma::Col<double> death_adult = this->mort_adult(temp);
  
  // larvae mortality with density dependence
  arma::Row<double> larvae_tot = arma::sum(this->L, 0);
  arma::Col<double> death_larvae_dd = death_larvae;
  for (auto i = 0u; i < this->p; ++i) {
    death_larvae_dd(i) += p0 * larvae_tot(i) / (p1 + larvae_tot(i));
  }
  
  // egg laying
  arma::Row<double> lambda = (ovi.t() % this->A) * this->dt;
  
  // survival
  int i{0};
  this->E.each_col([&i, &death_egg, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_egg(i) * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });
  
  i = 0;
  this->L.each_col([&death_larvae_dd, &i, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_larvae_dd(i) * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });
  
  i = 0;
  this->P.each_col([&i, &death_pupae, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_pupae(i) * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });

  i = 0;
  this->A.for_each([&i, &death_adult, dt = this->dt](arma::Row<double>::elem_type& val) {
    double surv = R::pexp(death_adult(i) * dt, 1.0, 0 , 0);
    val *= surv;
    i++;
  });
  
  // dispersal
  this->A = this->A * this->psi;
  
  // advancement
  arma::Row<double> E2L = this->E.row(0);
  this->E.row(0).zeros();
  this->E = this->shiftE * this->E;
  this->E.row(tau_E-1) += lambda;
  
  arma::Row<double> L2P = this->L.row(0);
  this->L.row(0).zeros();
  this->L = this->shiftL * this->L;
  this->L.row(tau_L-1) += E2L;
  
  arma::Row<double> P2A = this->P.row(0);
  this->P.row(0).zeros();
  this->P = this->shiftP * this->P;
  this->P.row(tau_P-1) += L2P;;
  
  this->A += P2A;
  
  this->step += 1;
  
};


#endif
