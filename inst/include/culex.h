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
  
  culex(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_, const arma::Mat<double>& psi_);
  ~culex() = default;
  
  void update(const Rcpp::List& parameters);
  
};

// constructor
template <typename T>
inline culex<T>::culex(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const double dt_, const arma::Mat<double>& psi_) :
  psi(psi_), step(0), p(p_), dt(dt_), tau_E(tau_E_), tau_L(tau_L_), tau_P(tau_P_)
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
  double temp = temperature(tnow, parameters);
  
  // photoperiod
  double pp = photoperiod(tnow, parameters);
  double pp_1 = photoperiod(tnow - 1.0, parameters);
  
  // gonotrophic cycle
  double gon = gonotrophic(temp, parameters);
  
  // mortality
  double death_egg = death_egg_rate(temp, parameters);
  double death_larvae = death_larvae_rate(temp, parameters);
  double death_pupae = death_pupae_rate(temp, parameters);
  double death_adult = death_adult_rate(temp, parameters);
  
  arma::Row<int> larvae_tot = arma::sum(this->L, 0);
  std::vector<double> death_larvae_tot(this->p, death_larvae);
  for (auto i = 0u; i < this->p; ++i) {
    death_larvae_tot[i] += p0 * larvae_tot(i) / (p1 + larvae_tot(i));
  }
  
  // diapause and egg laying
  double dia;
  if (pp > pp_1) {
    dia = diapause_spring(pp);
  } else {
    dia = diapause_autumn(pp);
  }
  
  // egg laying
  arma::Row<int> lambda(this->p, arma::fill::zeros);
  int i{0};
  lambda.for_each([&i, this, dia, gon, &parameters](arma::Row<int>::elem_type& val) {
    double lambda_mean = oviposition(dia, gon, parameters) * this->A(i) * this->dt;
    val = R::rpois(lambda_mean);
    i++;
  });
  
  // survival 
  this->E.for_each([death_egg = death_egg, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_egg * dt, 1.0, 0, 0));  
    }
  });
  
  i = 0;
  this->L.each_col([&death_larvae_tot, &i, dt = this->dt](arma::Col<int>& val) {
    double surv = R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0);
    val.for_each([surv](arma::Col<int>::elem_type& val) {
      if (val > 0) {
        val = R::rbinom(val, surv);
      }
    });
    i++;
  });
  
  this->P.for_each([death_pupae = death_pupae, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_pupae * dt, 1.0, 0, 0));  
    }
  });
  
  this->A.for_each([death_adult = death_adult, dt = this->dt](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, R::pexp(death_adult * dt, 1.0, 0, 0));  
    }
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
  double temp = temperature(tnow, parameters);
  
  // photoperiod
  double pp = photoperiod(tnow, parameters);
  double pp_1 = photoperiod(tnow - 1.0, parameters);
  
  // gonotrophic cycle
  double gon = gonotrophic(temp, parameters);
  
  // mortality
  double death_egg = death_egg_rate(temp, parameters);
  double death_larvae = death_larvae_rate(temp, parameters);
  double death_pupae = death_pupae_rate(temp, parameters);
  double death_adult = death_adult_rate(temp, parameters);
  
  arma::Row<double> larvae_tot = arma::sum(this->L, 0);
  std::vector<double> death_larvae_tot(this->p, death_larvae);
  for (auto i = 0u; i < this->p; ++i) {
    death_larvae_tot[i] += p0 * larvae_tot(i) / (p1 + larvae_tot(i));
  }
  
  // diapause and egg laying
  double dia;
  if (pp > pp_1) {
    dia = diapause_spring(pp);
  } else {
    dia = diapause_autumn(pp);
  }
  
  // egg laying
  arma::Row<double> lambda = oviposition(dia, gon, parameters) * this->A * this->dt;
  
  // survival
  double surv = R::pexp(death_egg * dt, 1.0, 0, 0);
  this->E *= surv;
  
  int i{0};
  this->L.each_col([&death_larvae_tot, &i, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });
  
  surv = R::pexp(death_pupae * dt, 1.0, 0, 0);
  this->P *= surv;
  
  surv = R::pexp(death_adult * dt, 1.0, 0, 0);
  this->A *= surv;
  
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
