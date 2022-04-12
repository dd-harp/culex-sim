/*
 * culex.hpp
 * model with infection states
 * mathematical model by David Ewing
 * stochastic adaption by Sean L. Wu
 * March 2022
 */

#ifndef CULEX_INF_HPP
#define CULEX_INF_HPP

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
// E_I: infected eggs (steps X patches)
// L_I: infected larvae (steps X patches)
// P_I: infected pupae (steps X patches)
// A_S: adults (patches)
// A_E: adults (steps X patches)
// A_I: adults (patches)
// psi: a movement matrix; rows should sum to 1 (patches X patches)
// shiftE: sparse matrix when multiplies E/E_I on the left to advance by 1 step
// shiftL: sparse matrix when multiplies L/L_I on the left to advance by 1 step
// shiftP: sparse matrix when multiplies P/P_I on the left to advance by 1 step
// shiftEIP: sparse matrix when multiplies A_E on the left to advance by 1 step
// step: current step
// p: number of patches
// dt: size of time step
// f: feeding rate (patches)
// q: proportion of bites on each host species (species X patches)
// kappa: net infectiousness of each host species (species X patches)
// tau_E: tau_E[t] gives the number of steps required for a new E at step t to mature
// tau_L: tau_L[t] gives the number of steps required for a new L at step t to mature
// tau_P: tau_P[t] gives the number of steps required for a new P at step t to mature
// tau_EIP: tau_EIP[t] gives the number of steps required for a new exposed mosquito A_E at step t to become infectious
template <typename T>
struct culex_inf {
  
  arma::Mat<T> E;
  arma::Mat<T> L;
  arma::Mat<T> P;
  
  arma::Mat<T> E_I;
  arma::Mat<T> L_I;
  arma::Mat<T> P_I;
  
  arma::Row<T> A_S;
  arma::Mat<T> A_E;
  arma::Row<T> A_I;
  
  arma::Mat<double> psi;
  
  arma::SpMat<int> shiftE;
  arma::SpMat<int> shiftL;
  arma::SpMat<int> shiftP;
  arma::SpMat<int> shiftEIP;
  
  int step;
  int p;
  double dt;
  
  // feeding behavior and transmission
  arma::Row<double> f;
  arma::Mat<double> q;
  arma::Mat<double> kappa;
  
  // delay vectors
  std::vector<int> tau_E;
  std::vector<int> tau_L;
  std::vector<int> tau_P;
  std::vector<int> tau_EIP;
  
  // ctor/dtor
  culex_inf(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const std::vector<int>& tau_EIP_, const double dt_, const arma::Mat<double>& psi_, const int n_species);
  ~culex_inf() = default;
  
  // methods
  void update(const Rcpp::List& parameters);
  
};

// constructor
template <typename T>
inline culex_inf<T>::culex_inf(const int p_, const std::vector<int>& tau_E_, const std::vector<int>& tau_L_, const std::vector<int>& tau_P_, const std::vector<int>& tau_EIP_, const double dt_, const arma::Mat<double>& psi_, const int n_species) :
  psi(psi_), step(0), p(p_), dt(dt_), 
  f(n_species, arma::fill::zeros), q(n_species, p_, arma::fill::zeros), kappa(n_species, p_, arma::fill::zeros),
  tau_E(tau_E_), tau_L(tau_L_), tau_P(tau_P_), tau_EIP(tau_EIP_)
{
  int maxE = *max_element(tau_E.begin(), tau_E.end());
  int maxL = *max_element(tau_L.begin(), tau_L.end());
  int maxP = *max_element(tau_P.begin(), tau_P.end());
  int maxEIP = *max_element(tau_EIP.begin(), tau_EIP.end());
  
  // state
  E = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P = arma::Mat<T>(maxP, p, arma::fill::zeros);
  
  E_I = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L_I = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P_I = arma::Mat<T>(maxP, p, arma::fill::zeros);
  
  A_S = arma::Row<T>(p, arma::fill::zeros);
  A_E = arma::Mat<T>(maxEIP, p, arma::fill::zeros);
  A_I = arma::Row<T>(p, arma::fill::zeros);
  
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
  
  arma::umat fillEIP(2, maxEIP);
  for (auto i = 0u; i < (maxEIP - 1); ++i) {
    fillEIP(0, i) = i;
    fillEIP(1, i) = i+1;
  }
  arma::Col<int> fillEIPvals(maxEIP, arma::fill::ones);
  shiftEIP = arma::SpMat<int>(fillEIP, fillEIPvals, maxEIP, maxEIP);
  
};

// stochastic update
template <>
inline void culex_inf<int>::update(const Rcpp::List& parameters) {

  // for stochastic update, better to work with psi transposed
  // due to artifact of how to retrieve pointer to memory needed in multinomial sample
  if (this->step == 0) {
    this->psi = this->psi.t();
  }

  double tnow = this->step * this->dt;
  int tau_E = this->tau_E[this->step];
  int tau_L = this->tau_L[this->step];
  int tau_P = this->tau_P[this->step];
  int tau_EIP = this->tau_EIP[this->step];

  double p0 = parameters["p0"];
  double p1 = parameters["p1"];
  double pvt = parameters["pvt"];

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

  // larval mortality depends on total larvae (both S and I)
  arma::Row<int> larvae_tot = arma::sum(this->L, 0) + arma::sum(this->L_I, 0);
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
  lambda.for_each([&i, this, dia, gon, pvt, &parameters](arma::Row<int>::elem_type& val) {
    double mosy = this->A_S(i) + arma::accu(this->A_E.col(i)) + (this->A_I(i) * (1.0 - pvt));
    double lambda_mean = oviposition(dia, gon, parameters) * mosy * this->dt;
    val = R::rpois(lambda_mean);
    i++;
  });
  
  arma::Row<int> lambda_I(this->p, arma::fill::zeros);
  i = 0;
  lambda_I.for_each([&i, this, dia, gon, pvt, &parameters](arma::Row<int>::elem_type& val) {
    double lambda_mean = oviposition(dia, gon, parameters) * (this->A_I(i) * pvt) * this->dt;
    val = R::rpois(lambda_mean);
    i++;
  });
  
  // infections
  arma::Row<double> h = this->f % arma::sum(this->q % this->kappa, 0);
  arma::Row<int> S2E(this->p, arma::fill::zeros);
  for (auto i = 0u; i < this->p; ++i) {
    S2E(i) = R::rbinom(this->A_S(i), R::pexp(h(i) * this->dt, 1.0, 1, 0));
  }
  this->A_S -= S2E;
  
  // survival
  double surv = R::pexp(death_egg * dt, 1.0, 0, 0);
  this->E.for_each([surv](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  this->E_I.for_each([surv](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
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
  
  i = 0;
  this->L_I.each_col([&death_larvae_tot, &i, dt = this->dt](arma::Col<int>& val) {
    double surv = R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0);
    val.for_each([surv](arma::Col<int>::elem_type& val) {
      if (val > 0) {
        val = R::rbinom(val, surv);
      }
    });
    i++;
  });
  
  surv = R::pexp(death_pupae * dt, 1.0, 0, 0);
  this->P.for_each([surv](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  this->P_I.for_each([surv](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  
  surv = R::pexp(death_adult * dt, 1.0, 0, 0);
  this->A_S.for_each([surv](arma::Row<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  this->A_E.for_each([surv](arma::Mat<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  this->A_I.for_each([surv](arma::Row<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  S2E.for_each([surv](arma::Row<int>::elem_type& val) {
    if (val > 0) {
      val = R::rbinom(val, surv);
    }
  });
  
  // dispersal
  arma::Row<int> A_move(this->p, arma::fill::zeros);
  arma::Row<int> tmp(this->p, arma::fill::zeros);
  for (auto i = 0u; i < this->p; ++i) {
    R::rmultinom(this->A_S(i), this->psi.colptr(i), this->p, tmp.memptr());
    A_move += tmp;
  }
  this->A_S = A_move;
  
  for (auto j = 0u; j < this->A_E.n_rows; ++j) {
    A_move.zeros();
    tmp.zeros();
    for (auto i = 0u; i < this->p; ++i) {
      R::rmultinom(this->A_E(j, i), this->psi.colptr(i), this->p, tmp.memptr());
      A_move += tmp;
    }
    this->A_E.row(j) = A_move;
  }

  A_move.zeros();
  tmp.zeros();
  for (auto i = 0u; i < this->p; ++i) {
    R::rmultinom(this->A_I(i), this->psi.colptr(i), this->p, tmp.memptr());
    A_move += tmp;
  }
  this->A_I = A_move;
  
  A_move.zeros();
  tmp.zeros();
  for (auto i = 0u; i < this->p; ++i) {
    R::rmultinom(S2E(i), this->psi.colptr(i), this->p, tmp.memptr());
    A_move += tmp;
  }
  S2E = A_move;

  // advancement
  
  // uninfected immatures
  arma::Row<int> E2L = this->E.row(0);
  this->E.row(0).zeros();
  this->E = this->shiftE * this->E;
  this->E.row(tau_E-1) = lambda;
  
  arma::Row<int> L2P = this->L.row(0);
  this->L.row(0).zeros();
  this->L = this->shiftL * this->L;
  this->L.row(tau_L-1) = E2L;
  
  arma::Row<int> P2A = this->P.row(0);
  this->P.row(0).zeros();
  this->P = this->shiftP * this->P;
  this->P.row(tau_P-1) = L2P;;
  
  this->A_S += P2A;
  
  // infected immatures
  E2L = this->E_I.row(0);
  this->E_I.row(0).zeros();
  this->E_I = this->shiftE * this->E_I;
  this->E_I.row(tau_E-1) = lambda_I;
  
  L2P = this->L_I.row(0);
  this->L_I.row(0).zeros();
  this->L_I = this->shiftL * this->L_I;
  this->L_I.row(tau_L-1) = E2L;
  
  P2A = this->P_I.row(0);
  this->P_I.row(0).zeros();
  this->P_I = this->shiftP * this->P_I;
  this->P_I.row(tau_P-1) = L2P;
  
  this->A_I += P2A;
  
  // extrinsic incubation period
  arma::Row<int> E2I = this->A_E.row(0);
  this->A_E.row(0).zeros();
  this->A_E = this->shiftEIP * this->A_E;
  this->A_E.row(tau_EIP-1) = S2E;
  
  this->A_I += E2I;
  
  // advance time step
  this->step += 1;

};

// deterministic update
template <>
inline void culex_inf<double>::update(const Rcpp::List& parameters) {
  
  double tnow = this->step * this->dt;
  int tau_E = this->tau_E[this->step];
  int tau_L = this->tau_L[this->step];
  int tau_P = this->tau_P[this->step];
  int tau_EIP = this->tau_EIP[this->step];
  
  double p0 = parameters["p0"];
  double p1 = parameters["p1"];
  double pvt = parameters["pvt"];
  
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
  
  // larval mortality depends on total larvae (both S and I)
  arma::Row<double> larvae_tot = arma::sum(this->L, 0) + arma::sum(this->L_I, 0);
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
  arma::Row<double> lambda = this->A_S + arma::sum(this->A_E, 0) + (this->A_I * (1.0 - pvt));
  lambda *= oviposition(dia, gon, parameters) * this->dt;
  
  arma::Row<double> lambda_I = (this->A_I * pvt) * oviposition(dia, gon, parameters) * this->dt;
  
  // infections
  arma::Row<double> h = this->f % arma::sum(this->q % this->kappa, 0);
  arma::Row<double> S2E(this->p, arma::fill::zeros);
  for (auto i = 0u; i < this->p; ++i) {
    S2E(i) = this->A_S(i) * R::pexp(h(i) * this->dt, 1.0, 1, 0);
  }
  this->A_S -= S2E;
  
  // survival
  double surv = R::pexp(death_egg * dt, 1.0, 0, 0);
  this->E *= surv;
  this->E_I *= surv;
  
  int i{0};
  this->L.each_col([&death_larvae_tot, &i, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });
  
  i = 0;
  this->L_I.each_col([&death_larvae_tot, &i, dt = this->dt](arma::Col<double>& val) {
    double surv = R::pexp(death_larvae_tot[i] * dt, 1.0, 0, 0);
    val *= surv;
    i++;
  });
  
  surv = R::pexp(death_pupae * dt, 1.0, 0, 0);
  this->P *= surv;
  this->P_I *= surv;
  
  surv = R::pexp(death_adult * dt, 1.0, 0, 0);
  this->A_S *= surv;
  this->A_E *= surv;
  this->A_I *= surv;
  S2E *= surv;
  
  // dispersal
  this->A_S = this->A_S * this->psi;
  this->A_E = this->A_E * this->psi;
  this->A_I = this->A_I * this->psi;
  S2E *= this->psi;
  
  // advancement
  
  // uninfected immatures
  arma::Row<double> E2L = this->E.row(0);
  this->E.row(0).zeros();
  this->E = this->shiftE * this->E;
  this->E.row(tau_E-1) = lambda;
  
  arma::Row<double> L2P = this->L.row(0);
  this->L.row(0).zeros();
  this->L = this->shiftL * this->L;
  this->L.row(tau_L-1) = E2L;
  
  arma::Row<double> P2A = this->P.row(0);
  this->P.row(0).zeros();
  this->P = this->shiftP * this->P;
  this->P.row(tau_P-1) = L2P;;
  
  this->A_S += P2A;
  
  // infected immatures
  E2L = this->E_I.row(0);
  this->E_I.row(0).zeros();
  this->E_I = this->shiftE * this->E_I;
  this->E_I.row(tau_E-1) = lambda_I;
  
  L2P = this->L_I.row(0);
  this->L_I.row(0).zeros();
  this->L_I = this->shiftL * this->L_I;
  this->L_I.row(tau_L-1) = E2L;
  
  P2A = this->P_I.row(0);
  this->P_I.row(0).zeros();
  this->P_I = this->shiftP * this->P_I;
  this->P_I.row(tau_P-1) = L2P;
  
  this->A_I += P2A;
  
  // extrinsic incubation period
  arma::Row<double> E2I = this->A_E.row(0);
  this->A_E.row(0).zeros();
  this->A_E = this->shiftEIP * this->A_E;
  this->A_E.row(tau_EIP-1) = S2E;
  
  this->A_I += E2I;
  
  // advance time step
  this->step += 1;
  
};


#endif
