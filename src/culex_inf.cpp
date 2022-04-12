#include "culex_types.h"

// ---------- stochastic model interface ----------

//' @title Create stochastic culex model object
//' @param p number of patches
//' @param tau_E vector of egg delays
//' @param tau_L vector of larvae delays
//' @param tau_P vector of pupae delays
//' @param dt size of time step
//' @param psi adult movement matrix
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<culex_infection_stochastic> create_culex_infection_stochastic(
    const int p, 
    const std::vector<int>& tau_E, const std::vector<int>& tau_L, const std::vector<int>& tau_P, const std::vector<int>& tau_EIP,
    const double dt, const arma::Mat<double>& psi, const int n_species)
{
  return Rcpp::XPtr<culex_infection_stochastic>(
    new culex_infection_stochastic(p, tau_E, tau_L, tau_P, tau_EIP, dt, psi, n_species),
    true
  );
};


// void set_f(const arma::Row<double>& f_) {f = f_;};
// void set_q(const arma::Mat<double> q_) {q = q_;};
// void set_kappa(const arma::Mat<double> kappa_) {kappa = kappa_;};
// 
// void set_E(const arma::Mat<T>& E_) {E = E_;};
// void set_L(const arma::Mat<T>& L_) {L = L_;};
// void set_P(const arma::Mat<T>& P_) {P = P_;};
// void get_E() {return E;};
// void get_L() {return L;};
// void get_P() {return P;};
// 
// void set_E_I(const arma::Mat<T>& E_I_) {E_I = E_I_;};
// void set_L_I(const arma::Mat<T>& L_I_) {L_I = L_I_;};
// void set_P_I(const arma::Mat<T>& P_I_) {P_I = P_I_;};
// void get_E_I() {return E_I;};
// void get_L_I() {return L_I;};
// void get_P_I() {return P_I;};
// 
// void set_A_S(const arma::Row<T>& A_S_) {A_S = A_S_;};
// void set_A_E(const arma::Mat<T>& A_E_) {A_E = A_E_;};
// void set_A_I(const arma::Row<T>& A_I_) {A_I = A_I_;};
// void get_A_S() {return A_S;};
// void get_A_E() {return A_E;};
// void get_A_I() {return A_I;};