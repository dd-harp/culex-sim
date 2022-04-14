#include "culex_types.h"

// ---------- deterministic model interface ----------

//' @title Create deterministic culex model object
//' @param p number of patches
//' @param tau_E vector of egg delays
//' @param tau_L vector of larvae delays
//' @param tau_P vector of pupae delays
//' @param dt size of time step
//' @param psi adult movement matrix
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<culex_deterministic> create_culex_deterministic(const int p, const std::vector<int>& tau_E, const std::vector<int>& tau_L, const std::vector<int>& tau_P, const double dt, const arma::Mat<double>& psi) {
  return Rcpp::XPtr<culex_deterministic>(
    new culex_deterministic(p, tau_E, tau_L, tau_P, dt, psi),
    true
  );
};

//' @title Step deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @param parameters a [list] of parameters
//' @export
// [[Rcpp::export]]
void step_culex_deterministic(Rcpp::XPtr<culex_deterministic> mod, const Rcpp::List& parameters) {
  mod->update(parameters);
}

//' @title Set adults for deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @param A a vector
//' @export
// [[Rcpp::export]]
void set_A_deterministic(Rcpp::XPtr<culex_deterministic> mod, arma::Row<double> A) {
  mod->A = A;
};

//' @title Return adults for deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_A_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return mod->A;
};

//' @title Return eggs for deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_E_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->E, 0);
};

//' @title Return larvae for deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_L_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->L, 0);
};

//' @title Return pupae for deterministic culex model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_P_deterministic(Rcpp::XPtr<culex_deterministic> mod) {
  return arma::sum(mod->P, 0);
};