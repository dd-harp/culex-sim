#include "culex_types.h"

// ---------- deterministic model interface ----------

//' @title Create deterministic culex infection model object
//' @param p number of patches
//' @param tau_E vector of egg delays
//' @param tau_L vector of larvae delays
//' @param tau_P vector of pupae delays
//' @param tau_EIP vector of extrinsic incubation period delays
//' @param dt size of time step
//' @param psi adult movement matrix
//' @param n_species number of host species
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<culex_infection_deterministic> create_culex_infection_deterministic(
    const int p, 
    const std::vector<int>& tau_E, const std::vector<int>& tau_L, const std::vector<int>& tau_P, const std::vector<int>& tau_EIP,
    const double dt, const arma::Mat<double>& psi, const int n_species)
{
  return Rcpp::XPtr<culex_infection_deterministic>(
    new culex_infection_deterministic(p, tau_E, tau_L, tau_P, tau_EIP, dt, psi, n_species),
    true
  );
};

//' @title Step deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param parameters a [list] of parameters
//' @export
// [[Rcpp::export]]
void step_culex_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, const Rcpp::List& parameters) {
  mod->update(parameters);
}

//' @title Set feeding rate for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param f a row vector (feeding rate by patch)
//' @export
// [[Rcpp::export]]
void set_f_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, const arma::Row<double>& f) {
  mod->f = f;
}

//' @title Get feeding rate for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_f_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->f;
} 

//' @title Set feeding habit for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param q a matrix (columns must sum to `1`; each column indicates proportion of bites allocated to each host species in that patch)
//' @export
// [[Rcpp::export]]
void set_q_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, const arma::Mat<double>& q) {
  mod->q = q;
}

//' @title Get feeding habit for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_q_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->q;
} 

//' @title Set kappa for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param kappa a matrix (each column indicates net infectiousness of each host species in that patch)
//' @export
// [[Rcpp::export]]
void set_kappa_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, const arma::Mat<double>& kappa) {
  mod->kappa = kappa;
}

//' @title Get kappa for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_kappa_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->kappa;
} 

//' @title Set susceptible adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a row vector
//' @export
// [[Rcpp::export]]
void set_AS_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Row<double>& A) {
  mod->A_S = A;
};

//' @title Get susceptible adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_AS_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->A_S;
};

//' @title Set incubating adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a matrix
//' @export
// [[Rcpp::export]]
void set_AE_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& A) {
  mod->A_E = A;
};

//' @title Get incubating adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_AE_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->A_E;
};

//' @title Set infectious adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a row vector
//' @export
// [[Rcpp::export]]
void set_AI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Row<double>& A) {
  mod->A_I = A;
};

//' @title Get infectious adults for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_AI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->A_I;
};

//' @title Set eggs for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param E a matrix
//' @export
// [[Rcpp::export]]
void set_E_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& E) {
  mod->E = E;
};

//' @title Get eggs for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_E_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->E;
};

//' @title Set infected eggs for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param E a matrix
//' @export
// [[Rcpp::export]]
void set_EI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& E) {
  mod->E_I = E;
};

//' @title Get infected eggs for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_EI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->E_I;
};

//' @title Set larvae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param L a matrix
//' @export
// [[Rcpp::export]]
void set_L_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& L) {
  mod->L = L;
};

//' @title Get larvae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_L_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->L;
};

//' @title Set infected larvae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param L a matrix
//' @export
// [[Rcpp::export]]
void set_LI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& L) {
  mod->L_I = L;
};

//' @title Get infected larvae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_LI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->L_I;
};

//' @title Set pupae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param P a matrix
//' @export
// [[Rcpp::export]]
void set_P_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& P) {
  mod->P = P;
};

//' @title Get pupae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_P_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->P;
};

//' @title Set infected pupae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param P a matrix
//' @export
// [[Rcpp::export]]
void set_PI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod, arma::Mat<double>& P) {
  mod->P_I = P;
};

//' @title Get infected pupae for deterministic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_PI_infection_deterministic(Rcpp::XPtr<culex_infection_deterministic> mod) {
  return mod->P_I;
};
