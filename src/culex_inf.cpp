#include "culex_types.h"

// ---------- stochastic model interface ----------

//' @title Create stochastic culex infection model object
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

//' @title Step stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param parameters a [list] of parameters
//' @export
// [[Rcpp::export]]
void step_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, const Rcpp::List& parameters) {
  mod->update(parameters);
}

//' @title Set feeding rate for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param f a row vector (feeding rate by patch)
//' @export
// [[Rcpp::export]]
void set_f_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, const arma::Row<double>& f) {
  mod->f = f;
}

//' @title Get feeding rate for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<double> get_f_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->f;
} 

//' @title Set feeding habit for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param q a matrix (columns must sum to `1`; each column indicates proportion of bites allocated to each host species in that patch)
//' @export
// [[Rcpp::export]]
void set_q_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, const arma::Mat<double>& q) {
  mod->q = q;
}

//' @title Get feeding habit for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_q_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->q;
} 

//' @title Set kappa for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param kappa a matrix (each column indicates net infectiousness of each host species in that patch)
//' @export
// [[Rcpp::export]]
void set_kappa_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, const arma::Mat<double>& kappa) {
  mod->kappa = kappa;
}

//' @title Get kappa for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<double> get_kappa_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->kappa;
} 

//' @title Set susceptible adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a row vector
//' @export
// [[Rcpp::export]]
void set_AS_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Row<int>& A) {
  mod->A_S = A;
};

//' @title Get susceptible adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<int> get_AS_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->A_S;
};

//' @title Set incubating adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a matrix
//' @export
// [[Rcpp::export]]
void set_AE_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& A) {
  mod->A_E = A;
};

//' @title Get incubating adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_AE_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->A_E;
};

//' @title Set infectious adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param A a row vector
//' @export
// [[Rcpp::export]]
void set_AI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Row<int>& A) {
  mod->A_I = A;
};

//' @title Get infectious adults for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Row<int> get_AI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->A_I;
};

//' @title Set eggs for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param E a matrix
//' @export
// [[Rcpp::export]]
void set_E_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& E) {
  mod->E = E;
};

//' @title Get eggs for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_E_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->E;
};

//' @title Set infected eggs for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param E a matrix
//' @export
// [[Rcpp::export]]
void set_EI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& E) {
  mod->E_I = E;
};

//' @title Get infected eggs for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_EI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->E_I;
};

//' @title Set larvae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param L a matrix
//' @export
// [[Rcpp::export]]
void set_L_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& L) {
  mod->L = L;
};

//' @title Get larvae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_L_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->L;
};

//' @title Set infected larvae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param L a matrix
//' @export
// [[Rcpp::export]]
void set_LI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& L) {
  mod->L_I = L;
};

//' @title Get infected larvae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_LI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->L_I;
};

//' @title Set pupae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param P a matrix
//' @export
// [[Rcpp::export]]
void set_P_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& P) {
  mod->P = P;
};

//' @title Get pupae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_P_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->P;
};

//' @title Set infected pupae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @param P a matrix
//' @export
// [[Rcpp::export]]
void set_PI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod, arma::Mat<int>& P) {
  mod->P_I = P;
};

//' @title Get infected pupae for stochastic culex infection model
//' @param mod an [methods::externalptr-class] object
//' @export
// [[Rcpp::export]]
arma::Mat<int> get_PI_culex_infection_stochastic(Rcpp::XPtr<culex_infection_stochastic> mod) {
  return mod->P_I;
};


// ---------- deterministic model interface ----------