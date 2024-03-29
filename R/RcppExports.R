# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Create deterministic culex model object
#' @param p number of patches
#' @param tau_E vector of egg delays
#' @param tau_L vector of larvae delays
#' @param tau_P vector of pupae delays
#' @param dt size of time step
#' @param parameters a [list] of parameters
#' @export
create_culex_deterministic <- function(tau_E, tau_L, tau_P, dt, parameters) {
    .Call('_culex_create_culex_deterministic', PACKAGE = 'culex', tau_E, tau_L, tau_P, dt, parameters)
}

#' @title Step deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @param parameters a [list] of parameters
#' @export
step_culex_deterministic <- function(mod, parameters) {
    invisible(.Call('_culex_step_culex_deterministic', PACKAGE = 'culex', mod, parameters))
}

#' @title Set adults for deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @param A a vector
#' @export
set_A_deterministic <- function(mod, A) {
    invisible(.Call('_culex_set_A_deterministic', PACKAGE = 'culex', mod, A))
}

#' @title Return adults for deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_A_deterministic <- function(mod) {
    .Call('_culex_get_A_deterministic', PACKAGE = 'culex', mod)
}

#' @title Return eggs for deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_E_deterministic <- function(mod) {
    .Call('_culex_get_E_deterministic', PACKAGE = 'culex', mod)
}

#' @title Return larvae for deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_L_deterministic <- function(mod) {
    .Call('_culex_get_L_deterministic', PACKAGE = 'culex', mod)
}

#' @title Return pupae for deterministic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_P_deterministic <- function(mod) {
    .Call('_culex_get_P_deterministic', PACKAGE = 'culex', mod)
}

#' @title Create deterministic culex infection model object
#' @param tau_E vector of egg delays
#' @param tau_L vector of larvae delays
#' @param tau_P vector of pupae delays
#' @param tau_EIP vector of extrinsic incubation period delays
#' @param dt size of time step
#' @param parameters a [list] of parameters
#' @param n_species number of host species
#' @export
create_culex_infection_deterministic <- function(tau_E, tau_L, tau_P, tau_EIP, dt, parameters, n_species) {
    .Call('_culex_create_culex_infection_deterministic', PACKAGE = 'culex', tau_E, tau_L, tau_P, tau_EIP, dt, parameters, n_species)
}

#' @title Step deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param parameters a [list] of parameters
#' @export
step_culex_infection_deterministic <- function(mod, parameters) {
    invisible(.Call('_culex_step_culex_infection_deterministic', PACKAGE = 'culex', mod, parameters))
}

#' @title Set feeding rate for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param f a row vector (feeding rate by patch)
#' @export
set_f_infection_deterministic <- function(mod, f) {
    invisible(.Call('_culex_set_f_infection_deterministic', PACKAGE = 'culex', mod, f))
}

#' @title Get feeding rate for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_f_infection_deterministic <- function(mod) {
    .Call('_culex_get_f_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set feeding habit for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param q a matrix (columns must sum to `1`; each column indicates proportion of bites allocated to each host species in that patch)
#' @export
set_q_infection_deterministic <- function(mod, q) {
    invisible(.Call('_culex_set_q_infection_deterministic', PACKAGE = 'culex', mod, q))
}

#' @title Get feeding habit for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_q_infection_deterministic <- function(mod) {
    .Call('_culex_get_q_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set kappa for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param kappa a matrix (each column indicates net infectiousness of each host species in that patch)
#' @export
set_kappa_infection_deterministic <- function(mod, kappa) {
    invisible(.Call('_culex_set_kappa_infection_deterministic', PACKAGE = 'culex', mod, kappa))
}

#' @title Get kappa for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_kappa_infection_deterministic <- function(mod) {
    .Call('_culex_get_kappa_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set susceptible adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a row vector
#' @export
set_AS_infection_deterministic <- function(mod, A) {
    invisible(.Call('_culex_set_AS_infection_deterministic', PACKAGE = 'culex', mod, A))
}

#' @title Get susceptible adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AS_infection_deterministic <- function(mod) {
    .Call('_culex_get_AS_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set incubating adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a matrix
#' @export
set_AE_infection_deterministic <- function(mod, A) {
    invisible(.Call('_culex_set_AE_infection_deterministic', PACKAGE = 'culex', mod, A))
}

#' @title Get incubating adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AE_infection_deterministic <- function(mod) {
    .Call('_culex_get_AE_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set infectious adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a row vector
#' @export
set_AI_infection_deterministic <- function(mod, A) {
    invisible(.Call('_culex_set_AI_infection_deterministic', PACKAGE = 'culex', mod, A))
}

#' @title Get infectious adults for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AI_infection_deterministic <- function(mod) {
    .Call('_culex_get_AI_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set eggs for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param E a matrix
#' @export
set_E_infection_deterministic <- function(mod, E) {
    invisible(.Call('_culex_set_E_infection_deterministic', PACKAGE = 'culex', mod, E))
}

#' @title Get eggs for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_E_infection_deterministic <- function(mod) {
    .Call('_culex_get_E_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set infected eggs for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param E a matrix
#' @export
set_EI_infection_deterministic <- function(mod, E) {
    invisible(.Call('_culex_set_EI_infection_deterministic', PACKAGE = 'culex', mod, E))
}

#' @title Get infected eggs for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_EI_infection_deterministic <- function(mod) {
    .Call('_culex_get_EI_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set larvae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param L a matrix
#' @export
set_L_infection_deterministic <- function(mod, L) {
    invisible(.Call('_culex_set_L_infection_deterministic', PACKAGE = 'culex', mod, L))
}

#' @title Get larvae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_L_infection_deterministic <- function(mod) {
    .Call('_culex_get_L_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set infected larvae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param L a matrix
#' @export
set_LI_infection_deterministic <- function(mod, L) {
    invisible(.Call('_culex_set_LI_infection_deterministic', PACKAGE = 'culex', mod, L))
}

#' @title Get infected larvae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_LI_infection_deterministic <- function(mod) {
    .Call('_culex_get_LI_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set pupae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param P a matrix
#' @export
set_P_infection_deterministic <- function(mod, P) {
    invisible(.Call('_culex_set_P_infection_deterministic', PACKAGE = 'culex', mod, P))
}

#' @title Get pupae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_P_infection_deterministic <- function(mod) {
    .Call('_culex_get_P_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Set infected pupae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param P a matrix
#' @export
set_PI_infection_deterministic <- function(mod, P) {
    invisible(.Call('_culex_set_PI_infection_deterministic', PACKAGE = 'culex', mod, P))
}

#' @title Get infected pupae for deterministic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_PI_infection_deterministic <- function(mod) {
    .Call('_culex_get_PI_infection_deterministic', PACKAGE = 'culex', mod)
}

#' @title Create stochastic culex infection model object
#' @param tau_E vector of egg delays
#' @param tau_L vector of larvae delays
#' @param tau_P vector of pupae delays
#' @param tau_EIP vector of extrinsic incubation period delays
#' @param dt size of time step
#' @param parameters a [list] of parameters
#' @param n_species number of host species
#' @export
create_culex_infection_stochastic <- function(tau_E, tau_L, tau_P, tau_EIP, dt, parameters, n_species) {
    .Call('_culex_create_culex_infection_stochastic', PACKAGE = 'culex', tau_E, tau_L, tau_P, tau_EIP, dt, parameters, n_species)
}

#' @title Step stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param parameters a [list] of parameters
#' @export
step_culex_infection_stochastic <- function(mod, parameters) {
    invisible(.Call('_culex_step_culex_infection_stochastic', PACKAGE = 'culex', mod, parameters))
}

#' @title Set feeding rate for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param f a row vector (feeding rate by patch)
#' @export
set_f_infection_stochastic <- function(mod, f) {
    invisible(.Call('_culex_set_f_infection_stochastic', PACKAGE = 'culex', mod, f))
}

#' @title Get feeding rate for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_f_infection_stochastic <- function(mod) {
    .Call('_culex_get_f_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set feeding habit for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param q a matrix (columns must sum to `1`; each column indicates proportion of bites allocated to each host species in that patch)
#' @export
set_q_infection_stochastic <- function(mod, q) {
    invisible(.Call('_culex_set_q_infection_stochastic', PACKAGE = 'culex', mod, q))
}

#' @title Get feeding habit for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_q_infection_stochastic <- function(mod) {
    .Call('_culex_get_q_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set kappa for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param kappa a matrix (each column indicates net infectiousness of each host species in that patch)
#' @export
set_kappa_infection_stochastic <- function(mod, kappa) {
    invisible(.Call('_culex_set_kappa_infection_stochastic', PACKAGE = 'culex', mod, kappa))
}

#' @title Get kappa for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_kappa_infection_stochastic <- function(mod) {
    .Call('_culex_get_kappa_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set susceptible adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a row vector
#' @export
set_AS_infection_stochastic <- function(mod, A) {
    invisible(.Call('_culex_set_AS_infection_stochastic', PACKAGE = 'culex', mod, A))
}

#' @title Get susceptible adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AS_infection_stochastic <- function(mod) {
    .Call('_culex_get_AS_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set incubating adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a matrix
#' @export
set_AE_infection_stochastic <- function(mod, A) {
    invisible(.Call('_culex_set_AE_infection_stochastic', PACKAGE = 'culex', mod, A))
}

#' @title Get incubating adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AE_infection_stochastic <- function(mod) {
    .Call('_culex_get_AE_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set infectious adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param A a row vector
#' @export
set_AI_infection_stochastic <- function(mod, A) {
    invisible(.Call('_culex_set_AI_infection_stochastic', PACKAGE = 'culex', mod, A))
}

#' @title Get infectious adults for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_AI_infection_stochastic <- function(mod) {
    .Call('_culex_get_AI_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set eggs for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param E a matrix
#' @export
set_E_infection_stochastic <- function(mod, E) {
    invisible(.Call('_culex_set_E_infection_stochastic', PACKAGE = 'culex', mod, E))
}

#' @title Get eggs for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_E_infection_stochastic <- function(mod) {
    .Call('_culex_get_E_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set infected eggs for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param E a matrix
#' @export
set_EI_infection_stochastic <- function(mod, E) {
    invisible(.Call('_culex_set_EI_infection_stochastic', PACKAGE = 'culex', mod, E))
}

#' @title Get infected eggs for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_EI_infection_stochastic <- function(mod) {
    .Call('_culex_get_EI_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set larvae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param L a matrix
#' @export
set_L_infection_stochastic <- function(mod, L) {
    invisible(.Call('_culex_set_L_infection_stochastic', PACKAGE = 'culex', mod, L))
}

#' @title Get larvae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_L_infection_stochastic <- function(mod) {
    .Call('_culex_get_L_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set infected larvae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param L a matrix
#' @export
set_LI_infection_stochastic <- function(mod, L) {
    invisible(.Call('_culex_set_LI_infection_stochastic', PACKAGE = 'culex', mod, L))
}

#' @title Get infected larvae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_LI_infection_stochastic <- function(mod) {
    .Call('_culex_get_LI_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set pupae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param P a matrix
#' @export
set_P_infection_stochastic <- function(mod, P) {
    invisible(.Call('_culex_set_P_infection_stochastic', PACKAGE = 'culex', mod, P))
}

#' @title Get pupae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_P_infection_stochastic <- function(mod) {
    .Call('_culex_get_P_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Set infected pupae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @param P a matrix
#' @export
set_PI_infection_stochastic <- function(mod, P) {
    invisible(.Call('_culex_set_PI_infection_stochastic', PACKAGE = 'culex', mod, P))
}

#' @title Get infected pupae for stochastic culex infection model
#' @param mod an [methods::externalptr-class] object
#' @export
get_PI_infection_stochastic <- function(mod) {
    .Call('_culex_get_PI_infection_stochastic', PACKAGE = 'culex', mod)
}

#' @title Create stochastic culex model object
#' @param p number of patches
#' @param tau_E vector of egg delays
#' @param tau_L vector of larvae delays
#' @param tau_P vector of pupae delays
#' @param dt size of time step
#' @param parameters a [list] of parameters
#' @export
create_culex_stochastic <- function(tau_E, tau_L, tau_P, dt, parameters) {
    .Call('_culex_create_culex_stochastic', PACKAGE = 'culex', tau_E, tau_L, tau_P, dt, parameters)
}

#' @title Step stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @param parameters a [list] of parameters
#' @export
step_culex_stochastic <- function(mod, parameters) {
    invisible(.Call('_culex_step_culex_stochastic', PACKAGE = 'culex', mod, parameters))
}

#' @title Set adults for stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @param A a vector
#' @export
set_A_stochastic <- function(mod, A) {
    invisible(.Call('_culex_set_A_stochastic', PACKAGE = 'culex', mod, A))
}

#' @title Return adults for stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_A_stochastic <- function(mod) {
    .Call('_culex_get_A_stochastic', PACKAGE = 'culex', mod)
}

#' @title Return eggs for stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_E_stochastic <- function(mod) {
    .Call('_culex_get_E_stochastic', PACKAGE = 'culex', mod)
}

#' @title Return larvae for stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_L_stochastic <- function(mod) {
    .Call('_culex_get_L_stochastic', PACKAGE = 'culex', mod)
}

#' @title Return pupae for stochastic culex model
#' @param mod an [methods::externalptr-class] object
#' @export
get_P_stochastic <- function(mod) {
    .Call('_culex_get_P_stochastic', PACKAGE = 'culex', mod)
}

#' @title Differential equations describing change in delays
#' @description This function is meant to be passed to [deSolve::ode].
#' @param t time
#' @param y state
#' @param params a list of parameters
#' @export
tau_diffeqn <- function(t, y, params) {
    .Call('_culex_tau_diffeqn', PACKAGE = 'culex', t, y, params)
}

