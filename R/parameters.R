#' @title Model parameters
#' @description Most of these parameters can be found at [https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90].
#' @param phi temperature phase
#' @param lambda temperature A
#' @param mu temperature M
#' @param gamma temperature power
#' @param L latitude for photoperiod
#' @param max_egg maximum egg raft size R
#' @param q1 gonotrophic cycle KG
#' @param q2 gonotrophic cycle QG
#' @param q3 gonotrophic cycle BG
#' @param nu_0E egg stage mortality U3
#' @param nu_1E egg stage mortality U4
#' @param nu_2E egg stage mortality U5
#' @param nu_0L larval stage mortality U3
#' @param nu_1L larval stage mortality U4
#' @param nu_2L larval stage mortality U5
#' @param nu_0P pupal stage mortality U3
#' @param nu_1P pupal stage mortality U4
#' @param nu_2P pupal stage mortality U5
#' @param alpha_A adult stage mortality ALPHA
#' @param beta_A adult stage mortality BETA
#' @param alpha_E egg stage development ALPHA
#' @param beta_E egg stage development BETA
#' @param alpha_L larval stage development ALPHA
#' @param beta_L larval stage development BETA
#' @param alpha_P pupal stage development ALPHA
#' @param beta_P pupal stage development BETA
#' @param a larval predation
#' @param h larval predation
#' @param r larval predation
#' @param V larval predation
#' @param death_max maximum death rate for aquatic stages
#' @param death_min_a minimum death rate for adult stage
#' @param gon_min minimum rate of gonotrophic cycle
#' @param maturation_min minimum rate of aquatic stage development
#' @param eip_q EIP parameter
#' @param eip_tmax EIP parameter
#' @param eip_tmin EIP parameter
#' @param eip_min minimum rate of EIP progression
#' @export
culex_parameters <- function(
  phi = 1.4,
  lambda = 6.3,
  mu = 10.3,
  gamma = 1.21,
  L = 51,
  max_egg = 200,
  q1 = 0.2024,
  q2 = 74.48,
  q3 = 0.2456,
  nu_0E = 0.0157,
  nu_1E = 20.5,
  nu_2E = 7,
  nu_0L = 0.0157,
  nu_1L = 20.5,
  nu_2L = 7,
  nu_0P = 0.0157,
  nu_1P = 20.5,
  nu_2P = 7,
  alpha_A = 2.166e-8,
  beta_A = 4.483,
  alpha_E = 0.0022,
  beta_E = 1.77,
  alpha_L = 0.00315,
  beta_L = 1.12,
  alpha_P = 0.0007109,
  beta_P = 1.8865648,
  a = 1,
  h = 0.002,
  r = 0.001,
  V = 200,
  death_max = 1.0,
  death_min_a = 0.01,
  gon_min = 0.0333,
  maturation_min = 0.016667,
  eip_q = 7.83e-5,
  eip_tmax = 45.2,
  eip_tmin = 11.4,
  eip_min = 0.005
) {
  pars <- list(
    phi = phi,
    lambda = lambda,
    mu = mu,
    gamma = gamma,
    L = L,
    max_egg = max_egg,
    q1 = q1,
    q2 = q2,
    q3 = q3,
    nu_0E = nu_0E,
    nu_1E = nu_1E,
    nu_2E = nu_2E,
    nu_0L = nu_0L,
    nu_1L = nu_1L,
    nu_2L = nu_2L,
    nu_0P = nu_0P,
    nu_1P = nu_1P,
    nu_2P = nu_2P,
    alpha_A = alpha_A,
    beta_A = beta_A,
    alpha_E = alpha_E,
    beta_E = beta_E,
    alpha_L = alpha_L,
    beta_L = beta_L,
    alpha_P = alpha_P,
    beta_P = beta_P,
    a = a,
    h = h,
    r = r,
    V = V,
    death_max = death_max,
    death_min_a = death_min_a,
    gon_min = gon_min,
    maturation_min = maturation_min,
    eip_q = eip_q,
    eip_tmax = eip_tmax,
    eip_tmin = eip_tmin,
    eip_min = eip_min
  )
  pars[["p0"]] <- pars[["r"]] / pars[["h"]]
  pars[["p1"]] <- pars[["V"]] / (pars[["a"]] * pars[["h"]])
  return(pars)
}