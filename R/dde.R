# translation of culex model by David Ewing (https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90)

# --------------------------------------------------------------------------------
#   state variable history
# --------------------------------------------------------------------------------

#' @title Calculate history for DDE system
#' @param t temperature
#' @param p a [list]
#' @export
calculate_history <- function(t, p){
  
  temp <- temperature(t,p)
  
  # history vector
  Y <- rep(0,13)
  
  Y[8] = 1.0 / egg_maturation_rate(temp,p) # tau_E
  Y[9] = 1.0 / larvae_maturation_rate(temp,p) # tau_L
  Y[10] = 1.0 / pupae_maturation_rate(temp,p) # tau_P
  
  Y[5] = exp(-death_egg_rate(temp,p) * Y[8]) # S_E
  Y[6] = exp(-death_larvae_rate(temp,p) * Y[9]) # S_L
  Y[7] = exp(-death_pupae_rate(temp,p) * Y[10]) #S_P
  
  temp_L = temperature(t - Y[9],p)
  temp_P = temperature(t - Y[10],p)
  
  Y[11] = 1.0 / egg_maturation_rate(temp_L,p) # tau_E(t - tau_L(t))
  Y[12] = 1.0 / larvae_maturation_rate(temp_P,p) # tau_L(t - tau_P(t))
  
  temp_LP = temperature(t - Y[10] - Y[12],p)
  
  Y[13] = 1.0 / egg_maturation_rate(temp_LP,p) # tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))
  
  return(Y)
}


# --------------------------------------------------------------------------------
#   initial conditions
# --------------------------------------------------------------------------------

#' @title Calculate initial conditions for DDE system
#' @param A0 initial number of adults
#' @param t0 initial temperature (assumed constant for `t<0`)
#' @param pars a [list]
#' @export
calculate_IC <- function(A0, t0, pars){
  
  u0 <- rep(0,13)
  u0[4] = A0
  
  # calculate initial lags first
  u0[8] = 1.0 / egg_maturation_rate(t0,pars) # tau_E
  u0[9] = 1.0 / larvae_maturation_rate(t0,pars) # tau_L
  u0[10] = 1.0 / pupae_maturation_rate(t0,pars) # tau_P
  
  u0[11] = u0[8] # tau_E(t - tau_L(t))
  u0[12] = u0[9] # tau_L(t - tau_P(t))
  u0[13] = u0[8] # tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))
  
  # survival probabilities
  u0[5] = exp(-u0[8]*death_egg_rate(t0,pars)) # S_E
  u0[6] = exp(-u0[9]*death_larvae_rate(t0,pars)) # S_L
  u0[7] = exp(-u0[10]*death_pupae_rate(t0,pars)) # S_P
  
  return(u0)
}



# --------------------------------------------------------------------------------
#   system of DDEs
#   13 equations
#   6 delays
# --------------------------------------------------------------------------------

# dxdt
# Y(1) = E
# Y(2) = L
# Y(3) = P
# Y(4) = A
# Y(5) = S_E
# Y(6) = S_L
# Y(7) = S_P
# Y(8) = tau_E
# Y(9) = tau_L
# Y(10) = tau_P
# Y(11) = tau_E(t - tau_L(t))
# Y(12) = tau_L(t - tau_P(t))
# Y(13) = tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))

# delays
# Z(x,1) = t - tau_E(t)
# Z(x,2) = t - tau_L(t) - tau_E(t - tau_L(t))
# Z(x,3) = t - tau_P(t) - tau_L(t - tau_P(t)) - tau_E(t - tau_P(t) - tau_L(t - tau_P(t))
# Z(x,4) = t - tau_L(t)
# Z(x,5) = t - tau_P(t) - tau_L(t - tau_P(t))
# Z(x,6) = t - tau_P(t)

# when you access values of Z in FORTRAN, corresponds to h in Julia
# Z(4,1) = A(t - tau_E(t))
# Z(4,2) = A(t - tau_L(t) - tau_E(t - tau_L(t)))
# Z(5,4) = S_E(t - tau_L(t))
# Z(4,3) = A(t - tau_P(t) - tau_L(t - tau_P(t)) - tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))
# Z(5,5) = S_E(t - tau_P(t) - tau_L(t - tau_P(t)))
# Z(6,6) = S_L(t - tau_P(t))
# Z(2,4) = L(t - tau_L(t))

#' @title DDE system
#' @description R translation of culex model by David Ewing [(code here)](https://github.com/davewi13/Temperate-Mosquito-DDE/blob/master/Chapter%202%20DDE%20code.f90)
#' @param t time
#' @param y state
#' @param params a [list]
#' @importFrom PBSddesolve pastvalue
#' @export
culex_dde <- function(t,y,params){
  
  # state variables
  E <- y[1]
  LAR <- y[2]
  PUP <- y[3]
  ADU <- y[4]
  SE <- y[5]
  SL <- y[6]
  SP <- y[7]
  DE <- y[8] # tau_E(t)
  DL <- y[9] # tau_L(t)
  DP <- y[10] # tau_P(t)
  DEL <- y[11] # tau_E(t - tau_L(t))
  DLP <- y[12] # tau_L(t - tau_P(t))
  DELP <- y[13] # tau_E(t - tau_P(t) - tau_L(t - taup_P(t)))
  
  # larval predation parameters
  p0 <- params[["p0"]]
  p1 <- params[["p1"]]
  
  # Z: state variables at each of the 6 lagged times (lags follow same order as Z/BETA in DDE_SOLVER)
  
  # Z(x,1): t - tau_E(t)
  d_z1 <- t - DE
  if(d_z1 > 0){
    Z1 <- pastvalue(d_z1)
  } else {
    Z1 <- calculate_history(d_z1, params)
  }
  
  # Z(x,2): t - tau_L(t) - tau_E(t - tau_L(t))
  d_z2 <- t - DL - DEL
  if(d_z2 > 0){
    Z2 <- pastvalue(d_z2)
  } else {
    Z2 <- calculate_history(d_z2, params)
  }
  
  # Z(x,3): t - tau_P(t) - tau_L(t - tau_P(t)) - tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))
  d_z3 <- t - DP - DLP - DELP
  if(d_z3 > 0){
    Z3 <- pastvalue(d_z3)
  } else {
    Z3 <- calculate_history(d_z3, params)
  }
  
  # Z(x,4): t - tau_L(t)
  d_z4 <- t - DL
  if(d_z4 > 0){
    Z4 <- pastvalue(d_z4)
  } else {
    Z4 <- calculate_history(d_z4, params)
  }
  
  # Z(x,5): t - tau_P(t) - tau_L(t - tau_P(t))
  d_z5 <- t - DP - DLP
  if(d_z5 > 0){
    Z5 <- pastvalue(d_z5)
  } else {
    Z5 <- calculate_history(d_z5, params)
  }
  
  # Z(x,6): t - tau_P(t)
  d_z6 <- t - DP
  if(d_z6 > 0){
    Z6 <- pastvalue(d_z6)
  } else {
    Z6 <- calculate_history(d_z6, params)
  }
  
  # (lagged) temperature
  temp <- temperature(t, params)
  temp_E <- temperature(t - DE, params)
  temp_L <- temperature(t - DL, params)
  temp_P <- temperature(t - DP, params)
  temp_EL <- temperature(t - DL - Z4[8], params)
  temp_ELP <- temperature(t - DP - Z6[9] - Z5[8], params)
  temp_LP <- temperature(t - DP - Z6[9], params)
  
  
  # (lagged) photoperiod
  pp <- photoperiod(t, params)
  pp_1 <- photoperiod(t - 1, params)
  pp_E <- photoperiod(t - DE, params)
  pp_EL <- photoperiod(t - DL - Z4[8], params)
  pp_ELP <- photoperiod(t - DP - Z6[9] - Z5[8], params)
  
  # (lagged) gonotrophic cycle
  gon <- gonotrophic(temp, params)
  gon_E <- gonotrophic(temp_E, params)
  gon_EL <- gonotrophic(temp_EL, params)
  gon_ELP <- gonotrophic(temp_ELP, params)
  
  # diapause and birth
  if(pp > pp_1){
    dia <- diapause_spring(pp)
    dia_E <- diapause_spring(pp_E)
    dia_EL <- diapause_spring(pp_EL)
    dia_ELP <- diapause_spring(pp_ELP)
  } else {
    dia <- diapause_autumn(pp)
    dia_E <- diapause_autumn(pp_E)
    dia_EL <- diapause_autumn(pp_EL)
    dia_ELP <- diapause_autumn(pp_ELP)
  }
  
  birth <- oviposition(dia,gon,params)
  birth_E <- oviposition(dia_E,gon_E,params)
  birth_EL <- oviposition(dia_EL,gon_EL,params)
  birth_ELP <- oviposition(dia_ELP,gon_ELP,params)
  
  # (lagged) death
  death_egg <- death_egg_rate(temp,params)
  death_egg_E <- death_egg_rate(temp_E,params)
  
  death_larvae <- death_larvae_rate(temp,params)
  death_larvae_L <- death_larvae_rate(temp_L,params)
  
  death_pupae <- death_pupae_rate(temp,params)
  death_pupae_P <- death_pupae_rate(temp_P,params)
  
  death_adult <- death_adult_rate(temp,params)
  
  # (lagged) development
  larvae_maturation <- larvae_maturation_rate(temp,params)
  larvae_maturation_L <- larvae_maturation_rate(temp_L,params)
  larvae_maturation_P <- larvae_maturation_rate(temp_P,params)
  larvae_maturation_LP <- larvae_maturation_rate(temp_LP,params)
  
  egg_maturation <- egg_maturation_rate(temp,params)
  egg_maturation_E <- egg_maturation_rate(temp_E,params)
  egg_maturation_L <- egg_maturation_rate(temp_L,params)
  egg_maturation_EL <- egg_maturation_rate(temp_EL,params)
  egg_maturation_LP <- egg_maturation_rate(temp_LP,params)
  egg_maturation_ELP <- egg_maturation_rate(temp_ELP,params)
  
  pupae_maturation <- pupae_maturation_rate(temp,params)
  pupae_maturation_P <- pupae_maturation_rate(temp_P,params)
  
  # DDEs describing change in state duration
  dDEdt = 1 - egg_maturation/egg_maturation_E
  dDLdt = 1 - larvae_maturation/larvae_maturation_L
  dDPdt = 1 - pupae_maturation/pupae_maturation_P
  dDELdt = (1 - dDLdt) * (1 - egg_maturation_L/egg_maturation_EL)
  dDLPdt = (1 - dDPdt) * (1 - larvae_maturation_P/larvae_maturation_LP)
  dDELPdt = (1 - dDPdt - dDLPdt) * (1 - egg_maturation_LP/egg_maturation_ELP)
  
  # stage recruitment
  R_E = birth * ADU
  R_L = birth_E * Z1[4] * SE * egg_maturation/egg_maturation_E
  R_P = birth_EL * Z2[4] * Z4[5] * SL * larvae_maturation/larvae_maturation_L * (1 - dDELdt)
  R_A = birth_ELP * Z3[4] * Z5[5] * Z6[6] * SP * pupae_maturation/pupae_maturation_P * (1 - dDLPdt) * (1 - dDELPdt)
  
  # maturation rates
  M_E = R_L
  M_L = R_P
  M_P = R_A
  
  # death rates
  D_E = death_egg * E
  D_L = ((p0*LAR/(p1+LAR)) + death_larvae) * LAR
  D_P = death_pupae * PUP
  D_A = death_adult * ADU
  
  # DDE system
  du <- rep(NaN,13)
  du[1] = R_E - M_E - D_E  # E
  du[2] = R_L - M_L - D_L  # L
  du[3] = R_P - M_P - D_P  # P
  du[4] = R_A - D_A        # A
  
  du[5] = SE * ((egg_maturation * death_egg_E / egg_maturation_E) - death_egg)
  du[6] = SL * (((p0*Z4[2] / (p1+Z4[2])) + death_larvae_L) * (1-dDLdt) - (p0*LAR / (p1+LAR)) - death_larvae)
  du[7] = SP * ((pupae_maturation * death_pupae_P / pupae_maturation_P) - death_pupae)
  
  du[8] = dDEdt # tau_E(t)
  du[9] = dDLdt # tau_L(t)
  du[10] = dDPdt # tau_P(t)
  du[11] = dDELdt # tau_E(t - tau_L(t))
  du[12] = dDLPdt # tau_L(t - tau_P(t))
  du[13] = dDELPdt # tau_E(t - tau_P(t) - tau_L(t - tau_P(t)))
  
  return(du)
}
