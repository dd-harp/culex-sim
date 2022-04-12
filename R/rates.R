# --------------------------------------------------------------------------------
#   external forcing
# --------------------------------------------------------------------------------

#' @title Temperature
#' @description A modified cosine function.
#' @param t time
#' @param pars a [list]
#' @export
temperature <- function(t, pars){
  
  phi = pars[["phi"]] # PHASE
  lambda = pars[["lambda"]] # A
  mu = pars[["mu"]] # M
  gamma = pars[["gamma"]] # POWER
  
  temp = 0.0
  
  if(t < 0.0){
    temp = (mu - lambda) + lambda * 2.0 * (0.5 * (1.0 + cos(2.0 * pi * (0.0 - phi) / 365.0)))^gamma
  } else {
    temp = (mu - lambda) + lambda * 2.0 * (0.5 * (1.0 + cos(2.0 * pi * (t - phi) / 365.0)))^gamma
  }
  
  return(temp)
}


# --------------------------------------------------------------------------------
#   lifecycle stage progression rates
# --------------------------------------------------------------------------------

#' @title Egg maturation rate
#' @param temp temperature
#' @param pars a list
#' @export
egg_maturation_rate <- function(temp, pars){
  
  alpha_E = pars[["alpha_E"]] # ALPHA
  beta_E = pars[["beta_E"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate egg development rate
  if(temp < 0.0){
    egg_maturation = 0.016667
  } else {
    egg_maturation = alpha_E * (temp^beta_E)
  }
  
  if(egg_maturation < maturation_min){
    egg_maturation = maturation_min
  }
  return(egg_maturation)
}

#' @title Larvae maturation rate
#' @param temp temperature
#' @param pars a list
#' @export
larvae_maturation_rate <- function(temp, pars){
  
  alpha_L = pars[["alpha_L"]] # ALPHA
  beta_L = pars[["beta_L"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate larvae development rate
  if(temp < 0.0){
    larvae_maturation = 0.016667
  } else {
    larvae_maturation = alpha_L * (temp^beta_L)
  }
  
  if(larvae_maturation < maturation_min){
    larvae_maturation = maturation_min
  }
  return(larvae_maturation)
}

#' @title Pupae maturation rate
#' @param temp temperature
#' @param pars a list
#' @export
pupae_maturation_rate <- function(temp, pars){
  
  alpha_P = pars[["alpha_P"]] # ALPHA
  beta_P = pars[["beta_P"]] # BETA
  maturation_min = pars[["maturation_min"]]
  
  # calculate larvae development rate
  if(temp < 0.0){
    pupae_maturation = 0.016667
  } else {
    pupae_maturation = alpha_P * (temp^beta_P)
  }
  
  if(pupae_maturation < maturation_min){
    pupae_maturation = maturation_min
  }
  return(pupae_maturation)
}


# --------------------------------------------------------------------------------
#   extrinsic incubation period
# --------------------------------------------------------------------------------
  
#' @title Extrinsic incubation development rate
#' @param temp temperature
#' @param pars a list
#' @export
eip_rate <- function(temp, pars) {
    
  q = pars[["eip_q"]]
  Tmax = pars[["eip_tmax"]]
  Tmin = pars[["eip_tmin"]]
  eip_min = pars[["eip_min"]]
    
  if (temp > Tmax) {
    return(eip_min)
  } else if (temp < Tmin) {
    return(eip_min)
  } else {
    eip = q*temp*(temp - Tmin) * sqrt(Tmax - temp)
    return(max(eip, eip_min))
  }
    
}

