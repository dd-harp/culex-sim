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

#' @title Photoperiod
#' @param t time
#' @param pars a [list]
#' @export
photoperiod <- function(t, pars){
  
  L = pars[["L"]] # latitude (51 in thesis)
  
  # define photoperiod values
  EPS = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (t - 3.5)))))
  NUM = sin(0.8333 * pi/ 180.0) + (sin(L * pi / 180.0) * sin(EPS))
  DEN = cos(L * pi / 180.0) * cos(EPS)
  DAYLIGHT = 24.0 - (24.0 / pi) * acos(NUM / DEN)
  
  return(DAYLIGHT)
}


# --------------------------------------------------------------------------------
#   diapause
# --------------------------------------------------------------------------------

#' @title Diapause (spring)
#' @param pp photoperiod
#' @export
diapause_spring <- function(pp){
  1.0 / (1.0 + exp(5.0 * (14.0 - pp)))
}

#' @title Diapause (autumn)
#' @param pp photoperiod
#' @export
diapause_autumn <- function(pp){
  1.0 / (1.0 + exp(5.0 * (13.0 - pp)))
}


# --------------------------------------------------------------------------------
#   mortality
# --------------------------------------------------------------------------------

#' @title Egg mortality rate
#' @param temp temperature
#' @param pars a [list]
#' @export
death_egg_rate <- function(temp, pars){
  
  nu_0E = pars[["nu_0E"]] # U3
  nu_1E = pars[["nu_1E"]] # U4
  nu_2E = pars[["nu_2E"]] # U5
  death_max <- pars[["death_max"]]
  
  # calculate egg death rate
  egg_d = nu_0E * exp(((temp - nu_1E) / nu_2E)^2)
  
  if(egg_d > death_max){
    egg_d = death_max
  }
  
  return(egg_d)
}

#' @title Larvae mortality rate
#' @param temp temperature
#' @param pars a [list]
#' @export
death_larvae_rate <- function(temp, pars){
  
  nu_0L = pars[["nu_0L"]] # U3
  nu_1L = pars[["nu_1L"]] # U4
  nu_2L = pars[["nu_2L"]] # U5
  death_max <- pars[["death_max"]]
  
  # calculate egg death rate
  larvae_d = nu_0L * exp(((temp - nu_1L) / nu_2L)^2)
  
  if(larvae_d > death_max){
    larvae_d = death_max
  }
  
  return(larvae_d)
}

#' @title Pupae mortality rate
#' @param temp temperature
#' @param pars a [list]
#' @export
death_pupae_rate <- function(temp, pars){
  
  nu_0P = pars[["nu_0P"]] # U3
  nu_1P = pars[["nu_1P"]] # U4
  nu_2P = pars[["nu_2P"]] # U5
  death_max <- pars[["death_max"]]
  
  # calculate egg death rate
  pupal_d = nu_0P * exp(((temp - nu_1P)/nu_2P)^2)
  
  if(pupal_d > death_max){
    pupal_d = death_max
  }
  return(pupal_d)
}

#' @title Adult mortality rate
#' @param temp temperature
#' @param pars a [list]
#' @export
death_adult_rate <- function(temp, pars){
  
  alpha_A = pars[["alpha_A"]] # ALPHA
  beta_A = pars[["beta_A"]] # BETA
  death_min_a = pars[["death_min_a"]]
  
  # calculate adult death rate
  adult_d = alpha_A * (temp^beta_A)
  
  if(adult_d < death_min_a){
    adult_d = death_min_a
  }
  return(adult_d)
}


# --------------------------------------------------------------------------------
#   lifecycle stage progression rates
# --------------------------------------------------------------------------------

#' @title Gonotrophic cycle length
#' @param temp temperature
#' @param pars a [list]
#' @export
gonotrophic <- function(temp, pars){
  
  q1 = pars[["q1"]] # KG
  q2 = pars[["q2"]] # QG
  q3 = pars[["q3"]] # BG
  gon_min = pars[["gon_min"]]
  
  # calculate gonotrophic cycle length
  if(temp < 0.0){
    grate = 0.0333
  } else {
    grate = q1 / (1 + q2*exp(-q3*temp))
  }
  
  if(grate < gon_min){
    grate = gon_min
  }
  
  return(1.0 / grate)
}

#' @title Egg laying rate
#' @param d diapause
#' @param G duration of gonotrophic cycle
#' @param pars a [list]
#' @export
oviposition <- function(d, G, pars){
  
  max_egg = pars[["max_egg"]]
  
  egg_raft = d * max_egg * 0.5
  ovi = egg_raft / G
  
  return(ovi)
}

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

