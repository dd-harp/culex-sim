test_that("stochastic model object is created correctly", {
  parameters <- culex_parameters()
  p <- 2
  tau_E <- c(5,6,7,6,5)
  tau_L <- c(2,3,2,3,3)
  tau_P <- c(4,5,5,5,4)
  tau_EIP <- c(6,5,4,6,7)
  dt <- 0.5
  psi <- diag(2)
  n_species <- 2
  mod <- create_culex_infection_stochastic(p = p, tau_E = tau_E, tau_L = tau_L, tau_P = tau_P, tau_EIP = tau_EIP, dt = dt, psi = psi, n_species = n_species)
  
  expect_true(all(get_AS_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AS_infection_stochastic(mod)) == 1L)
  expect_true(ncol(get_AS_infection_stochastic(mod)) == p)
  
  expect_true(all(get_AE_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AE_infection_stochastic(mod)) == max(tau_EIP))
  expect_true(ncol(get_AE_infection_stochastic(mod)) == p)
  
  expect_true(all(get_AI_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AI_infection_stochastic(mod)) == 1L)
  expect_true(ncol(get_AI_infection_stochastic(mod)) == p)
  
})


test_that("stochastic model object has correct stage advancement", {
  parameters <- culex_parameters()
  p <- 2
  tau_E <- c(5,6,7,6,5)
  tau_L <- c(2,3,2,3,3)
  tau_P <- c(4,5,5,5,4)
  tau_EIP <- c(6,5,4,6,7)
  dt <- 0.5
  psi <- diag(2)
  n_species <- 2
  
  # check AE advances through EIP properly
  mod <- create_culex_infection_stochastic(p = p, tau_E = tau_E, tau_L = tau_L, tau_P = tau_P, tau_EIP = tau_EIP, dt = dt, psi = psi, n_species = n_species)
  
  newE <- matrix(0, nrow = max(tau_EIP), ncol = p)
  newE[max(tau_EIP), ] <- 1e3
  set_AE_infection_stochastic(mod, newE)
  
  step_culex_infection_stochastic(mod = mod, parameters = parameters)
  outE <- get_AE_infection_stochastic(mod)
  
  expect_true(all(outE[max(tau_EIP) -1 ,] > 0))
  expect_true(sum(outE) <= sum(newE))
  expect_true(all(outE[-(max(tau_EIP) -1), ] == 0))
  
  # check AE goes to AI at end of their incubation period
  mod <- create_culex_infection_stochastic(p = p, tau_E = tau_E, tau_L = tau_L, tau_P = tau_P, tau_EIP = tau_EIP, dt = dt, psi = psi, n_species = n_species)
  
  newE <- newE*0
  newE[1, ] <- 1e3
  set_AE_infection_stochastic(mod, newE)
  
  step_culex_infection_stochastic(mod = mod, parameters = parameters)
  
  outE <- get_AE_infection_stochastic(mod)
  outI <- get_AI_infection_stochastic(mod)
  
  expect_true(sum(outE) == 0)
  expect_true(sum(outI) <= sum(newE))
  
  # check new infecteds go to the right spot in EIP
  mod <- create_culex_infection_stochastic(p = p, tau_E = tau_E, tau_L = tau_L, tau_P = tau_P, tau_EIP = tau_EIP, dt = dt, psi = psi, n_species = n_species)
  
  set_AS_infection_stochastic(mod, c(1e3, 1e3))
  set_f_infection_stochastic(mod, c(5, 5))
  set_q_infection_stochastic(mod, matrix(0.5, 2, 2))
  set_kappa_infection_stochastic(mod, matrix(0.5, 2, 2))
  
  step_culex_infection_stochastic(mod = mod, parameters = parameters)
  outAS <- get_AS_infection_stochastic(mod)
  outAE <- get_AE_infection_stochastic(mod)
  
  expect_true(sum(outAS) < 2e3)
  expect_true(sum(outAE[-tau_EIP[1], ]) == 0)
  expect_true(sum(outAE) > 0)
  
  # check that infected pupae go to I adults
  mod <- create_culex_infection_stochastic(p = p, tau_E = tau_E, tau_L = tau_L, tau_P = tau_P, tau_EIP = tau_EIP, dt = dt, psi = psi, n_species = n_species)
  
  newP <- matrix(0, nrow = max(tau_P), ncol = p)
  newP[1, ] <- 1e3
  set_PI_infection_stochastic(mod, newP)
  set_AI_infection_stochastic(mod, c(0, 0))
  set_AS_infection_stochastic(mod, c(0, 0))
  
  step_culex_infection_stochastic(mod = mod, parameters = parameters)
  
  expect_true(sum(get_PI_infection_stochastic(mod)) == 0)
  expect_true(sum(get_AI_infection_stochastic(mod)) <= sum(newP))
  expect_true(sum(get_AS_infection_stochastic(mod)) == 0)
})
