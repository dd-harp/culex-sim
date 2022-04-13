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
  
  expect_true(all(get_AS_culex_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AS_culex_infection_stochastic(mod)) == 1L)
  expect_true(ncol(get_AS_culex_infection_stochastic(mod)) == p)
  
  expect_true(all(get_AE_culex_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AE_culex_infection_stochastic(mod)) == max(tau_EIP))
  expect_true(ncol(get_AE_culex_infection_stochastic(mod)) == p)
  
  expect_true(all(get_AI_culex_infection_stochastic(mod) == 0))
  expect_true(nrow(get_AI_culex_infection_stochastic(mod)) == 1L)
  expect_true(ncol(get_AI_culex_infection_stochastic(mod)) == p)
  
  newE <- matrix(0, nrow = max(tau_EIP), ncol = p)
  newE[max(tau_EIP), ] <- 1e3
  set_AE_culex_infection_stochastic(mod, newE)
  # testAE <- newE*0
  # testAE[max(tau_EIP)-1, ] <- 10
  
  get_AE_culex_infection_stochastic(mod)
  step_culex_infection_stochastic(mod = mod, parameters = parameters)
  get_AE_culex_infection_stochastic(mod)
  
  expect_true()
})
