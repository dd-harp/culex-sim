#ifndef CULEX_TYPES
#define CULEX_TYPES

#include "culex.h"
#include "culex_inf.h"

using culex_stochastic = culex<int>;
using culex_deterministic = culex<double>;
using culex_infection_stochastic = culex_inf<int>;
using culex_infection_deterministic = culex_inf<double>;

#endif