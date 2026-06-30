#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "Helper.h"

#define _USE_MATH_DEFINES

#ifndef _MAT_LIB_
#define _MAT_LIB_

double rho = 2.77e3; // [kg/m^3]
double E = 71e9;    // [Pa]
double nu = 0.33;

/// Lame Coefficients
double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
double mu = E / (2.0 * (1.0 + nu));

/// Elasticity Matrix [D]
double D[directions][directions] = {
    {lambda + 2 * mu, lambda, lambda, 0, 0, 0},
    {lambda, lambda + 2 * mu, lambda, 0, 0, 0},
    {lambda, lambda, lambda + 2 * mu, 0, 0, 0},
    {0, 0, 0, mu, 0, 0},
    {0, 0, 0, 0, mu, 0},
    {0, 0, 0, 0, 0, mu},
};

#endif _MAT_LIB_