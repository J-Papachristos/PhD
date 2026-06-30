#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "Helper.h"

#define _USE_MATH_DEFINES

#ifndef _MAT_LIB_
#define _MAT_LIB_

/// Elasticity Matrix [D]
// double D[directions][directions] = {
//     {lambda + 2 * mu, lambda, lambda, 0, 0, 0},
//     {lambda, lambda + 2 * mu, lambda, 0, 0, 0},
//     {lambda, lambda, lambda + 2 * mu, 0, 0, 0},
//     {0, 0, 0, mu, 0, 0},
//     {0, 0, 0, 0, mu, 0},
//     {0, 0, 0, 0, 0, mu},
// };

enum materialType {
    linear_elastic,
    neo_hookean,
    matTypes
};

class Material {
  private:
    int type;

  public:
    /// Physical Properties
    double rho; // Density [kg/m^3]

    /// Linear Elastic
    double E;  // Young's Modulus [Pa]
    double nu; // Poisson's Ratio [-]

    // Lame Parameters
    double lambda, mu;

    /// Elasticity Matrix [D]
    double D[directions][directions];

    /// Methods

    /// @brief Sets the Physical Properties of the Material
    /// @param rho Density [kg/m^3]
    void setPhysical(double rho) {
        this->rho = rho;
    }

    /// @brief Sets the Linear Elastic Parameters and Calculates the Lame Parameters
    /// @param E Young's Modulus [Pa]
    /// @param nu Poisson's Ratio [-]
    void setLinearElastic(double E, double nu) {
        this->E = E;
        this->nu = nu;

        this->lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        this->mu = E / (2.0 * (1.0 + nu));

        for (int i = 0; i < directions; i++) {
            for (int j = 0; j < directions; j++) {
                (i < directions / 2 && j < directions / 2) ? this->D[i][j] = this->lambda : this->D[i][j] = 0;
            }
            (i < directions / 2) ? this->D[i][i] = this->lambda + 2 * this->mu : this->D[i][i] = this->mu;
        }
    }
};

#endif // !_MAT_LIB_