#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <stdarg.h>

#include "Helper.h"

#define _USE_MATH_DEFINES

#ifndef _MAT_LIB_
#define _MAT_LIB_

enum materialType {
    linear_elastic,
    neo_hookean,
    mooney_rivlin,
    ogden,
    yeoh,
    matTypes
};

class Material {
  private:
    int type; // Material Model Type

    /// @brief Sets the Material Model
    void setType(int type) {
        this->type = type;
    }

  public:
    /// Physical Properties
    double rho; // Density [kg/m^3]

    /// Linear Elastic
    double E;  // Young's Modulus [Pa]
    double nu; // Poisson's Ratio [-]

    // Lame Parameters
    double lambda, mu;

    /// Methods

    /// @brief Gets the Model Type of the Material
    /// @return Material Model Type
    int getType() {
        return this->type;
    }

    /// @brief Sets the Physical Properties of the Material
    /// @param rho Density [kg/m^3]
    void setPhysical(double rho) {
        this->rho = rho;
    }

    /// @brief Sets the Elasticity Parameters of the Material, depending on the Model Type
    /// @param type Model Type
    /// @param
    void setElasticity(int type, ...) {
        this->setType(type);
        va_list args;

        switch (type) {
            default:
            case linear_elastic:
                va_start(args, type);

                this->E = va_arg(args, double);
                this->nu = va_arg(args, double);
                va_end(args);

                this->lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
                this->mu = E / (2.0 * (1.0 + nu));
                break;
            case neo_hookean:
                va_start(args, type);

                this->lambda = va_arg(args, double);
                this->mu = va_arg(args, double);
                va_end(args);
                break;
        }
    }
};

#endif // !_MAT_LIB_