#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "Helper.h"
#include "Material.h"

#define _USE_MATH_DEFINES

#ifndef _HEX_LIB_
#define _HEX_LIB_

#define FREE_DOF 0    // pGroup of Free DOFs
#define NOT_PARSED -1 // pGroup for unparsed DOFs

typedef double (*funcPtr)(double, double, double);

double ksi_ind[quadratic] = {-1, +1, +1, -1, -1, +1, +1, -1, +0, +1, +0, -1, +0, +1, +0, -1, -1, +1, +1, -1};
double eta_ind[quadratic] = {-1, -1, +1, +1, -1, -1, +1, +1, -1, +0, +1, +0, -1, +0, +1, +0, -1, -1, +1, +1};
double zeta_ind[quadratic] = {-1, -1, -1, -1, +1, +1, +1, +1, -1, -1, -1, -1, +1, +1, +1, +1, +0, +0, +0, +0};

template <int elemNodes>
class HexElem {
  private:
    uint64_t nelem; // ID
  public:
    uint64_t nDOF[elemNodes]; // First DOF of each node
    bool isBoundary;          // True, if element is on a boundary

    // Nodal Values
    double x[elemNodes]; // Real Coordinate x per Node (Array)
    double y[elemNodes]; // Real Coordinate y per Node (Array)
    double z[elemNodes]; // Real Coordinate z per Node (Array)

    int nodeIndex[elemNodes];  // Node Indices (Array)
    int pGroupElem[elemNodes]; // Node Group Indices (Array)

    /// Jacobian
    double J[nHexDOFs][nHexDOFs];     // Jacobian (3x3)
    double J_inv[nHexDOFs][nHexDOFs]; // Inverse Jacobian (3x3)
    double detJ = 0;                  // Jacobian Determinant

    /// Deformation Matrix
    double B_L[directions][nHexDOFs * elemNodes];
    double B_NL[nHexDOFs * nHexDOFs][nHexDOFs * elemNodes];

    /// Deformation Gradient
    double F[nHexDOFs][nHexDOFs];     // Deformation Gradient (3x3)
    double F_inv[nHexDOFs][nHexDOFs]; // Inverse Deformation Gradient (3x3)
    double detF = 0;                  // Deformation Gradient Determinant

    // Cauchy-Green Strain Tensor
    double C[nHexDOFs][nHexDOFs];     // Right Cauchy-Green Strain Tensor (3x3)
    double C_inv[nHexDOFs][nHexDOFs]; // Inverse Right Cauchy-Green Strain Tensor (3x3)
    double detC = 0;                  // Right Cauchy-Green Strain Tensor Determinant

    // Elasticity Matrix
    double D[directions][directions];

    // Nodal Displacement Vector [u1_x, u1_y, u1_z, ..., un_x, un_y, un_z]
    double d[nHexDOFs * elemNodes];

    // Strains
    double epsilon[nHexDOFs][nHexDOFs]; // Nodal Strain Matrix
    double epsilon_v[directions];       // Nodal Strain Vector [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]
    double eps_GL[nHexDOFs][nHexDOFs];  // Green-Lagrange Strain Tensor
    double eps_GL_v[directions];        // Green-Lagrange Strain Tensor [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]

    // Stresses
    double sigma[nHexDOFs][nHexDOFs]; // Nodal Stress Matrix
    double sigma_v[directions];       // Nodal Stress Vector
    double S2[nHexDOFs][nHexDOFs];    // 2nd Piola Stress Matrix
    double S2_v[directions];          // 2nd Piola Stress Vector

    /// Methods

    /// @brief Default Constructor |
    /// Sets nelem to 0
    HexElem() {
        this->nelem = 0;
    }

    /// @brief Default Constructor
    /// @param nelem ID of Element
    /// @param x[elemNodes] Element's Nodal x Coordinates
    /// @param y[elemNodes] Element's Nodal y Coordinates
    /// @param z[elemNodes] Element's Nodal z Coordinates
    /// @param nodeIndex[elemNodes] Array of Node Indices
    /// @param pGroupElem[elemNodes] Array of Node's Group Index
    HexElem(int nelem,
            double x[elemNodes],
            double y[elemNodes],
            double z[elemNodes],
            int nodeIndex[elemNodes],
            int pGroupElem[elemNodes]) {
        this->nelem = nelem;
        for (int i = 0; i < elemNodes; i++) {
            this->x[i] = x[i];
            this->y[i] = y[i];
            this->z[i] = z[i];

            this->nodeIndex[i] = nodeIndex[i];
            this->pGroupElem[i] = pGroupElem[i];
        }

        for (int i = 0; i < nHexDOFs * elemNodes; i++) {
            this->d[i] = 0.0;
        }
    }

    int getNelem() {
        return this->nelem;
    }

    /// @brief Helper Function, Prints 3x3 [J]
    void printJacobian() {
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[0][0], this->J[0][1], this->J[0][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[1][0], this->J[1][1], this->J[1][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[2][0], this->J[2][1], this->J[2][2]);
        printf("\n");
    }
    /// @brief Helper Function, Prints 3x3 [J]^{-1}
    void printInvJacobian() {
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[0][0], this->J_inv[0][1], this->J_inv[0][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[1][0], this->J_inv[1][1], this->J_inv[1][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[2][0], this->J_inv[2][1], this->J_inv[2][2]);
        printf("\n");
    }
};

class hex8 : public HexElem<linear> {
  private:
    int elemNodes = linear; // Order of Element (linear,quadratic)
  public:
    using HexElem<linear>::HexElem;

    int getElemNodes() {
        return this->elemNodes;
    }

    /// @brief Shape Functions for hex8 Element
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return N_i - hex8 Shape Function i at {ξ,η,ζ}
    double N(int i, double ksi, double eta, double zeta) {
        return 0.125 * (1 + ksi * ksi_ind[i]) * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]);
    }

    /// @brief Shape Function Derivatives for hex8 Element with respect to ξ at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dξ
    double dNdksi(int i, double ksi, double eta, double zeta) {
        return 0.125 * ksi_ind[i] * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]);
    }

    /// @brief Shape Function Derivatives for hex8 Element with respect to η at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dη
    double dNdeta(int i, double ksi, double eta, double zeta) {
        return 0.125 * eta_ind[i] * (1 + ksi * ksi_ind[i]) * (1 + zeta * zeta_ind[i]);
    }

    /// @brief Shape Function Derivatives for hex8 Element with respect to ζ at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dζ
    double dNdzeta(int i, double ksi, double eta, double zeta) {
        return 0.125 * zeta_ind[i] * (1 + ksi * ksi_ind[i]) * (1 + eta * eta_ind[i]);
    }

    /// @brief Calculates the Jacobian [J] and
    /// the Inverse Jacobian [J]^{-1} for the
    /// hex8 Element
    /// @param ksi  ξ Isoparametric Variable
    /// @param eta  η Isoparametric Variable
    /// @param zeta ζ Isoparametric Variable
    void jacobian(double ksi, double eta, double zeta) {
        // Init
        double J_11 = 0, J_12 = 0, J_13 = 0;
        double J_21 = 0, J_22 = 0, J_23 = 0;
        double J_31 = 0, J_32 = 0, J_33 = 0;

        // Calculate Jacobian
        for (int i = hex8_1; i < linear; i++) {
            J_11 += dNdksi(i, ksi, eta, zeta) * x[i];
            J_21 += dNdeta(i, ksi, eta, zeta) * x[i];
            J_31 += dNdzeta(i, ksi, eta, zeta) * x[i];

            J_12 += dNdksi(i, ksi, eta, zeta) * y[i];
            J_22 += dNdeta(i, ksi, eta, zeta) * y[i];
            J_32 += dNdzeta(i, ksi, eta, zeta) * y[i];

            J_13 += dNdksi(i, ksi, eta, zeta) * z[i];
            J_23 += dNdeta(i, ksi, eta, zeta) * z[i];
            J_33 += dNdzeta(i, ksi, eta, zeta) * z[i];
        }
        this->J[0][0] = J_11, this->J[0][1] = J_12, this->J[0][2] = J_13;
        this->J[1][0] = J_21, this->J[1][1] = J_22, this->J[1][2] = J_23;
        this->J[2][0] = J_31, this->J[2][1] = J_32, this->J[2][2] = J_33;

        // Calculate Inverse Jacobian
        inverse3(this->J, this->J_inv, &this->detJ);
    }

    double getVolume() {
        if (this->detJ == 0) {
            this->jacobian(0, 0, 0);
        }
        return (2.0 * 2.0 * 2.0) * this->detJ; // Fix!
    }

    /// @brief Transform from Cartesian {x,y,z}
    /// Coordinates to Isoparametric {ξ,η,ζ} Coordinates
    /// @param p Point p, with members x,y,z defined
    void cart2iso(Point *p) {
        p->ksi = 0, p->eta = 0, p->zeta = 0;

        // Newton-Raphson Method
        double err = INFINITY;
        while (err >= 1e-6) {
            double x_t = -p->x, y_t = -p->y, z_t = -p->z; // Target Values / Functions
            for (int i = 0; i < this->elemNodes; i++) {
                x_t += this->N(i, p->ksi, p->eta, p->zeta) * this->x[i];
                y_t += this->N(i, p->ksi, p->eta, p->zeta) * this->y[i];
                z_t += this->N(i, p->ksi, p->eta, p->zeta) * this->z[i];
                this->jacobian(p->ksi, p->eta, p->zeta);
            }
            p->ksi -= this->J_inv[0][0] * x_t;
            p->ksi -= this->J_inv[0][1] * y_t;
            p->ksi -= this->J_inv[0][2] * z_t;

            p->eta -= this->J_inv[1][0] * x_t;
            p->eta -= this->J_inv[1][1] * y_t;
            p->eta -= this->J_inv[1][2] * z_t;

            p->zeta -= this->J_inv[2][0] * x_t;
            p->zeta -= this->J_inv[2][1] * y_t;
            p->zeta -= this->J_inv[2][2] * z_t;

            err = sqrt(x_t * x_t + y_t * y_t + z_t * z_t);
        }
    }

    /// @brief Transform from Isoparametric {ξ,η,ζ}
    /// Coordinates to Cartesian {x,y,z} Coordinates
    /// @param p Point p, with members ksi,eta,zeta defined
    void iso2cart(Point *p) {
        p->x = 0, p->y = 0, p->z = 0;
        for (int i = 0; i < this->elemNodes; i++) {
            p->x += this->N(i, p->ksi, p->eta, p->zeta) * this->x[i];
            p->y += this->N(i, p->ksi, p->eta, p->zeta) * this->y[i];
            p->z += this->N(i, p->ksi, p->eta, p->zeta) * this->z[i];
        }
    }

    void calcD(Material mat) {
        for (int i = 0; i < directions; i++) {
            for (int j = 0; j < directions; j++) {
                (i < directions / 2 && j < directions / 2) ? this->D[i][j] = mat.lambda : this->D[i][j] = 0;
            }
            (i < directions / 2) ? this->D[i][i] = mat.lambda + 2 * mat.mu : this->D[i][i] = mat.mu;
        }
    }

    /// @brief Calculates the Deformation Gradient [F]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getDeformGradient(double ksi, double eta, double zeta) {
        this->jacobian(ksi, eta, zeta);
        for (int i = 0; i < nHexDOFs; i++) {
            for (int j = 0; j < nHexDOFs; j++) {
                double duidxj = 0;
                for (int node = 0; node < elemNodes; node++) {
                    // duidxj += (dNdksi(node, ksi, eta, zeta) * this->J_inv[j][0] +
                    //            dNdeta(node, ksi, eta, zeta) * this->J_inv[j][1] +
                    //            dNdzeta(node, ksi, eta, zeta) * this->J_inv[j][2]) *
                    //           this->d[node * nHexDOFs + i];
                    switch (j) {
                        case DOF_ux:
                            duidxj += dNdksi(node, ksi, eta, zeta) * this->d[node * nHexDOFs + i];
                            break;
                        case DOF_uy:
                            duidxj += dNdeta(node, ksi, eta, zeta) * this->d[node * nHexDOFs + i];
                            break;
                        case DOF_uz:
                            duidxj += dNdzeta(node, ksi, eta, zeta) * this->d[node * nHexDOFs + i];
                            break;
                    }
                }

                // F_ij = δ_ij + du_i/dx_j
                this->F[i][j] = ((i == j) ? 1.0 : 0.0) + duidxj;
            }
        }

        // Calculate Inverse Deformation Gradient
        inverse3(this->F, this->F_inv, &this->detF);

        if (detF <= 1e-10) {
            printf("det{F} < 0\n");
            this->printDeformGradient();
            exit(1);
        } else if (isnan(detF)) {
            printf("det{F} is NaN!\n");
            exit(1);
        }
    }

    /// @brief Calculates the Right Cauchy-Green Strain Tensor [C] = [F]^T * [F]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getRightCauchy(double ksi, double eta, double zeta) {
        this->getDeformGradient(ksi, eta, zeta);
        for (int i = 0; i < nHexDOFs; i++) {
            for (int j = 0; j < nHexDOFs; j++) {
                this->C[i][j] = 0.0;
                for (int k = 0; k < nHexDOFs; k++) {
                    this->C[i][j] += F[k][i] * F[k][j];
                }
            }
        }

        // Calculate Inverse Cauchy
        inverse3(this->C, this->C_inv, &this->detC);
    }

    /// @brief Calculates the Linear Deformation Matrix [B_L]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getLinearDeformMatrix(double ksi, double eta, double zeta) {
        this->getDeformGradient(ksi, eta, zeta);
        // memset(this->B_L, 0, sizeof(this->B_L));
        for (int i = 0; i < directions; i++) {
            for (int j = 0; j < nHexDOFs * elemNodes; j++) {
                this->B_L[i][j] = 0;
            }
        }
        double dNdx = 0;
        double dNdy = 0;
        double dNdz = 0;

        for (int node = 0; node < this->elemNodes; node++) {
            dNdx = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[0][0];
            dNdx += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[0][1];
            dNdx += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[0][2];

            dNdy = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[1][0];
            dNdy += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[1][1];
            dNdy += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[1][2];

            dNdz = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[2][0];
            dNdz += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[2][1];
            dNdz += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[2][2];

            this->B_L[DOF_ux][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_ux] * dNdx;
            this->B_L[DOF_uy][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_uy] * dNdy;
            this->B_L[DOF_uz][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_uz] * dNdz;

            this->B_L[DOF_ux][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_ux] * dNdx;
            this->B_L[DOF_uy][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_uy] * dNdy;
            this->B_L[DOF_uz][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_uz] * dNdz;

            this->B_L[DOF_ux][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_ux] * dNdx;
            this->B_L[DOF_uy][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_uy] * dNdy;
            this->B_L[DOF_uz][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_uz] * dNdz;

            this->B_L[DOF_ux + nHexDOFs][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_ux] * dNdy + this->F[DOF_ux][DOF_uy] * dNdx;
            this->B_L[DOF_uy + nHexDOFs][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_uy] * dNdz + this->F[DOF_ux][DOF_uz] * dNdy;
            this->B_L[DOF_uz + nHexDOFs][node * nHexDOFs + DOF_ux] = this->F[DOF_ux][DOF_uz] * dNdx + this->F[DOF_ux][DOF_ux] * dNdz;

            this->B_L[DOF_ux + nHexDOFs][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_ux] * dNdy + this->F[DOF_uy][DOF_uy] * dNdx;
            this->B_L[DOF_uy + nHexDOFs][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_uy] * dNdz + this->F[DOF_uy][DOF_uz] * dNdy;
            this->B_L[DOF_uz + nHexDOFs][node * nHexDOFs + DOF_uy] = this->F[DOF_uy][DOF_uz] * dNdx + this->F[DOF_uy][DOF_ux] * dNdz;

            this->B_L[DOF_ux + nHexDOFs][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_ux] * dNdy + this->F[DOF_uz][DOF_uy] * dNdx;
            this->B_L[DOF_uy + nHexDOFs][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_uy] * dNdz + this->F[DOF_uz][DOF_uz] * dNdy;
            this->B_L[DOF_uz + nHexDOFs][node * nHexDOFs + DOF_uz] = this->F[DOF_uz][DOF_uz] * dNdx + this->F[DOF_uz][DOF_ux] * dNdz;
        }
    }

    /// @brief Calculates the Non-Linear Deformation Matrix [B_NL]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getNonLinearDeformMatrix(double ksi, double eta, double zeta) {
        this->jacobian(ksi, eta, zeta);
        // memset(this->B_NL, 0, sizeof(this->B_NL));
        for (int i = 0; i < nHexDOFs * nHexDOFs; i++) {
            for (int j = 0; j < nHexDOFs * elemNodes; j++) {
                this->B_NL[i][j] = 0;
            }
        }

        double dNdx = 0;
        double dNdy = 0;
        double dNdz = 0;

        for (int node = 0; node < this->elemNodes; node++) {
            dNdx = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[0][0];
            dNdx += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[0][1];
            dNdx += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[0][2];

            dNdy = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[1][0];
            dNdy += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[1][1];
            dNdy += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[1][2];

            dNdz = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[2][0];
            dNdz += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[2][1];
            dNdz += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[2][2];

            for (int i = 0; i < nHexDOFs; i++) {
                this->B_NL[i * nHexDOFs + DOF_ux][node * nHexDOFs + i] = dNdx;
                this->B_NL[i * nHexDOFs + DOF_uy][node * nHexDOFs + i] = dNdy;
                this->B_NL[i * nHexDOFs + DOF_uz][node * nHexDOFs + i] = dNdz;
            }
        }
    }

    void calculateStrain() {
        for (int dir = 0; dir < directions; dir++) {
            this->epsilon_v[dir] = 0;
            for (int cols = 0; cols < nHexDOFs * this->elemNodes; cols++) {
                this->epsilon_v[dir] += (dir > nHexDOFs ? 2.0 : 1.0) * this->B_L[dir][cols] * this->d[cols];
            }
        }
    }

    void calculateGreenLagrangeStrain(double ksi, double eta, double zeta) {
        this->getRightCauchy(ksi, eta, zeta);
        for (int i = 0; i < nHexDOFs; i++) {
            for (int j = 0; j < nHexDOFs; j++) {
                this->eps_GL[i][j] = 0.5 * (this->C[i][j] - ((i == j) ? 1. : 0.));
            }
        }

        this->eps_GL_v[DOF_ux] = this->eps_GL[0][0];
        this->eps_GL_v[DOF_uy] = this->eps_GL[1][1];
        this->eps_GL_v[DOF_uz] = this->eps_GL[2][2];
        this->eps_GL_v[nHexDOFs + DOF_ux] = 2.0 * this->eps_GL[0][1];
        this->eps_GL_v[nHexDOFs + DOF_uy] = 2.0 * this->eps_GL[1][2];
        this->eps_GL_v[nHexDOFs + DOF_uz] = 2.0 * this->eps_GL[0][2];
    }

    void calculateStress(Material mat) {
        this->calculateStrain();
        for (int dir = 0; dir < directions; dir++) {
            this->sigma_v[dir] = 0;
            for (int cols = 0; cols < directions; cols++) {
                this->sigma_v[dir] += this->D[dir][cols] * this->epsilon_v[cols];
            }
        }

        this->sigma[0][0] = this->sigma_v[0]; // σ_xx
        this->sigma[1][1] = this->sigma_v[1]; // σ_yy
        this->sigma[2][2] = this->sigma_v[2]; // σ_zz

        this->sigma[0][1] = this->sigma_v[3]; // σ_xy
        this->sigma[1][0] = this->sigma_v[3];
        this->sigma[1][2] = this->sigma_v[4]; // σ_yz
        this->sigma[2][1] = this->sigma_v[4];
        this->sigma[2][0] = this->sigma_v[5]; // σ_zx
        this->sigma[0][2] = this->sigma_v[5];
    }

    /// @brief Calculates the 2nd Piola Stress
    /// Matrix τ(3x3) = |F| * F^{-1} * σ * F^{-T}
    void calculatePiola2(double ksi, double eta, double zeta, Material mat) {
        this->calculateGreenLagrangeStrain(ksi, eta, zeta);

        for (int i = 0; i < directions; i++) {
            this->S2_v[i] = 0;
            for (int j = 0; j < directions; j++) {
                this->S2_v[i] += this->D[i][j] * this->eps_GL_v[j];
            }
        }

        this->S2[DOF_ux][DOF_ux] = this->S2_v[DOF_ux];
        this->S2[DOF_uy][DOF_uy] = this->S2_v[DOF_uy];
        this->S2[DOF_uz][DOF_uz] = this->S2_v[DOF_uz];

        this->S2[DOF_ux][DOF_uy] = this->S2_v[nHexDOFs + DOF_ux];
        this->S2[DOF_uy][DOF_ux] = this->S2_v[nHexDOFs + DOF_ux];

        this->S2[DOF_uy][DOF_uz] = this->S2_v[nHexDOFs + DOF_uy];
        this->S2[DOF_uz][DOF_uy] = this->S2_v[nHexDOFs + DOF_uy];

        this->S2[DOF_uz][DOF_ux] = this->S2_v[nHexDOFs + DOF_uz];
        this->S2[DOF_ux][DOF_uz] = this->S2_v[nHexDOFs + DOF_uz];
    }

    /// @brief Debug Method, Prints Cauchy Stress and Engineering Strain Matrices
    void printStressStrain() {
        for (int dir = 0; dir < directions; dir++) {
            printf("%+012.8lf, %+012.8lf [GPa]\n", this->epsilon_v[dir], this->sigma_v[dir] / 1e9);
        }
        printf("--------------------------------\n");
    }

    /// @brief Debug Method, Prints 2nd Piola Stress and Engineering Strain Matrices
    void printStressStrain2() {
        for (int dir = 0; dir < directions; dir++) {
            printf("%+012.8lf, %+012.8lf [GPa]\n", this->eps_GL_v[dir], this->S2_v[dir] / 1e9);
        }
        printf("--------------------------------\n");
    }

    void printDeformGradient() {
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->F[0][0], this->F[0][1], this->F[0][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->F[1][0], this->F[1][1], this->F[1][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->F[2][0], this->F[2][1], this->F[2][2]);
        printf("\n");
    }

    void printDisplacements() {
        printf("Element %d:\n", this->getNelem());
        for (int j = 0; j < elemNodes; j++) {
            printf("\tNode[%d] = {%lf, %lf, %lf} [mm]\n", j,
                   this->d[j * nHexDOFs + DOF_ux],
                   this->d[j * nHexDOFs + DOF_uy],
                   this->d[j * nHexDOFs + DOF_uz]);
        }
    }
};

class hex20 : public HexElem<quadratic> {
  private:
  public:
    int elemNodes = quadratic; // Order of Element
    using HexElem<quadratic>::HexElem;

    /// @brief Shape Functions for hex20 Element
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return N_i - hex20 Shape Function i at {ξ,η,ζ}
    double N(int i, double ksi, double eta, double zeta) {
        switch (i + 1) {
            // Edge Nodes
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                return 0.125 * (1 + ksi * ksi_ind[i]) * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]) *
                       (ksi * ksi_ind[i] + eta * eta_ind[i] + zeta * zeta_ind[i] - 2);

            // ξ_i = 0
            case 9:
            case 11:
            case 13:
            case 15:
                return 0.250 * (1 - ksi * ksi) * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]);

            // η_i = 0
            case 10:
            case 12:
            case 14:
            case 16:
                return 0.250 * (1 - eta * eta) * (1 + ksi * ksi_ind[i]) * (1 + zeta * zeta_ind[i]);

            // ζ_i = 0
            case 17:
            case 18:
            case 19:
            case 20:
                return 0.250 * (1 - zeta * zeta) * (1 + ksi * ksi_ind[i]) * (1 + eta * eta_ind[i]);
            default:
                return 0;
        }
    }

    //// Shape Function Derivatives

    /// @brief Shape Function Derivatives for hex20 Element with respect to ξ at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dξ
    double dNdksi(int i, double ksi, double eta, double zeta) {
        switch (i + 1) {
            // Edge Nodes
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                return 0.125 * ksi_ind[i] * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]) * (2 * ksi * ksi_ind[i] + eta * eta_ind[i] + zeta * zeta_ind[i] - 1);

            // ξ_i = 0
            case 9:
            case 11:
            case 13:
            case 15:
                return -0.50 * ksi * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]);

            // η_i = 0
            case 10:
            case 12:
            case 14:
            case 16:
                return 0.250 * ksi_ind[i] * (1 - eta * eta) * (1 + zeta * zeta_ind[i]);

            // ζ_i = 0
            case 17:
            case 18:
            case 19:
            case 20:
                return 0.250 * ksi_ind[i] * (1 - zeta * zeta) * (1 + eta * eta_ind[i]);
            default:
                return 0;
        }
    }

    /// @brief Shape Function Derivatives for hex20 Element with respect to η at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dη
    double dNdeta(int i, double ksi, double eta, double zeta) {
        switch (i + 1) {
            // Edge Nodes
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                return 0.125 * eta_ind[i] * (1 + ksi * ksi_ind[i]) * (1 + zeta * zeta_ind[i]) * (ksi * ksi_ind[i] + 2 * eta * eta_ind[i] + zeta * zeta_ind[i] - 1);

            // ξ_i = 0
            case 9:
            case 11:
            case 13:
            case 15:
                return 0.250 * eta_ind[i] * (1 - ksi * ksi) * (1 + zeta * zeta_ind[i]);

            // η_i = 0
            case 10:
            case 12:
            case 14:
            case 16:
                return -0.50 * eta * (1 + ksi * ksi_ind[i]) * (1 + zeta * zeta_ind[i]);

            // ζ_i = 0
            case 17:
            case 18:
            case 19:
            case 20:
                return 0.250 * eta_ind[i] * (1 - zeta * zeta) * (1 + ksi * ksi_ind[i]);
            default:
                return 0;
        }
    }

    /// @brief Shape Function Derivatives for hex20 Element with respect to ζ at {ξ,η,ζ}
    /// @param i Selection of Shape Function
    /// @param ksi  ξ Coordinate [-1,1]
    /// @param eta  η Coordinate [-1,1]
    /// @param zeta ζ Coordinate [-1,1]
    /// @return dN_i/dζ
    double dNdzeta(int i, double ksi, double eta, double zeta) {
        switch (i + 1) {
            // Edge Nodes
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
                return 0.125 * zeta_ind[i] * (1 + ksi * ksi_ind[i]) * (1 + eta * eta_ind[i]) * (ksi * ksi_ind[i] + eta * eta_ind[i] + 2 * zeta * zeta_ind[i] - 1);

            // ξ_i = 0
            case 9:
            case 11:
            case 13:
            case 15:
                return 0.250 * zeta_ind[i] * (1 - ksi * ksi) * (1 + eta * eta_ind[i]);

            // η_i = 0
            case 10:
            case 12:
            case 14:
            case 16:
                return 0.250 * zeta_ind[i] * (1 - eta * eta) * (1 + ksi * ksi_ind[i]);

            // ζ_i = 0
            case 17:
            case 18:
            case 19:
            case 20:
                return -0.50 * zeta * (1 + eta * eta_ind[i]) * (1 + ksi * ksi_ind[i]);
            default:
                return 0;
        }
    }

    /*
        //// Shape Function Second Derivatives
        // d^2N/dξ^2
        double d2Ndksi2(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return 2 * ksi_ind[i] * ksi_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                    return -0.50 * (1 + eta * eta_ind[i]) * (1 + zeta * zeta_ind[i]);

                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return 0;
                default:
                    return 0;
            }
        }

        // d^2N/dη^2
        double d2Ndeta2(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return 2 * eta_ind[i] * eta_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                    return 0;

                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                    return -0.50 * (1 + ksi * ksi_ind[i]) * (1 + zeta * zeta_ind[i]);

                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return 0;
                default:
                    return 0;
            }
        }

        // d^2N/dζ^2
        double d2Ndzeta2(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return 2 * zeta_ind[i] * zeta_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                    return 0;

                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return -0.50 * (1 + eta * eta_ind[i]) * (1 + ksi * ksi_ind[i]);
                default:
                    return 0;
            }
        }

        // d^2N/dξη=d^2N/dηξ
        double d2Ndksi_eta(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return eta_ind[i] * ksi_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                    return -0.50 * ksi * eta_ind[i] * (1 + zeta * zeta_ind[i]);

                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                    return -0.50 * ksi_ind[i] * eta * (1 + zeta * zeta_ind[i]);

                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return 0.250 * ksi_ind[i] * eta_ind[i] * (1 - zeta * zeta);
                default:
                    return 0;
            }
        }

        // d^2N/dξζ=d^2N/dζξ
        double d2Ndksi_zeta(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return zeta_ind[i] * ksi_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                    return -0.50 * ksi * zeta_ind[i] * (1 + eta * eta_ind[i]);

                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                    return 0.250 * ksi_ind[i] * zeta_ind[i] * (1 - eta * eta);

                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return -0.50 * ksi_ind[i] * zeta * (1 + eta * eta_ind[i]);
                default:
                    return 0;
            }
        }

        // d^2N/dηζ=d^2N/dζη
        double d2Ndeta_zeta(int i, double ksi, double eta, double zeta) {
            switch (i + 1) {
                // Edge Nodes
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                    return zeta_ind[i] * eta_ind[i];

                // ξ_i = 0
                case 9:
                case 11:
                case 13:
                case 15:
                    return 0.250 * eta_ind[i] * zeta_ind[i] * (1 - ksi * ksi);

                // η_i = 0
                case 10:
                case 12:
                case 14:
                case 16:
                    return -0.50 * eta * (1 + ksi * ksi_ind[i]) * zeta_ind[i];

                // ζ_i = 0
                case 17:
                case 18:
                case 19:
                case 20:
                    return -0.50 * zeta * (1 + ksi * ksi_ind[i]) * (1 + eta_ind[i]);
                default:
                    return 0;
            }
        }
        */

    /// @brief Calculates the Jacobian [J] and
    /// the Inverse Jacobian [J]^{-1} for the
    /// Hexahedral Element
    /// @param ksi ξ Isoparametric Variable
    /// @param eta η Isoparametric Variable
    /// @param zeta ζ Isoparametric Variable
    void jacobian(double ksi, double eta, double zeta) {
        // Init
        double J_11 = +0, J_12 = +0, J_13 = 0;
        double J_21 = +0, J_22 = +0, J_23 = 0;
        double J_31 = +0, J_32 = +0, J_33 = 0;

        // Calculate Jacobian
        for (int i = 0; i < quadratic; i++) {
            J_11 += dNdksi(i, ksi, eta, zeta) * x[i];
            J_21 += dNdeta(i, ksi, eta, zeta) * x[i];
            J_31 += dNdzeta(i, ksi, eta, zeta) * x[i];

            J_12 += dNdksi(i, ksi, eta, zeta) * y[i];
            J_22 += dNdeta(i, ksi, eta, zeta) * y[i];
            J_32 += dNdzeta(i, ksi, eta, zeta) * y[i];

            J_13 += dNdksi(i, ksi, eta, zeta) * z[i];
            J_23 += dNdeta(i, ksi, eta, zeta) * z[i];
            J_33 += dNdzeta(i, ksi, eta, zeta) * z[i];
        }
        this->J[0][0] = J_11, this->J[0][1] = J_12, this->J[0][2] = J_13;
        this->J[1][0] = J_21, this->J[1][1] = J_22, this->J[1][2] = J_23;
        this->J[2][0] = J_31, this->J[2][1] = J_32, this->J[2][2] = J_33;

        // Calculate Inverse Jacobian
        inverse3(this->J, this->J_inv, &this->detJ);
    }

    double getVolume() {
        if (this->detJ == 0) {
            this->jacobian(0, 0, 0);
        }
        return (2.0 * 2.0 * 2.0) * this->detJ;
    }
};

/// @brief Body of Elements (type hexType)
/// @tparam hexType Type of Element (hex8/hex20)
template <typename hexType>
class Body {
  private:
    bool isMasterBody;

  public:
    hexType *elemArray; // Array of Body Elements
    int *elemNumber;    // Global Numbering of each Element

    /// Constructors
    Body(hexType *elemArray, int nElems) {
        elemArray = (hexType *) malloc(nElems * sizeof(hexType));
        elemNumber = (int *) malloc(nElems * sizeof(int));
        for (int i = 0; i < nElems; i++) {
            this->elemArray[i] = elemArray[nElems];
        }
    }

    /// Destructor
    ~Body() {
        free(elemArray);
        free(elemNumber);
    }
    /// Methods
    inline bool isMaster() { return this->isMasterBody; }
    inline void setMaster() { this->isMasterBody = true; }
    inline void setSlave() { this->isMasterBody = false; }
};

#endif // !_HEX_LIB_