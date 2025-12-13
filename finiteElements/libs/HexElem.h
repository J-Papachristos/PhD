#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define _USE_MATH_DEFINES

#ifndef _HEX_LIB_
#define _HEX_LIB_

#define FREE_DOF 0    // pGroup of Free DOFs
#define NOT_PARSED -1 // pGroup for unparsed DOFs

typedef double (*funcPtr)(double, double, double);

/// Gauss-Quadrature
#define GQ_POINTS 2
double points_GQ[GQ_POINTS] = {-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)}; // Gauss-Quadrature Isoparametric Points
double w_GQ[GQ_POINTS] = {1.0, 1.0};                               // Gauss-Quadrature Weights

// #define GQ_POINTS 3
// double points_GQ[GQ_POINTS] = {-sqrt(0.6), +0, sqrt(0.6)};        // Gauss-Quadrature Isoparametric Points
// double w_GQ[GQ_POINTS] = {(5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0)}; // Gauss-Quadrature Weights

enum directionEnum {
    xx,
    yy,
    zz,
    xy,
    yz,
    zx,
    directions
};

double rho = 2.7e3; // [kg/m^3]
double E = 68e9;    // [Pa]
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

enum hex8_node {
    hex8_1,
    hex8_2,
    hex8_3,
    hex8_4,
    hex8_5,
    hex8_6,
    hex8_7,
    hex8_8,
    linear
};

enum hex20_node {
    hex20_1,
    hex20_2,
    hex20_3,
    hex20_4,
    hex20_5,
    hex20_6,
    hex20_7,
    hex20_8,
    hex20_9,
    hex20_10,
    hex20_11,
    hex20_12,
    hex20_13,
    hex20_14,
    hex20_15,
    hex20_16,
    hex20_17,
    hex20_18,
    hex20_19,
    hex20_20,
    serendipity
};

enum hex_nodeDOF {
    DOF_ux,
    DOF_uy,
    DOF_uz,
    nHexDOFs
};

/*
3D Linear Hex Element (Hex8) :
            z/ζ
             |
      {5}---------{8}
      /|     |    /|
     / |     |   / |
    /  |     +--/ -|--y/η
  {6}--+------{7}  |
   |  {1}--/---+--{4}
   |  /   /    |  /
   | /  x/ξ    | /
   |/          |/
  {2}---------{3}

3D Quadratic Hex Element (Hex20) :
          z/ζ
           |
     5-----18----8
    /|     |    /|
   17|     |   20|
  /  11    |  /  16
 6---|-19----7   |
 |   |     +-|---|--->y/η
 |   1----/10|---4
 13 /    /   15 /
 | 9    /    | 14
 |/    /     |/
 2----/12----3
     /
   x/ξ

Coefficients
|------|-Corner Nodes ( Same for Hex8 Element)-|-[ξ]--[η]--[ζ]--[η]--[ζ]--[ξ]--[ζ]--[ζ]--[ξ]--[η]--[η]--[ξ]|
|[Node]|[01]|[02]|[03]|[04]|[05]|[06]|[07]|[08]|[09]|[10]|[11]|[12]|[13]|[14]|[15]|[16]|[17]|[18]|[19]|[20]|
|  [ξ] | -1 | +1 | +1 | -1 | -1 | +1 | +1 | -1 | +0 | -1 | -1 | +1 | +1 | +0 | +1 | -1 | +0 | -1 | +1 | +0 |
|  [η] | -1 | -1 | +1 | +1 | -1 | -1 | +1 | +1 | -1 | +0 | -1 | +0 | -1 | +1 | +1 | +1 | -1 | +0 | +0 | +1 |
|  [ζ] | -1 | -1 | -1 | -1 | +1 | +1 | +1 | +1 | -1 | -1 | +0 | -1 | +0 | -1 | +0 | +0 | +1 | +1 | +1 | +1 |
*/

// double ksi_ind[serendipity] = {-1, +1, +1, -1, -1, +1, +1, -1, +0, -1, -1, +1, +1, +0, +1, -1, +0, -1, +1, +0};
// double eta_ind[serendipity] = {-1, -1, +1, +1, -1, -1, +1, +1, -1, +0, -1, +0, -1, +1, +1, +1, -1, +0, +0, +1};
// double zeta_ind[serendipity] = {-1, -1, -1, -1, +1, +1, +1, +1, -1, -1, +0, -1, +0, -1, +0, +0, +1, +1, +1, +1};

double ksi_ind[serendipity] = {-1, +1, +1, -1, -1, +1, +1, -1, +0, +1, +0, -1, +0, +1, +0, -1, -1, +1, +1, -1};
double eta_ind[serendipity] = {-1, -1, +1, +1, -1, -1, +1, +1, -1, +0, +1, +0, -1, +0, +1, +0, -1, -1, +1, +1};
double zeta_ind[serendipity] = {-1, -1, -1, -1, +1, +1, +1, +1, -1, -1, -1, -1, +1, +1, +1, +1, +0, +0, +0, +0};

class Point {
  private:
  public:
    double x;
    double y;
    double z;

    double ksi;
    double eta;
    double zeta;
};

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

    double ux_dd[elemNodes]; // x-Axis Nodal Acceleration (Array)
    double uy_dd[elemNodes]; // y-Axis Nodal Acceleration (Array)
    double uz_dd[elemNodes]; // z-Axis Nodal Acceleration (Array)

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

    // Nodal Displacement Vector [u1_x, u1_y, u1_z, ..., un_x, un_y, un_z]
    double d[nHexDOFs * elemNodes];

    // Strains
    double epsilon[nHexDOFs][nHexDOFs]; // Nodal Strain Matrix
    double epsilon_v[directions];       // Nodal Strain Vector [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx]

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
    }

    /// @brief Sets Element Values {ux_dd,uy_dd,uz_dd} to 0
    void zeroElementValues() {
        for (int i = 0; i < elemNodes; i++) {
            this->ux_dd[i] = 0;
            this->uy_dd[i] = 0;
            this->uz_dd[i] = 0;
        }
    }

    /// @brief Sets Element Values {u_x,u_y,u_z}
    /// @param u_x  Matrix (Size elemNodes) of u_x values
    /// @param u_y Matrix (Size elemNodes) of u_y values
    /// @param u_z  Matrix (Size elemNodes) of u_z values
    void setElementValues(double ux_dd[elemNodes],
                          double uy_dd[elemNodes],
                          double uz_dd[elemNodes]) {
        for (int i = 0; i < elemNodes; i++) {
            this->ux_dd[i] = ux_dd[i];
            this->uy_dd[i] = uy_dd[i];
            this->uz_dd[i] = uz_dd[i];
        }
    }

    /// @brief Helper Function, Prints 3x3 [J]
    void printJacobian() {
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[0][0], this->J[0][1], this->J[0][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[1][0], this->J[1][1], this->J[1][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J[2][0], this->J[2][1], this->J[2][2]);
    }
    /// @brief Helper Function, Prints 3x3 [J]^{-1}
    void printInvJacobian() {
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[0][0], this->J_inv[0][1], this->J_inv[0][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[1][0], this->J_inv[1][1], this->J_inv[1][2]);
        printf("|%+.8lf %+.8lf %+.8lf|\n", this->J_inv[2][0], this->J_inv[2][1], this->J_inv[2][2]);
    }
};

class hex8 : public HexElem<linear> {
  private:
    int elemNodes = linear; // Order of Element (linear,serendipity)
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
        double J_11 = +0, J_12 = +0, J_13 = 0;
        double J_21 = +0, J_22 = +0, J_23 = 0;
        double J_31 = +0, J_32 = +0, J_33 = 0;

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
        double A = +((J_22 * J_33) - (J_23 * J_32));
        double B = -((J_21 * J_33) - (J_23 * J_31));
        double C = +((J_21 * J_32) - (J_22 * J_31));
        double D = -((J_12 * J_33) - (J_13 * J_32));
        double E = +((J_11 * J_33) - (J_13 * J_31));
        double F = -((J_11 * J_32) - (J_12 * J_31));
        double G = +((J_12 * J_23) - (J_13 * J_22));
        double H = -((J_11 * J_23) - (J_13 * J_21));
        double I = +((J_11 * J_22) - (J_12 * J_21));

        this->detJ = J_11 * A + J_12 * B + J_13 * C;
        if (detJ != 0) {
            this->J_inv[0][0] = A / detJ;
            this->J_inv[1][0] = B / detJ;
            this->J_inv[2][0] = C / detJ;
            this->J_inv[0][1] = D / detJ;
            this->J_inv[1][1] = E / detJ;
            this->J_inv[2][1] = F / detJ;
            this->J_inv[0][2] = G / detJ;
            this->J_inv[1][2] = H / detJ;
            this->J_inv[2][2] = I / detJ;
        }
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
        while (err <= 1e-6) {
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

    /// @brief Calculates the Deformation Gradient [F]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getDeformGradient(double ksi, double eta, double zeta) {
        for (int i = 0; i < nHexDOFs; i++) {
            for (int j = 0; j < nHexDOFs; j++) {
                double duidxj = 0;
                for (int k = 0; k < elemNodes; k++) {
                    duidxj += (dNdksi(k, ksi, eta, zeta) * this->J_inv[j][0] +
                               dNdeta(k, ksi, eta, zeta) * this->J_inv[j][1] +
                               dNdzeta(k, ksi, eta, zeta) * this->J_inv[j][2]) *
                              this->d[k * nHexDOFs + i];
                }

                // F_ij = δ_ij + du_i/dx_j
                this->F[i][j] = ((i == j) ? 1.0 : 0.0) + duidxj;
            }
        }

        // Calculate Inverse Deformation Gradient
        double A = +((this->F[1][1] * this->F[2][2]) - (this->F[1][2] * this->F[2][1]));
        double B = -((this->F[1][0] * this->F[2][2]) - (this->F[1][2] * this->F[2][0]));
        double C = +((this->F[1][0] * this->F[2][1]) - (this->F[1][1] * this->F[2][0]));
        double D = -((this->F[0][1] * this->F[2][2]) - (this->F[0][2] * this->F[2][1]));
        double E = +((this->F[0][0] * this->F[2][2]) - (this->F[0][2] * this->F[2][0]));
        double F = -((this->F[0][0] * this->F[2][1]) - (this->F[0][1] * this->F[2][0]));
        double G = +((this->F[0][1] * this->F[1][2]) - (this->F[0][2] * this->F[1][1]));
        double H = -((this->F[0][0] * this->F[1][2]) - (this->F[0][2] * this->F[1][0]));
        double I = +((this->F[0][0] * this->F[1][1]) - (this->F[0][1] * this->F[1][0]));

        this->detF = this->F[0][0] * A + this->F[0][1] * B + this->F[0][2] * C;
        if (detF != 0) {
            this->F_inv[0][0] = A / detF;
            this->F_inv[1][0] = B / detF;
            this->F_inv[2][0] = C / detF;
            this->F_inv[0][1] = D / detF;
            this->F_inv[1][1] = E / detF;
            this->F_inv[2][1] = F / detF;
            this->F_inv[0][2] = G / detF;
            this->F_inv[1][2] = H / detF;
            this->F_inv[2][2] = I / detF;
        }
    }

    /// @brief Calculates the Linear Deformation Matrix [B_L]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getLinearDeformMatrix(double ksi, double eta, double zeta) {
        memset(this->B_L, 0, sizeof(this->B_L));
        double dNdx[this->elemNodes];
        double dNdy[this->elemNodes];
        double dNdz[this->elemNodes];

        for (int node = 0; node < this->elemNodes; node++) {
            dNdx[node] = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[0][0];
            dNdx[node] += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[0][1];
            dNdx[node] += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[0][2];

            dNdy[node] = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[1][0];
            dNdy[node] += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[1][1];
            dNdy[node] += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[1][2];

            dNdz[node] = this->dNdksi(node, ksi, eta, zeta) * this->J_inv[2][0];
            dNdz[node] += this->dNdeta(node, ksi, eta, zeta) * this->J_inv[2][1];
            dNdz[node] += this->dNdzeta(node, ksi, eta, zeta) * this->J_inv[2][2];
        }

        for (int node = 0; node < this->elemNodes; node++) {
            this->B_L[DOF_ux][node * nHexDOFs + DOF_ux] = dNdx[node];
            this->B_L[DOF_uy][node * nHexDOFs + DOF_uy] = dNdy[node];
            this->B_L[DOF_uz][node * nHexDOFs + DOF_uz] = dNdz[node];

            this->B_L[nHexDOFs + DOF_ux][node * nHexDOFs + DOF_ux] = dNdy[node];
            this->B_L[nHexDOFs + DOF_ux][node * nHexDOFs + DOF_uy] = dNdx[node];

            this->B_L[nHexDOFs + DOF_uy][node * nHexDOFs + DOF_uy] = dNdz[node];
            this->B_L[nHexDOFs + DOF_uy][node * nHexDOFs + DOF_uz] = dNdy[node];

            this->B_L[nHexDOFs + DOF_uz][node * nHexDOFs + DOF_ux] = dNdz[node];
            this->B_L[nHexDOFs + DOF_uz][node * nHexDOFs + DOF_uz] = dNdx[node];
        }
    }

    /// @brief Calculates the Non-Linear Deformation Matrix [B_NL]
    /// @param ksi  Isoparametric Coordinate ξ
    /// @param eta  Isoparametric Coordinate η
    /// @param zeta Isoparametric Coordinate ζ
    void getNonLinearDeformMatrix(double ksi, double eta, double zeta) {
        memset(this->B_NL, 0, sizeof(this->B_NL));
        double dNdx;
        double dNdy;
        double dNdz;

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
                B_NL[i * nHexDOFs + DOF_ux][node * nHexDOFs + DOF_ux] = dNdx;
                B_NL[i * nHexDOFs + DOF_uy][node * nHexDOFs + DOF_ux] = dNdy;
                B_NL[i * nHexDOFs + DOF_uz][node * nHexDOFs + DOF_ux] = dNdz;
            }
        }
    }

    void calculateStrain() {
        for (int dir = 0; dir < directions; dir++) {
            this->epsilon_v[dir] = 0;
            for (int cols = 0; cols < nHexDOFs * this->elemNodes; cols++) {
                this->epsilon_v[dir] += this->B_L[dir][cols] * this->d[cols];
            }
        }
    }

    void calculateStress() {
        this->calculateStrain();
        for (int dir = 0; dir < directions; dir++) {
            this->sigma_v[dir] = 0;
            for (int cols = 0; cols < directions; cols++) {
                this->sigma_v[dir] += D[dir][cols] * this->epsilon_v[cols];
            }
        }

        this->sigma[0][0] = this->sigma_v[0];
        this->sigma[1][1] = this->sigma_v[1];
        this->sigma[2][2] = this->sigma_v[2];

        this->sigma[0][1] = this->sigma_v[3];
        this->sigma[1][0] = this->sigma_v[3];
        this->sigma[1][2] = this->sigma_v[4];
        this->sigma[2][1] = this->sigma_v[4];
        this->sigma[2][0] = this->sigma_v[5];
        this->sigma[0][2] = this->sigma_v[5];
    }

    /// @brief Calculates the 2nd Piola Stress
    /// Matrix τ(3x3) = |F| * F^{-1} * σ * F^{-T}
    void calculatePiola2() {
        this->calculateStress();
        for (int i = 0; i < nHexDOFs; i++) {
            for (int j = i; j < nHexDOFs; j++) {
                this->S2[i][j] = 0;
                for (int l = 0; l < nHexDOFs; l++) {
                    double sum = 0;
                    for (int k = 0; k < nHexDOFs; k++) {
                        sum += this->F_inv[i][k] * this->sigma[k][l];
                    }
                    this->S2[i][j] += sum * F_inv[j][l];
                }
                this->S2[i][j] *= detF;
                (i != j) ? this->S2[j][i] = this->S2[i][j] : 0;
            }
        }

        this->S2_v[0] = this->S2[0][0];
        this->S2_v[1] = this->S2[1][1];
        this->S2_v[2] = this->S2[2][2];
        this->S2_v[3] = this->S2[0][1];
        this->S2_v[4] = this->S2[1][2];
        this->S2_v[5] = this->S2[0][2];
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
            printf("%+012.8lf, %+012.8lf [GPa]\n", this->epsilon_v[dir], this->S2_v[dir] / 1e9);
        }
        printf("--------------------------------\n");
    }
};

class hex20 : public HexElem<serendipity> {
  private:
  public:
    int elemNodes = serendipity; // Order of Element
    using HexElem<serendipity>::HexElem;

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
        for (int i = 0; i < serendipity; i++) {
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
        double A = +((J_22 * J_33) - (J_23 * J_32));
        double B_L = -((J_21 * J_33) - (J_23 * J_31));
        double C = +((J_21 * J_32) - (J_22 * J_31));
        double D = -((J_12 * J_33) - (J_13 * J_32));
        double E = +((J_11 * J_33) - (J_13 * J_31));
        double F = -((J_11 * J_32) - (J_12 * J_31));
        double G = +((J_12 * J_23) - (J_13 * J_22));
        double H = -((J_11 * J_23) - (J_13 * J_21));
        double I = +((J_11 * J_22) - (J_12 * J_21));

        this->detJ = J_11 * A + J_12 * B_L + J_13 * C;
        if (detJ != 0) {
            this->J_inv[0][0] = A / detJ;
            this->J_inv[1][0] = B_L / detJ;
            this->J_inv[2][0] = C / detJ;
            this->J_inv[0][1] = D / detJ;
            this->J_inv[1][1] = E / detJ;
            this->J_inv[2][1] = F / detJ;
            this->J_inv[0][2] = G / detJ;
            this->J_inv[1][2] = H / detJ;
            this->J_inv[2][2] = I / detJ;
        }
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