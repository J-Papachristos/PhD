#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define _USE_MATH_DEFINES

#ifndef _HELPER_LIB_
#define _HELPER_LIB_

/// Gauss-Quadrature
#define GQ_POINTS 2
double points_GQ[GQ_POINTS] = {-sqrt(1.0 / 3.0), +sqrt(1.0 / 3.0)}; // Gauss-Quadrature Isoparametric Points
double w_GQ[GQ_POINTS] = {+1.0, +1.0};                              // Gauss-Quadrature Weights

// #define GQ_POINTS 3
// double points_GQ[GQ_POINTS] = {-sqrt(0.6), +0, +sqrt(0.6)};        // Gauss-Quadrature Isoparametric Points
// double w_GQ[GQ_POINTS] = {+(5.0 / 9.0), +(8.0 / 9.0), +(5.0 / 9.0)}; // Gauss-Quadrature Weights

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
    quadratic
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

enum directionEnum {
    xx,
    yy,
    zz,
    xy,
    yz,
    zx,
    directions
};

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

/// @brief Helper Function to Invert a 3x3 Matrix A
void inverse3(double A[3][3], double A_inv[3][3], double *detA) {
    // Calculate Inverse
    double _A_ = +((A[1][1] * A[2][2]) - (A[1][2] * A[2][1]));
    double _B_ = -((A[1][0] * A[2][2]) - (A[1][2] * A[2][0]));
    double _C_ = +((A[1][0] * A[2][1]) - (A[1][1] * A[2][0]));
    double _D_ = -((A[0][1] * A[2][2]) - (A[0][2] * A[2][1]));
    double _E_ = +((A[0][0] * A[2][2]) - (A[0][2] * A[2][0]));
    double _F_ = -((A[0][0] * A[2][1]) - (A[0][1] * A[2][0]));
    double _G_ = +((A[0][1] * A[1][2]) - (A[0][2] * A[1][1]));
    double _H_ = -((A[0][0] * A[1][2]) - (A[0][2] * A[1][0]));
    double _I_ = +((A[0][0] * A[1][1]) - (A[0][1] * A[1][0]));

    *detA = A[0][0] * _A_ + A[0][1] * _B_ + A[0][2] * _C_;
    if (fabs(*detA) >= 1e-10) {
        A_inv[0][0] = _A_ / (*detA);
        A_inv[1][0] = _B_ / (*detA);
        A_inv[2][0] = _C_ / (*detA);
        A_inv[0][1] = _D_ / (*detA);
        A_inv[1][1] = _E_ / (*detA);
        A_inv[2][1] = _F_ / (*detA);
        A_inv[0][2] = _G_ / (*detA);
        A_inv[1][2] = _H_ / (*detA);
        A_inv[2][2] = _I_ / (*detA);
    }
}

/// @brief Kronecker Delta
int delta(int i, int j) {
    return (i == j) ? 1 : 0;
}

/// @brief Helper Function to convert from 4th
/// Order Minor-Symmetric Tensor to a 6x6 Matrix
/// (Vectorization Notation)
int fourthOrder2matrix(int i, int j) {
    if (i == j)
        return i;
    else if (i + j == xx + yy)
        return nHexDOFs + DOF_ux;
    else if (i + j == yy + zz)
        return nHexDOFs + DOF_uy;
    else if (i + j == zz + xx)
        return nHexDOFs + DOF_uz;
    return -1;
}

#endif // !_HELPER_LIB_