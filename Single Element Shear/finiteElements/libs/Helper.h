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

#endif // !_HELPER_LIB_