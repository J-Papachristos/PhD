#ifndef _NUMERICAL_LIB_
#define _NUMERICAL_LIB_

#define _USE_MATH_DEFINES
#include "Sparse.h"
#include <math.h>

/// @brief Sparse Matrix Implementation of the Bi-Conjugate
/// Gradient Stabilized Algorithm for Linear Systems
/// @param A Sparse Matrix [A] (pointer to Sparse Row)
/// @param b RHS Vector {b}
/// @param x Solution Vector {x}
/// @param n Size of [A]
/// @return # Iterations
int BiCGSTAB(Sparse *A, double *b, double *x, int n) {
    double rs_err = 1e-12; // Tolerance

    // Initialize Vectors
    double *r0 = (double *) malloc(sizeof(double) * n);
    double *r0_hat = (double *) malloc(sizeof(double) * n);
    double *p0 = (double *) malloc(sizeof(double) * n);
    double *s = (double *) malloc(sizeof(double) * n);
    double *h = (double *) malloc(sizeof(double) * n);
    double rs_old = 0, rs_new = 0, rho0 = 0;

    double *Ap = (double *) malloc(sizeof(double) * n);
    double *t = (double *) malloc(sizeof(double) * n);

    // {r0} = {b} - [A]*{x}
    memcpy(r0, b, sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < A[i].n_nz - 1; j++) {
            r0[i] -= A[i].v[j] * x[A[i].c[j]];
        }
        rs_old += (r0[i] * r0[i]);

        // ρ0 = {r0_hat} * {r0}, Dot Product
        rho0 += r0[i] * r0[i];
    }

    // {r0_hat} = {r0}
    memcpy(r0_hat, r0, sizeof(double) * n);

    // {p0} = {r0}
    memcpy(p0, r0, sizeof(double) * n);

    int noIters = 0;
    while (++noIters < n) {
        // Save prev loop
        double rho0_old = rho0;

        double r0_Ap = 0;
        // {Ap} or {v} = [A] * {p_i-1}
        for (int i = 0; i < n; i++) {
            Ap[i] = 0;
            for (int j = 0; j < A[i].n_nz - 1; j++) {
                Ap[i] += A[i].v[j] * p0[A[i].c[j]];
            }
            // α = ρ0 / ({r0_hat} * {v}) = ρ0 / ({r0} * {Ap})
            r0_Ap += r0_hat[i] * Ap[i];
        }

        if (r0_Ap == 0) {
            return noIters;
        }
        double alpha = rho0_old / r0_Ap;

        // {x} = {x_i-1} + α * p_i-1
        // {r} = {r_i-1} - α * Ap
        for (int i = 0; i < n; i++) {
            h[i] = x[i] + alpha * p0[i];
            s[i] = r0[i] - alpha * Ap[i];
        }

        // {t} = [A] * {r_new}
        // ω = {t} * {s} / ({t} * {t})
        double ts = 0, tt = 0;
        for (int i = 0; i < n; i++) {
            t[i] = 0;
            for (int j = 0; j < A[i].n_nz - 1; j++) {
                t[i] += A[i].v[j] * s[A[i].c[j]];
            }
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }
        double omega = ts / tt;

        // {x_i} = {x_i-1} + ω*{s}
        // {r_i} = {s} - ω*{t}
        rs_new = 0;
        for (int i = 0; i < n; i++) {
            x[i] = h[i] + omega * s[i];
            r0[i] = s[i] - omega * t[i];
            rs_new += r0[i] * r0[i];
        }

        if (rs_old != 0) {
            // printf("res_norm = %lf\n", fabs((rs_new - rs_old) / rs_old + 1));
            if ((fabs((rs_new - rs_old) / rs_old + 1) <= rs_err) || (sqrt(rs_new) <= rs_err)) {
                // printf("\nIterations of BiCGSTAB: %d\n", noIters);
                break;
            }
        } else {
            printf("Residual = 0, Iterations of BiCGSTAB: %d\n", noIters);
            return noIters;
        }

        rho0 = 0;
        for (int i = 0; i < n; i++) {
            rho0 += r0_hat[i] * r0[i];
        }
        double beta = (rho0 / rho0_old) * (alpha / omega);

        for (int i = 0; i < n; i++) {
            p0[i] = r0[i] + beta * (p0[i] - omega * Ap[i]);
        }
    }

    free(r0);
    free(r0_hat);
    free(p0);
    free(s);
    free(h);

    free(Ap);
    free(t);
    return noIters;
}

#endif // !_NUMERICAL_LIB_
