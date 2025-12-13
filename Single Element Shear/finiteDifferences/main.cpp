#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./libs/Sparse.h"

double E = 68e9; // [Pa]
double nu = 0.33;
double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
double mu = E / (2.0 * (1.0 + nu));

enum boundaryPos_enum {
    LEFT,
    RIGHT,
    TOP,
    BOTTOM,
    CORNER_LT,
    CORNER_RT,
    CORNER_LB,
    CORNER_RB,
    NOT_BOUNDARY
};

enum boundaryType_enum {
    FREE,
    CONST_U,
    CONST_UX,
    CONST_UY,
    CONST_F,
    CONST_FX,
    CONST_FY,
};

int BiCGSTAB(Sparse *A, double *b, double *x, int n);

struct __boundaryGroup__ {
    int num;
    double ux;
    double uy;
    double fx;
    double fy;

    /// @brief Constructor - Sets number of boundary group and
    /// boundary conditions
    /// @param num Group Number
    /// @param ux Displacement Boundary Condition for the x-Axis
    /// @param uy Displacement Boundary Condition for the y-Axis
    /// @param fx External Force Boundary Condition for the x-Axis
    /// @param fy External Force Boundary Condition for the y-Axis
    __boundaryGroup__(int num, double ux, double uy, double fx, double fy) {
        this->num = num;

        this->ux = ux;
        this->uy = uy;
        this->fx = fx;
        this->fy = fy;
    }
} typedef bGroup;

struct __gridP__ {
    // Position
    double x;
    double y;

    int num;

    // State
    double u_x;
    double u_y;
    double eps_x;
    double eps_y;
    double sigma_x;
    double sigma_y;

    // Flags
    int boundaryPos;
    int boundaryCond;
    int boundaryGroup;

    /// @brief Constructor - Sets state to 0, defines position
    /// @param x X Position - i in grid
    /// @param y Y Position - j in grid
    /// @param boundaryPos Which boundary does the point belong to (=NOT_BOUNDARY if inner point)
    /// @param boundaryCond Which boundary condition does the point have (=NOT_BOUNDARY if none)
    /// @param boundaryCond Which boundary group this point belongs to (=0 if none)
    /// @param num Number of point in grid
    __gridP__(double x, double y, int boundaryPos, int boundaryCond, int boundaryGroup, int num) {
        this->x = x;
        this->y = y;

        this->boundaryPos = boundaryPos;
        this->boundaryCond = boundaryCond;
        this->boundaryGroup = boundaryGroup;

        this->u_x = 0;
        this->u_y = 0;
        this->eps_x = 0;
        this->eps_y = 0;
        this->sigma_x = 0;
        this->sigma_y = 0;
    }

    /// @brief Constructor - Sets state to 0,
    /// defines position, assumes non boundary
    /// point with no boundary condition
    /// @param x X Position - i in grid
    /// @param y Y Position - j in grid
    /// @param num Number of point in grid
    __gridP__(double x, double y, int num) {
        this->x = x;
        this->y = y;

        this->boundaryPos = NOT_BOUNDARY;
        this->boundaryCond = NOT_BOUNDARY;
        this->boundaryGroup = 0;

        this->u_x = 0;
        this->u_y = 0;
        this->eps_x = 0;
        this->eps_y = 0;
        this->sigma_x = 0;
        this->sigma_y = 0;
    }
} typedef gridPoint;

int main(int argc, char const *argv[]) {
    FILE *fp_header = fopen("./header.txt", "r+");
    FILE *fp_data = fopen("./data.txt", "r+");
    int Nx, Ny, Npoints;
    double dx, dy;
    int bGroupCount;

    // Get #Points and Nx,Ny
    fscanf(fp_header, "%d\n%d\n%d\n%lf\n%lf\n%d\n", &Npoints, &Nx, &Ny, &dx, &dy, &bGroupCount);
    bGroup *bGroups = (bGroup *) malloc(bGroupCount * sizeof(bGroup));
    for (int i = 0; i < bGroupCount; i++) {
        int num;
        double ux;
        double uy;
        double fx;
        double fy;
        fscanf(fp_header, "%d,%lf,%lf,%lf,%lf\n", &num, &ux, &uy, &fx, &fy);
        bGroups[i] = bGroup(num, ux, uy, fx, fy);
    }

    double dx2 = dx * dx, dy2 = dy * dy;
    double idx = 1.0 / dx, idx2 = idx * idx;
    double idy = 1.0 / dy, idy2 = idy * idy;

    // Grid (G) / Sparse Array Grid (Gs) Setup
    gridPoint *G = (gridPoint *) malloc(Npoints * sizeof(gridPoint));
    Sparse *Gs = (Sparse *) malloc(Ny * sizeof(Sparse));
    for (int i = 0; i < Ny; i++) {
        Gs[i] = {NNZ_INIT, Nx};
    }

    // Move data to Grid (G) / Sparse Array Grid (Gs)
    for (int i = 0; i < Npoints; i++) {
        double x;
        double y;
        int isBoundary;
        int boundaryCond;
        int boundaryGroup;
        fscanf(fp_data, "%lf,%lf,%d,%d,%d\n",
               &x, &y, &isBoundary, &boundaryCond, &boundaryGroup);

        int row = round(fabs(y) / dy);
        int col = round(x / dx);
        G[i] = gridPoint(x, y, isBoundary, boundaryCond, boundaryGroup - 1, i + 1);
        Gs[row].set(col, i + 1);
    }
    fclose(fp_header);
    fclose(fp_data);

    int size = 2 * Npoints;
    Sparse *A = (Sparse *) malloc(size * sizeof(Sparse));
    double *b = (double *) malloc(size * sizeof(double));
    double *u = (double *) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        A[i] = {NNZ_INIT, size};
        u[i] = 1;
        b[i] = 0.0;
    }

    double k1x = (lambda + 2.0 * mu) * idx2;        // (λ+2μ)/Δx^2
    double k1y = (lambda + 2.0 * mu) * idy2;        // (λ+2μ)/Δy^2
    double k2x = mu * idx2;                         // μ/Δx^2
    double k2y = mu * idy2;                         // μ/Δy^2
    double k3xy = 0.25 * (lambda + mu) * idx * idy; // (λ+μ)/4ΔxΔy

    for (int row = 0; row < Ny; row++) {
        for (int col = 0; col < Nx; col++) {
            int ind_ij = Gs[row][col];
            if (!ind_ij) // Not indexed == Not part of Grid
                continue;
            ind_ij--;

            int ind_ip1j = 0.0;
            int ind_ip1jp1 = 0.0;
            int ind_ip1jp2 = 0.0;
            int ind_ip1jm1 = 0.0;
            int ind_ip1jm2 = 0.0;
            int ind_ip2j = 0.0;
            int ind_ip2jp1 = 0.0;
            int ind_ip2jm1 = 0.0;
            int ind_im1j = 0.0;
            int ind_im1jp1 = 0.0;
            int ind_im1jp2 = 0.0;
            int ind_im1jm1 = 0.0;
            int ind_im1jm2 = 0.0;
            int ind_im2j = 0.0;
            int ind_im2jp1 = 0.0;
            int ind_im2jm1 = 0.0;
            int ind_ijp1 = 0.0;
            int ind_ijp2 = 0.0;
            int ind_ijm1 = 0.0;
            int ind_ijm2 = 0.0;

            switch (G[ind_ij].boundaryPos) {
                case LEFT:
                    // Indexing :
                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_ip2j = Gs[row][col + 2] - 1;

                    ind_ijp1 = Gs[row + 1][col] - 1;
                    ind_ijm1 = Gs[row - 1][col] - 1;
                    ind_ip1jp1 = Gs[row + 1][col + 1] - 1;
                    ind_ip1jm1 = Gs[row - 1][col + 1] - 1;
                    ind_ip2jp1 = Gs[row + 1][col + 2] - 1;
                    ind_ip2jm1 = Gs[row - 1][col + 2] - 1;

                    // x-Axis :
                    A[ind_ij].set(ind_ij, k1x - 2.0 * k2y); // K{i,j}

                    A[ind_ij].set(ind_ip1j, -2.0 * k1x); // K{i+1,j}
                    A[ind_ij].set(ind_ip2j, +k1x);       // K{i+2,j}
                    A[ind_ij].set(ind_ijp1, +k2y);       // K{i,j+1}
                    A[ind_ij].set(ind_ijm1, +k2y);       // K{i,j-1}

                    // A[ind_ij].set(ind_ijp1 + Npoints, -2.0 * k3xy);   // K{i,j+1}
                    // A[ind_ij].set(ind_ijm1 + Npoints, +2.0 * k3xy);   // K{i,j-1}
                    // A[ind_ij].set(ind_ip1jp1 + Npoints, +2.0 * k3xy); // K{i+1,j+1}
                    // A[ind_ij].set(ind_ip1jm1 + Npoints, -2.0 * k3xy); // K{i+1,j-1}

                    A[ind_ij].set(ind_ip2jp1 + Npoints, -k3xy);
                    A[ind_ij].set(ind_ip2jm1 + Npoints, +k3xy);
                    A[ind_ij].set(ind_ip1jp1 + Npoints, +4.0 * k3xy);
                    A[ind_ij].set(ind_ip1jm1 + Npoints, -4.0 * k3xy);
                    A[ind_ij].set(ind_ijp1 + Npoints, -3.0 * k3xy);
                    A[ind_ij].set(ind_ijm1 + Npoints, +3.0 * k3xy);

                    // y-Axis :
                    A[ind_ij + Npoints].set(ind_ij + Npoints, -2.0 * k1y + k2x); // K{i,j}

                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, -2.0 * k2x); // K{i+1,j}
                    A[ind_ij + Npoints].set(ind_ip2j + Npoints, +k2x);       // K{i+2,j}
                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, +k1y);       // K{i,j+1}
                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, +k1y);       // K{i,j-1}

                    // A[ind_ij + Npoints].set(ind_ijp1, -2.0 * k3xy);   // K{i,j+1}
                    // A[ind_ij + Npoints].set(ind_ijm1, +2.0 * k3xy);   // K{i,j-1}
                    // A[ind_ij + Npoints].set(ind_ip1jp1, +2.0 * k3xy); // K{i+1,j+1}
                    // A[ind_ij + Npoints].set(ind_ip1jm1, -2.0 * k3xy); // K{i+1,j-1}

                    A[ind_ij + Npoints].set(ind_ip2jp1, -k3xy);
                    A[ind_ij + Npoints].set(ind_ip2jm1, +k3xy);
                    A[ind_ij + Npoints].set(ind_ip1jp1, +4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ip1jm1, -4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ijp1, -3.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ijm1, +3.0 * k3xy);
                    break;
                case RIGHT:
                    // Indexing :
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_im2j = Gs[row][col - 2] - 1;

                    ind_ijp1 = Gs[row + 1][col] - 1;
                    ind_ijm1 = Gs[row - 1][col] - 1;
                    ind_im1jp1 = Gs[row + 1][col - 1] - 1;
                    ind_im1jm1 = Gs[row - 1][col - 1] - 1;
                    ind_im2jp1 = Gs[row + 1][col - 2] - 1;
                    ind_im2jm1 = Gs[row - 1][col - 2] - 1;

                    // x-Axis :
                    A[ind_ij].set(ind_ij, k1x - 2.0 * k2y); // K{i,j}

                    A[ind_ij].set(ind_im1j, -2.0 * k1x); // K{i-1,j}
                    A[ind_ij].set(ind_im2j, +k1x);       // K{i-2,j}
                    A[ind_ij].set(ind_ijp1, +k2y);       // K{i,j+1}
                    A[ind_ij].set(ind_ijm1, +k2y);       // K{i,j-1}

                    // A[ind_ij].set(ind_ijp1 + Npoints, +2.0 * k3xy);   // K{i,j+1}
                    // A[ind_ij].set(ind_ijm1 + Npoints, -2.0 * k3xy);   // K{i,j-1}
                    // A[ind_ij].set(ind_im1jp1 + Npoints, -2.0 * k3xy); // K{i-1,j+1}
                    // A[ind_ij].set(ind_im1jm1 + Npoints, +2.0 * k3xy); // K{i-1,j-1}

                    A[ind_ij].set(ind_im2jp1 + Npoints, +k3xy);
                    A[ind_ij].set(ind_im2jm1 + Npoints, -k3xy);
                    A[ind_ij].set(ind_im1jp1 + Npoints, -4.0 * k3xy);
                    A[ind_ij].set(ind_im1jm1 + Npoints, +4.0 * k3xy);
                    A[ind_ij].set(ind_ijp1 + Npoints, +3.0 * k3xy);
                    A[ind_ij].set(ind_ijm1 + Npoints, -3.0 * k3xy);

                    // y-Axis :
                    A[ind_ij + Npoints].set(ind_ij + Npoints, -2.0 * k1y + k2x); // K{i,j}

                    A[ind_ij + Npoints].set(ind_im1j + Npoints, -2.0 * k2x); // K{i-1,j}
                    A[ind_ij + Npoints].set(ind_im2j + Npoints, +k2x);       // K{i-2,j}
                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, +k1y);       // K{i,j+1}
                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, +k1y);       // K{i,j-1}

                    // A[ind_ij + Npoints].set(ind_ijp1, +2.0 * k3xy);   // K{i,j+1}
                    // A[ind_ij + Npoints].set(ind_ijm1, -2.0 * k3xy);   // K{i,j-1}
                    // A[ind_ij + Npoints].set(ind_im1jp1, -2.0 * k3xy); // K{i-1,j+1}
                    // A[ind_ij + Npoints].set(ind_im1jm1, +2.0 * k3xy); // K{i-1,j-1}

                    A[ind_ij + Npoints].set(ind_im2jp1, +k3xy);
                    A[ind_ij + Npoints].set(ind_im2jm1, -k3xy);
                    A[ind_ij + Npoints].set(ind_im1jp1, -4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_im1jm1, +4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ijp1, +3.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ijm1, -3.0 * k3xy);
                    break;
                case TOP:
                    // Indexing :
                    ind_ijp1 = Gs[row + 1][col] - 1;
                    ind_ijp2 = Gs[row + 2][col] - 1;

                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_ip1jp1 = Gs[row + 1][col + 1] - 1;
                    ind_im1jp1 = Gs[row + 1][col - 1] - 1;
                    ind_ip1jp2 = Gs[row + 2][col + 1] - 1;
                    ind_im1jp2 = Gs[row + 2][col - 1] - 1;

                    // x-Axis :
                    A[ind_ij].set(ind_ij, -2.0 * k1x + k2y); // K{i,j}

                    A[ind_ij].set(ind_ijp1, -2.0 * k2y); // K{i,j+1}
                    A[ind_ij].set(ind_ijp2, +k2y);       // K{i,j+2}
                    A[ind_ij].set(ind_ip1j, +k1x);       // K{i+1,j}
                    A[ind_ij].set(ind_im1j, +k1x);       // K{i-1,j}

                    // A[ind_ij].set(ind_ip1j + Npoints, -2.0 * k3xy);   // K{i+1,j}
                    // A[ind_ij].set(ind_im1j + Npoints, +2.0 * k3xy);   // K{i-1,j}
                    // A[ind_ij].set(ind_ip1jp1 + Npoints, +2.0 * k3xy); // K{i+1,j+1}
                    // A[ind_ij].set(ind_im1jp1 + Npoints, -2.0 * k3xy); // K{i-1,j+1}

                    A[ind_ij].set(ind_ip1jp2 + Npoints, -k3xy);
                    A[ind_ij].set(ind_im1jp2 + Npoints, +k3xy);
                    A[ind_ij].set(ind_ip1jp1 + Npoints, +4.0 * k3xy);
                    A[ind_ij].set(ind_im1jp1 + Npoints, -4.0 * k3xy);
                    A[ind_ij].set(ind_ip1j + Npoints, -3.0 * k3xy);
                    A[ind_ij].set(ind_im1j + Npoints, +3.0 * k3xy);

                    // y-Axis :
                    A[ind_ij + Npoints].set(ind_ij + Npoints, k1y - 2.0 * k2x); // K{i,j}

                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, -2.0 * k1y); // K{i,j+1}
                    A[ind_ij + Npoints].set(ind_ijp2 + Npoints, +k1y);       // K{i,j+2}
                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, +k2x);       // K{i+1,j}
                    A[ind_ij + Npoints].set(ind_im1j + Npoints, +k2x);       // K{i-1,j}

                    // A[ind_ij + Npoints].set(ind_ip1j, -2.0 * k3xy);   // K{i+1,j}
                    // A[ind_ij + Npoints].set(ind_im1j, +2.0 * k3xy);   // K{i-1,j}
                    // A[ind_ij + Npoints].set(ind_ip1jp1, +2.0 * k3xy); // K{i+1,j+1}
                    // A[ind_ij + Npoints].set(ind_im1jp1, -2.0 * k3xy); // K{i-1,j+1}

                    A[ind_ij + Npoints].set(ind_ip1jp2, -k3xy);
                    A[ind_ij + Npoints].set(ind_im1jp2, +k3xy);
                    A[ind_ij + Npoints].set(ind_ip1jp1, +4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_im1jp1, -4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ip1j, -3.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_im1j, +3.0 * k3xy);

                    break;
                case BOTTOM:
                    // Indexing :
                    ind_ijm1 = Gs[row - 1][col] - 1;
                    ind_ijm2 = Gs[row - 2][col] - 1;

                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_ip1jm1 = Gs[row - 1][col + 1] - 1;
                    ind_im1jm1 = Gs[row - 1][col - 1] - 1;
                    ind_ip1jm2 = Gs[row - 2][col + 1] - 1;
                    ind_im1jm2 = Gs[row - 2][col - 1] - 1;

                    // x-Axis :
                    A[ind_ij].set(ind_ij, -2.0 * k1x + k2y); // K{i,j}

                    A[ind_ij].set(ind_ijm1, -2.0 * k2y); // K{i,j-1}
                    A[ind_ij].set(ind_ijm2, +k2y);       // K{i,j-2}
                    A[ind_ij].set(ind_ip1j, +k1x);       // K{i+1,j}
                    A[ind_ij].set(ind_im1j, +k1x);       // K{i-1,j}

                    // A[ind_ij].set(ind_ip1j + Npoints, +2.0 * k3xy);   // K{i+1,j}
                    // A[ind_ij].set(ind_im1j + Npoints, -2.0 * k3xy);   // K{i-1,j}
                    // A[ind_ij].set(ind_ip1jm1 + Npoints, +2.0 * k3xy); // K{i+1,j-1}
                    // A[ind_ij].set(ind_im1jm1 + Npoints, -2.0 * k3xy); // K{i-1,j-1}

                    A[ind_ij].set(ind_ip1jm2 + Npoints, +k3xy);
                    A[ind_ij].set(ind_im1jm2 + Npoints, -k3xy);
                    A[ind_ij].set(ind_ip1jm1 + Npoints, -4.0 * k3xy);
                    A[ind_ij].set(ind_im1jm1 + Npoints, +4.0 * k3xy);
                    A[ind_ij].set(ind_ip1j + Npoints, +3.0 * k3xy);
                    A[ind_ij].set(ind_im1j + Npoints, -3.0 * k3xy);

                    // y-Axis :
                    A[ind_ij + Npoints].set(ind_ij + Npoints, k1y - 2.0 * k2x); // K{i,j}

                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, -2.0 * k1y); // K{i,j-1}
                    A[ind_ij + Npoints].set(ind_ijm2 + Npoints, +k1y);       // K{i,j-2}
                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, +k2x);       // K{i+1,j}
                    A[ind_ij + Npoints].set(ind_im1j + Npoints, +k2x);       // K{i-1,j}

                    // A[ind_ij + Npoints].set(ind_ip1j, +2.0 * k3xy);   // K{i+1,j}
                    // A[ind_ij + Npoints].set(ind_im1j, -2.0 * k3xy);   // K{i-1,j}
                    // A[ind_ij + Npoints].set(ind_ip1jm1, -2.0 * k3xy); // K{i+1,j-1}
                    // A[ind_ij + Npoints].set(ind_im1jm1, +2.0 * k3xy); // K{i-1,j-1}

                    A[ind_ij + Npoints].set(ind_ip1jm2, +k3xy);
                    A[ind_ij + Npoints].set(ind_im1jm2, -k3xy);
                    A[ind_ij + Npoints].set(ind_ip1jm1, -4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_im1jm1, +4.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_ip1j, +3.0 * k3xy);
                    A[ind_ij + Npoints].set(ind_im1j, -3.0 * k3xy);
                    break;
                case CORNER_LT:
                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_ijp1 = Gs[row + 1][col] - 1;

                    A[ind_ij].set(ind_ij, +1.0);
                    A[ind_ij].set(ind_ip1j, -0.5);
                    A[ind_ij].set(ind_ijp1, -0.5);

                    A[ind_ij + Npoints].set(ind_ij + Npoints, +1.0);
                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, -0.5);
                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, -0.5);
                    break;
                case CORNER_RT:
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_ijp1 = Gs[row + 1][col] - 1;

                    A[ind_ij].set(ind_ij, +1.0);
                    A[ind_ij].set(ind_im1j, -0.5);
                    A[ind_ij].set(ind_ijp1, -0.5);

                    A[ind_ij + Npoints].set(ind_ij + Npoints, +1.0);
                    A[ind_ij + Npoints].set(ind_im1j + Npoints, -0.5);
                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, -0.5);
                    break;
                case CORNER_LB:
                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_ijm1 = Gs[row - 1][col] - 1;

                    A[ind_ij].set(ind_ij, +1.0);
                    A[ind_ij].set(ind_ip1j, -0.5);
                    A[ind_ij].set(ind_ijm1, -0.5);

                    A[ind_ij + Npoints].set(ind_ij + Npoints, +1.0);
                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, -0.5);
                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, -0.5);
                    break;
                case CORNER_RB:
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_ijm1 = Gs[row - 1][col] - 1;

                    A[ind_ij].set(ind_ij, +1.0);
                    A[ind_ij].set(ind_im1j, -0.5);
                    A[ind_ij].set(ind_ijm1, -0.5);

                    A[ind_ij + Npoints].set(ind_ij + Npoints, +1.0);
                    A[ind_ij + Npoints].set(ind_im1j + Npoints, -0.5);
                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, -0.5);
                    break;
                case NOT_BOUNDARY: // Inner Points
                    // Indexing :
                    ind_ip1j = Gs[row][col + 1] - 1;
                    ind_im1j = Gs[row][col - 1] - 1;
                    ind_ijp1 = Gs[row + 1][col] - 1;
                    ind_ijm1 = Gs[row - 1][col] - 1;
                    ind_ip1jp1 = Gs[row + 1][col + 1] - 1;
                    ind_ip1jm1 = Gs[row - 1][col + 1] - 1;
                    ind_im1jp1 = Gs[row + 1][col - 1] - 1;
                    ind_im1jm1 = Gs[row - 1][col - 1] - 1;

                    // x-Axis :
                    A[ind_ij].set(ind_ij, -2.0 * (k1x + k2y)); // K{i,j}

                    A[ind_ij].set(ind_ip1j, +k1x); // K{i+1,j}
                    A[ind_ij].set(ind_im1j, +k1x); // K{i-1,j}
                    A[ind_ij].set(ind_ijp1, +k2y); // K{i,j+1}
                    A[ind_ij].set(ind_ijm1, +k2y); // K{i,j-1}

                    A[ind_ij].set(ind_ip1jp1 + Npoints, +k3xy); // K{i+1,j+1}
                    A[ind_ij].set(ind_ip1jm1 + Npoints, -k3xy); // K{i+1,j-1}
                    A[ind_ij].set(ind_im1jp1 + Npoints, -k3xy); // K{i-1,j+1}
                    A[ind_ij].set(ind_im1jm1 + Npoints, +k3xy); // K{i-1,j-1}

                    // y-Axis :
                    A[ind_ij + Npoints].set(ind_ij + Npoints, -2.0 * (k1y + k2x)); // K{i,j}

                    A[ind_ij + Npoints].set(ind_ip1j + Npoints, k2x); // K{i+1,j}
                    A[ind_ij + Npoints].set(ind_im1j + Npoints, k2x); // K{i-1,j}
                    A[ind_ij + Npoints].set(ind_ijp1 + Npoints, k1y); // K{i,j+1}
                    A[ind_ij + Npoints].set(ind_ijm1 + Npoints, k1y); // K{i,j-1}

                    A[ind_ij + Npoints].set(ind_ip1jp1, +k3xy); // K{i+1,j+1}
                    A[ind_ij + Npoints].set(ind_ip1jm1, -k3xy); // K{i+1,j-1}
                    A[ind_ij + Npoints].set(ind_im1jp1, -k3xy); // K{i-1,j+1}
                    A[ind_ij + Npoints].set(ind_im1jm1, +k3xy); // K{i-1,j-1}
                    break;
                default:
                    break;
            }

            switch (G[ind_ij].boundaryCond) {
                case CONST_U:
                    // ux{i,j}
                    A[ind_ij].clear();
                    A[ind_ij].set(ind_ij, 1.0);
                    b[ind_ij] = bGroups[G[ind_ij].boundaryGroup].ux;

                    // uy{i,j}
                    A[ind_ij + Npoints].clear();
                    A[ind_ij + Npoints].set(ind_ij + Npoints, 1.0);
                    b[ind_ij + Npoints] = bGroups[G[ind_ij].boundaryGroup].uy;
                    break;
                case CONST_UX:
                    // ux{i,j}
                    A[ind_ij].clear();
                    A[ind_ij].set(ind_ij, 1.0);
                    b[ind_ij] = bGroups[G[ind_ij].boundaryGroup].ux;
                    break;
                case CONST_UY:
                    // uy{i,j}
                    A[ind_ij + Npoints].clear();
                    A[ind_ij + Npoints].set(ind_ij + Npoints, 1.0);
                    b[ind_ij + Npoints] = bGroups[G[ind_ij].boundaryGroup].uy;
                    break;
                case CONST_F:
                    // Fx{i,j}
                    b[ind_ij] = bGroups[G[ind_ij].boundaryGroup].fx;
                    // Fy{i,j}
                    b[ind_ij + Npoints] = bGroups[G[ind_ij].boundaryGroup].fy;
                    break;
                case CONST_FX:
                    // Fx{i,j}
                    b[ind_ij] = bGroups[G[ind_ij].boundaryGroup].fx;
                    break;
                case CONST_FY:
                    // Fy{i,j}
                    b[ind_ij + Npoints] = bGroups[G[ind_ij].boundaryGroup].fy;
                    break;
                case FREE:
                default:
                    // Fx{i,j}
                    b[ind_ij] = 0;

                    // Fy{i,j}
                    b[ind_ij + Npoints] = 0;
                    break;
            }
        }
    }
    free(bGroups);

    /// Divide by Diagonal Element
    for (int row = 0; row < size; row++) {
        if (A[row][row] != 0) {
            double diag = A[row][row];
            for (int j = 0; j < A[row].n_nz; j++) {
                A[row].set(A[row].c[j], (A[row].v[j] / diag));
            }
            b[row] = b[row] / diag;
        }
    }

    FILE *fp_A = fopen("A_Sparse.txt", "w");
    fprintf(fp_A, "%d\n", size + 1);
    fprintf(fp_A, "Row(int),Col(int),Data(float)\n");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < A[i].n_nz - 1; j++) {
            if (A[i].c[j] != -1) {
                fprintf(fp_A, "%d, %d, %.16lf\n", i + 1, A[i].c[j] + 1, A[i].v[j]);
            }
        }
    }
    fclose(fp_A);

    FILE *fp_b = fopen("b_vector.txt", "w+");
    fprintf(fp_b, "Row(int),Data(float)\n");
    for (int i = 0; i < size; i++) {
        fprintf(fp_b, "%d, %.16lf\n", i + 1, b[i]);
    }
    fclose(fp_b);

    // BiCGSTAB(A, b, u, size);
    // FILE *fp_x = fopen("u.txt", "w+");
    // for (int i = 0; i < size; i++) {
    //     fprintf(fp_x, "%.16lf\n", u[i]);
    // }
    // fclose(fp_x);

    free(A);
    free(G);
    free(Gs);
    return 0;
}

int BiCGSTAB(Sparse *A, double *b, double *x, int n) {
    double rs_err = 1e-5;

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
        for (int j = 0; j < A[i].n_nz; j++) {
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
            for (int j = 0; j < A[i].n_nz; j++) {
                Ap[i] += A[i].v[j] * p0[A[i].c[j]];
            }
            // α = ρ0 / ({r0_hat} * {v}) = ρ0 / ({r0} * {Ap})
            r0_Ap += r0_hat[i] * Ap[i];
        }
        double alpha = rho0_old / r0_Ap;

        // {x} = {x_i-1} + α * p_i-1
        // {r} = {r_i-1} - α * Ap
        rs_new = 0;
        for (int i = 0; i < n; i++) {
            h[i] = x[i] + alpha * p0[i];
            s[i] = r0[i] - alpha * Ap[i];
            rs_new += s[i] * s[i];
        }

        printf("res_norm = %lf\n", sqrt(rs_new));
        if (sqrt(rs_new) <= rs_err) {
            // printf("\nIterations of BiCGSTAB: %d\n", noIters);
            break;
        }

        // {t} = [A] * {r_new}
        // ω = {t} * {s} / ({t} * {t})
        double ts = 0, tt = 0;
        for (int i = 0; i < n; i++) {
            t[i] = 0;
            for (int j = 0; j < A[i].n_nz; j++) {
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

        // printf("res_norm = %.12lf\n", sqrt(rs_new));
        if (sqrt(rs_new) <= rs_err) {
            // printf("\nIterations of BiCGSTAB: %d\n", noIters);
            break;
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