#include "./libs/HexElem.h"

#include <cholmod.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define g -9.81 // Acceleration of Gravity [m/sec^2]

double rho = 2.7e3; // [kg/m^3]
double E = 68e9;    // [Pa]
double nu = 0.33;

/// 3D Elasticity Constitutive Matrix [D]
double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
double mu = E / (2.0 * (1.0 + nu));

/// @brief Calculates the local stiffness matrix [k]
/// for the element elem
/// @tparam hexType Type of Hexahedral Element (hex8/hex20)
/// @param elem Element of type hexType
/// @param k_local Pointer to local stiffness matrix
template <typename hexType>
void localStiff(hexType *elem, double **k_local);

/// @brief Calculates the local force for the given node
/// of the element elem
/// @tparam hexType Type of Hexahedral Element (hex8/hex20)
/// @param elem Element of type hexType
/// @param Fi Force to be distributed
/// @param node Node on which to calculate the distributed force
template <typename hexType>
double localForce(hexType elem, double Fi, int node);

/// @brief Calculates the non-linear local stiffness
/// matrix component [k_NL] for the element elem
/// @tparam hexType Type of Hexahedral Element (hex8/hex20)
/// @param elem Element of type hexType
/// @param k_local_NL Pointer to non-linear local stiffness matrix
template <typename hexType>
double localStiff_NL(hexType *elem, double **k_local_NL);

void addElem(cholmod_triplet *T, int row, int col, double data);

int main(int argc, char const *argv[]) {
    int elemNodes; // Nodes (Order) of chosen Hex Element (hex8,hex20)

    int typeArgPos = 0;
    bool calcStressStrain = false;

    double t_stop = 1;
    double dt = 1;

    if (argc <= 2 || !strcmp(argv[1], "?")) {
        printf("Syntax :\n\t%s --elemType type\n", argv[0]);
        printf("Optional Parameters : \n");
        printf("\t--remesh : Runs the python mesher. If not specified, existing node coordinate and connectivity files are used\n");
        printf("\t--transient [t_stop] [dt] : Runs the transient problem, for the specified time [t_stop], with timestep [dt]\n");
        printf("\t--stress : Calculates the Stress and Strain Vectors for each element\n");
        return 0;
    }

    if (argc > 2) {
        for (int i = 0; i < argc; i++) {
            if (!strcmp(argv[i], "--elemType")) {
                if (!strcmp(argv[i + 1], "hex8")) {
                    elemNodes = linear;
                    typeArgPos = i + 1;
                    i++;
                    continue;
                } else if (!strcmp(argv[i + 1], "hex20")) {
                    elemNodes = serendipity;
                    typeArgPos = i + 1;
                    i++;
                    continue;
                }
            }
            if (!strcmp(argv[i], "--remesh")) {
                /// Mesh Generation
                char cmd[256];
                strcpy(cmd, "C:\\Miniforge3\\python.exe .\\mesher.py ");
                strcat(cmd, argv[typeArgPos]);
                strcat(cmd, " > nul"); // Send .py output to void
                system(cmd);
            }
            if (!strcmp(argv[i], "--transient")) {
                t_stop = atof(argv[i + 1]);
                dt = atof(argv[i + 2]);
                i += 2;
            }
            if (!strcmp(argv[i], "--stress")) {
                calcStressStrain = true;
            }
        }
    }

    int nBodies = 1;
    Body<hex8> *bodyArray = (Body<hex8> *) malloc(nBodies * sizeof(Body<hex8>));
    free(bodyArray); // Valgrind check

    /// Load Mesh
    // Load Node Coordinates
    int nNodes;
    FILE *fp = fopen("nodeCoords.txt", "r+");
    fscanf(fp, "%d\n", &nNodes);

    double *x = (double *) malloc(nNodes * sizeof(double));
    double *y = (double *) malloc(nNodes * sizeof(double));
    double *z = (double *) malloc(nNodes * sizeof(double));
    int *pGroup = (int *) malloc(nNodes * sizeof(int));

    for (int i = 0; i < nNodes; i++) {
        fscanf(fp, "%*d,%lf,%lf,%lf,%d\n",
               &x[i], &y[i], &z[i], &pGroup[i]);
    }
    fclose(fp);

    // Load Element Connectivity
    int nElems;
    int nx; // Number of Elements in x Direction
    int ny; // Number of Elements in x Direction
    int nz; // Number of Elements in x Direction
    fp = fopen("hexaCon.txt", "r+");
    fscanf(fp, "%d,%d,%d,%d,%d\n", &nElems, &elemNodes, &nx, &ny, &nz);

    hex8 *elemArray = (hex8 *) malloc(nElems * sizeof(hex8));
    // hex20 *elemArray = (hex20 *) malloc(nElems * sizeof(hex20));

    for (int i = 0; i < nElems; i++) {
        int nodeIndex[elemNodes];
        double xElem[elemNodes];
        double yElem[elemNodes];
        double zElem[elemNodes];
        int pGroupElem[elemNodes];

        for (int j = 0; j < elemNodes; j++) {
            if (j == elemNodes - 1) {
                fscanf(fp, "%d\n", &nodeIndex[elemNodes - 1]);
            } else {
                fscanf(fp, "%d,", &nodeIndex[j]);
            }
            int ind = nodeIndex[j] - 1;
            xElem[j] = x[ind];
            yElem[j] = y[ind];
            zElem[j] = z[ind];
            pGroupElem[j] = pGroup[ind];
            nodeIndex[j] = nodeIndex[j] - 1;
        }

        elemArray[i] = {i, xElem, yElem, zElem, nodeIndex, pGroupElem};
    }
    fclose(fp);
    free(x);
    free(y);
    free(z);
    free(pGroup);

    int nDOFs = nHexDOFs * nNodes;

    /// DOF Reordering
    int globalDOF = 0;
    int *globalDOFs = (int *) malloc(nNodes * sizeof(int));
    for (int i = 0; i < nNodes; i++)
        globalDOFs[i] = NOT_PARSED;

    for (int elem = 0; elem < nElems; elem++) {
        for (int node = 0; node < elemNodes; node++) {
            if (elemArray[elem].pGroupElem[node] == FREE_DOF) { // First, order the free DOFs
                int ind = elemArray[elem].nodeIndex[node];
                if (globalDOFs[ind] == NOT_PARSED) {
                    globalDOFs[ind] = globalDOF;
                    globalDOF += nHexDOFs;
                }
            }
        }
    }
    int freeDOFs = globalDOF;
    int fixedDOFs = nDOFs - freeDOFs;
    for (int elem = 0; elem < nElems; elem++) {
        for (int node = 0; node < elemNodes; node++) {
            if (elemArray[elem].pGroupElem[node] != FREE_DOF) { // Then, order the fixed DOFs
                int ind = elemArray[elem].nodeIndex[node];
                if (globalDOFs[ind] == NOT_PARSED) {
                    globalDOFs[ind] = globalDOF;
                    globalDOF += nHexDOFs;
                }
            }
        }
    }

    /// Start CHOLMOD
    cholmod_common c;
    cholmod_start(&c);

    /// Global [K] Matrix Definition, f: Free, s: Fixed(Supported)
    // | K_ff | K_fs |
    // | - - -|- - - |
    // | K_sf | K_ss |
    // K_ff * u_f + K_fs * u_s = b_f
    // K_sf * u_f + K_ss * u_s = b_s

    /// Max non-zero elements for each sub-matrix
    size_t nnz_max_ff = freeDOFs * elemNodes * nHexDOFs * elemNodes / 3;
    size_t nnz_max_fs = freeDOFs * fixedDOFs / 3;
    size_t nnz_max_sf = freeDOFs * fixedDOFs / 3;
    size_t nnz_max_ss = fixedDOFs * elemNodes * nHexDOFs * elemNodes / 3;

    // Solution {u} Vector
    cholmod_dense *u_free = cholmod_allocate_dense(freeDOFs, 1, freeDOFs,
                                                   CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);
    cholmod_dense *u_fixed = cholmod_allocate_dense(fixedDOFs, 1, freeDOFs,
                                                    CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);

    // RHS {b} Vector
    cholmod_dense *b_free = cholmod_allocate_dense(freeDOFs, 1, freeDOFs,
                                                   CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);
    cholmod_dense *b_fixed = cholmod_allocate_dense(fixedDOFs, 1, freeDOFs,
                                                    CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);

    /// Boundary Velocities (Fix later)
    double ux_d = 3.0; // [m/s]
    // double uy_d = 0.0; // [m/s]
    // double uz_d = 0.0; // [m/s]

    /// Boundary Accelerations (Fix later)
    // double ux_dd = 0.0; // [m/s^2]
    // double uy_dd = 0.0; // [m/s^2]
    // double uz_dd = 0.0; // [m/s^2]

    /// Define Local Stiffness Matrix [k]_local
    double **k_local = (double **) malloc(nHexDOFs * elemNodes * sizeof(double *));
    double **k_local_NL = (double **) malloc(nHexDOFs * elemNodes * sizeof(double *));
    for (int i = 0; i < nHexDOFs * elemNodes; i++) {
        k_local[i] = (double *) malloc(nHexDOFs * elemNodes * sizeof(double));
        k_local_NL[i] = (double *) malloc(nHexDOFs * elemNodes * sizeof(double));
    }

    /// Export Displacements for Visual
    FILE *fp_u_free = fopen("u_free.txt", "w");
    FILE *fp_u_fixed = fopen("u_fixed.txt", "w");
    FILE *fp_dofVector = fopen("dof_vector.txt", "w");

    // Vector of all DOFs
    fprintf(fp_dofVector, "%d,%d\n", -1, freeDOFs);
    for (int i = 0; i < nNodes; i++) {
        fprintf(fp_dofVector, "%d,%d\n", i, globalDOFs[i]);
    }

    /// Transient Simulation
    // If Static, t_stop = dt = 1 -> |u_dot| = |u|
    for (double t = 0; t <= t_stop; t += dt) {
        // k_ff : [K] Matrix of Free DOFs
        // Dimensions : (# Free DOFs) x (# Free DOFs)
        cholmod_sparse *k_ff;
        // Triplet for Free DOFs
        cholmod_triplet *T_ff = cholmod_allocate_triplet(freeDOFs, freeDOFs, nnz_max_ff, 1,
                                                         CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);

        // k_fs : [K] Matrix of Free/Fixed DOF Dependence
        // Dimensions : (# Free DOFs) x (# Fixed DOFs)
        cholmod_sparse *k_fs;
        // Triplet for Free/Fixed DOFs
        cholmod_triplet *T_fs = cholmod_allocate_triplet(freeDOFs, fixedDOFs, nnz_max_fs, 0,
                                                         CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);

        // k_sf : [K] Matrix of Fixed/Free DOF Dependence
        // Dimensions : (# Fixed DOFs) x (# Free DOFs)
        cholmod_sparse *k_sf;
        // Triplet for Fixed/Free DOFs
        cholmod_triplet *T_sf = cholmod_allocate_triplet(fixedDOFs, freeDOFs, nnz_max_sf, 0,
                                                         CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);

        // k_ss : [K] Matrix of Fixed DOFs
        // Dimensions : (# Fixed DOFs) x (# Free DOFs)
        cholmod_sparse *k_ss;
        // Triplet for Fixed DOFs
        cholmod_triplet *T_ss = cholmod_allocate_triplet(fixedDOFs, fixedDOFs, nnz_max_ss, 1,
                                                         CHOLMOD_DOUBLE + CHOLMOD_REAL, &c);
        // Reset Dense Vectors
        for (int i = 0; i < freeDOFs; i++) {
            ((double *) u_free->x)[i] = 0;
            ((double *) u_fixed->x)[i] = 0;

            ((double *) b_free->x)[i] = 0;
            ((double *) b_fixed->x)[i] = 0;
        }

        /// Construct Global Stiffness Matrix [K]_global
        for (int elem = 0; elem < nElems; elem++) {
            // Clear Local Stiffness Matrix
            for (int i = 0; i < nHexDOFs * elemNodes; i++) {
                for (int j = 0; j < nHexDOFs * elemNodes; j++) {
                    k_local[i][j] = 0;
                    k_local_NL[i][j] = 0;
                }
            }
            // Calculate Local Stiffness Matrix
            localStiff(&elemArray[elem], k_local);
            localStiff_NL(&elemArray[elem], k_local_NL);
            for (int i = 0; i < elemNodes; i++) {
                int row = globalDOFs[elemArray[elem].nodeIndex[i]];
                for (int j = 0; j < elemNodes; j++) {
                    int col = globalDOFs[elemArray[elem].nodeIndex[j]];

                    /// Global [K] Matrix is split into 4 segments
                    /// Pointer k_ptr picks each segment, dependent on position
                    cholmod_triplet *T_ptr = nullptr;

                    // Move Row,Col relative to segment
                    int row_tmp = row, col_tmp = col;
                    if (row < freeDOFs && col < freeDOFs) {
                        T_ptr = T_ff;
                    }
                    if (row < freeDOFs && col >= freeDOFs) {
                        col_tmp = col - freeDOFs;
                        T_ptr = T_fs;
                    }
                    if (row >= freeDOFs && col < freeDOFs) {
                        row_tmp = row - freeDOFs;
                        T_ptr = T_sf;
                        continue; /// (!) Skip this matrix, not needed
                    }
                    if (row >= freeDOFs && col >= freeDOFs) {
                        row_tmp = row - freeDOFs;
                        col_tmp = col - freeDOFs;
                        T_ptr = T_ss;
                        continue; /// (!) Skip this matrix, not needed
                    }

                    /// Add Local Stiffness Matrix to corresponding Global Matrix group
                    addElem(T_ptr, row_tmp + DOF_ux, col_tmp + DOF_ux,
                            k_local[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_ux] +
                                k_local_NL[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_ux]);
                    addElem(T_ptr, row_tmp + DOF_ux, col_tmp + DOF_uy,
                            k_local[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_uy] +
                                k_local_NL[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_uy]);
                    addElem(T_ptr, row_tmp + DOF_ux, col_tmp + DOF_uz,
                            k_local[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_uz] +
                                k_local_NL[(i * nHexDOFs) + DOF_ux][(j * nHexDOFs) + DOF_uz]);

                    addElem(T_ptr, row_tmp + DOF_uy, col_tmp + DOF_ux,
                            k_local[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_ux] +
                                k_local_NL[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_ux]);
                    addElem(T_ptr, row_tmp + DOF_uy, col_tmp + DOF_uy,
                            k_local[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_uy] +
                                k_local_NL[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_uy]);
                    addElem(T_ptr, row_tmp + DOF_uy, col_tmp + DOF_uz,
                            k_local[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_uz] +
                                k_local_NL[(i * nHexDOFs) + DOF_uy][(j * nHexDOFs) + DOF_uz]);

                    addElem(T_ptr, row_tmp + DOF_uz, col_tmp + DOF_ux,
                            k_local[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_ux] +
                                k_local_NL[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_ux]);
                    addElem(T_ptr, row_tmp + DOF_uz, col_tmp + DOF_uy,
                            k_local[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_uy] +
                                k_local_NL[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_uy]);
                    addElem(T_ptr, row_tmp + DOF_uz, col_tmp + DOF_uz,
                            k_local[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_uz] +
                                k_local_NL[(i * nHexDOFs) + DOF_uz][(j * nHexDOFs) + DOF_uz]);
                }

                /// External Forces on Free DOFs
                if (elemArray[elem].pGroupElem[i] == FREE_DOF) {
                    // b_free[row + DOF_uz] += localForce(elemArray[elem], rho * g, i); // Gravity
                }

                /// Fixed DOFs (Boundary Conditions)
                if (elemArray[elem].pGroupElem[i] == 4) {
                    ((double *) u_fixed->x)[row - freeDOFs + DOF_ux] = ux_d * t;
                    ((double *) u_fixed->x)[row - freeDOFs + DOF_uy] = 0;
                    ((double *) u_fixed->x)[row - freeDOFs + DOF_uz] = 0;
                }
            }
        }

        /// Transfer Triplets to Sparse Matrices
        k_ff = cholmod_triplet_to_sparse(T_ff, 1, &c);
        k_fs = cholmod_triplet_to_sparse(T_fs, 1, &c);
        k_sf = cholmod_triplet_to_sparse(T_sf, 1, &c);
        k_ss = cholmod_triplet_to_sparse(T_ss, 1, &c);

        /// Construct Free DOF Equation {RHS}
        /// {RHS} = b_free - K_fs * u_fixed
        double one[2] = {1.0, 0}, m1[2] = {-1.0, 0}; // Basic Scalars
        cholmod_sdmult(k_fs, 0, m1, one, u_fixed, b_free, &c);

        // FILE *f_kff = fopen("k_ff.txt", "w");
        // cholmod_write_sparse(f_kff, k_ff, NULL, NULL, &c);
        // fclose(f_kff);

        // FILE *f_bf = fopen("b.txt", "w");
        // cholmod_write_dense(f_kff, b_free, NULL, &c);
        // fclose(f_bf);

        /// Solve Free DOF Equation
        cholmod_factor *L = cholmod_analyze(k_ff, &c);
        cholmod_factorize(k_ff, L, &c);
        u_free = cholmod_solve(CHOLMOD_A, L, b_free, &c);
        cholmod_free_factor(&L, &c);

        for (int elem = 0; elem < nElems; elem++) {
            for (int i = 0; i < elemNodes; i++) {
                int row = globalDOFs[elemArray[elem].nodeIndex[i]];
                if (row < freeDOFs) {
                    // Assign Displacements to Element
                    elemArray[elem].d[i * nHexDOFs + DOF_ux] = ((double *) u_free->x)[row + DOF_ux];
                    elemArray[elem].d[i * nHexDOFs + DOF_uy] = ((double *) u_free->x)[row + DOF_ux];
                    elemArray[elem].d[i * nHexDOFs + DOF_uz] = ((double *) u_free->x)[row + DOF_ux];

                    // Move Nodes to new positions
                    // elemArray[elem].x[i] += ((double *) u_free->x)[row + DOF_ux];
                    // elemArray[elem].y[i] += ((double *) u_free->x)[row + DOF_uy];
                    // elemArray[elem].z[i] += ((double *) u_free->x)[row + DOF_uz];
                } else {
                    elemArray[elem].d[i * nHexDOFs + DOF_ux] = ((double *) u_fixed->x)[row - freeDOFs + DOF_ux];
                    elemArray[elem].d[i * nHexDOFs + DOF_uy] = ((double *) u_fixed->x)[row - freeDOFs + DOF_ux];
                    elemArray[elem].d[i * nHexDOFs + DOF_uz] = ((double *) u_fixed->x)[row - freeDOFs + DOF_ux];

                    // elemArray[elem].x[i] += ((double *) u_fixed->x)[row - freeDOFs + DOF_ux];
                    // elemArray[elem].y[i] += ((double *) u_fixed->x)[row - freeDOFs + DOF_uy];
                    // elemArray[elem].z[i] += ((double *) u_fixed->x)[row - freeDOFs + DOF_uz];
                }
            }
            if (calcStressStrain) {
                elemArray[elem].calculateStrain();
                elemArray[elem].calculateStress();
            }
        }

        // Free DOFs
        fprintf(fp_u_free, "%2.12lf,", t);
        for (int i = 0; i < freeDOFs - 1; i++) {
            fprintf(fp_u_free, "%.12lf,", ((double *) u_free->x)[i]);
        }
        fprintf(fp_u_free, "%.12lf\n", ((double *) u_free->x)[freeDOFs - 1]);

        // Indices of Fixed DOFs
        fprintf(fp_u_fixed, "%2.12lf,", t);
        for (int i = 0; i < fixedDOFs - 1; i++) {
            fprintf(fp_u_fixed, "%.12lf,", ((double *) u_fixed->x)[i]);
        }
        fprintf(fp_u_fixed, "%.12lf\n", ((double *) u_fixed->x)[fixedDOFs - 1]);

        /// Free Sparse Matrices
        cholmod_free_sparse(&k_ff, &c);
        cholmod_free_sparse(&k_fs, &c);
        cholmod_free_sparse(&k_sf, &c);
        cholmod_free_sparse(&k_ss, &c);

        /// Free Triplets
        cholmod_free_triplet(&T_ff, &c);
        cholmod_free_triplet(&T_fs, &c);
        cholmod_free_triplet(&T_sf, &c);
        cholmod_free_triplet(&T_ss, &c);
    }

    // Close Files
    fclose(fp_u_free);
    fclose(fp_u_fixed);
    fclose(fp_dofVector);

    /// Free Local Stiffness Matrix
    for (int i = 0; i < nHexDOFs * elemNodes; i++) {
        free(k_local[i]);
        free(k_local_NL[i]);
    }
    free(k_local);
    free(k_local_NL);

    /// Export Stress Vectors
    if (calcStressStrain) {
        FILE *fp_stress = fopen("stress.txt", "w");
        for (int elem = 0; elem < nElems; elem++) {
            for (int dir = 0; dir < directions; dir++) {
                fprintf(fp_stress, "%lf,", elemArray[elem].sigma[dir]);
            }
            fprintf(fp_stress, "\n");
        }
        fclose(fp_stress);
    }

    /// Free Dense Matrices
    cholmod_free_dense(&u_free, &c);
    cholmod_free_dense(&u_fixed, &c);
    cholmod_free_dense(&b_free, &c);
    cholmod_free_dense(&b_fixed, &c);

    free(elemArray);
    free(globalDOFs);

    cholmod_print_common("common", &c);
    cholmod_gpu_stats(&c);
    cholmod_finish(&c);
    return 0;
}

void addElem(cholmod_triplet *T, int row, int col, double data) {
    if (T->stype == 1) { // TRI U
        if (row > col)
            return;
    } else if (T->stype == -1) { // TRI L
        if (row < col)
            return;
    }
    ((int *) T->i)[T->nnz] = row;
    ((int *) T->j)[T->nnz] = col;
    ((double *) T->x)[T->nnz] = data;
    T->nnz++;
}

template <typename hexType>
void localStiff(hexType *elem, double **k_local) {
    int elemNodes = elem->getElemNodes();

    /// Gauss Quadrature
    for (int iKsi = 0; iKsi < GQ_POINTS; iKsi++) {
        for (int iEta = 0; iEta < GQ_POINTS; iEta++) {
            for (int iZeta = 0; iZeta < GQ_POINTS; iZeta++) {
                elem->jacobian(points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]);
                elem->getLinearDeformMatrix(points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]);

                for (int i = 0; i < nHexDOFs * elemNodes; i++) {
                    double BT_D[directions];
                    for (int col = 0; col < directions; col++) {
                        BT_D[col] = 0;
                        for (int row = 0; row < directions; row++) {
                            BT_D[col] += elem->B_L[row][i] * D[row][col];
                        }
                    }
                    for (int j = 0; j < nHexDOFs * elemNodes; j++) {
                        double BT_D_B_ij = 0;
                        for (int row = 0; row < directions; row++) {
                            BT_D_B_ij += BT_D[row] * elem->B_L[row][j];
                        }

                        BT_D_B_ij *= elem->detJ * w_GQ[iKsi] * w_GQ[iEta] * w_GQ[iZeta];
                        k_local[i][j] += 0.5 * BT_D_B_ij;
                        k_local[j][i] += 0.5 * BT_D_B_ij;
                    }
                }
            }
        }
    }
}

template <typename hexType>
double localStiff_NL(hexType *elem, double **k_local_NL) {
    int elemNodes = elem->getElemNodes();

    /// Gauss Quadrature
    for (int iKsi = 0; iKsi < GQ_POINTS; iKsi++) {
        for (int iEta = 0; iEta < GQ_POINTS; iEta++) {
            for (int iZeta = 0; iZeta < GQ_POINTS; iZeta++) {
                double weight = w_GQ[iKsi] * w_GQ[iEta] * w_GQ[iZeta];
                elem->jacobian(points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]);
                elem->getNonLinearDeformMatrix(points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]);

                for (int i = 0; i < nHexDOFs * elemNodes; i++) {
                    for (int l = 0; l < nHexDOFs * nHexDOFs; l++) {
                        double sum = 0;
                        for (int k = 0; k < nHexDOFs * nHexDOFs; k++) {
                            if (k / nHexDOFs == l / nHexDOFs)
                                sum += elem->B_NL[k][i] * elem->S2[k % nHexDOFs][l % nHexDOFs];
                        }
                        for (int j = 0; j < nHexDOFs * elemNodes; j++) {
                            k_local_NL[i][j] += sum * elem->B_NL[l][j] * elem->detJ * weight;
                        }
                    }
                }
            }
        }
    }
}

template <typename hexType>
double localForce(hexType elem, double Fi, int node) {
    double F = 0;
    for (int iKsi = 0; iKsi < GQ_POINTS; iKsi++) {
        for (int iEta = 0; iEta < GQ_POINTS; iEta++) {
            for (int iZeta = 0; iZeta < GQ_POINTS; iZeta++) {
                double weight = w_GQ[iKsi] * w_GQ[iEta] * w_GQ[iZeta];
                elem.jacobian(points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]);
                F += Fi * elem.N(node, points_GQ[iKsi], points_GQ[iEta], points_GQ[iZeta]) * fabs(elem.detJ) * weight;
            }
        }
    }
    return F;
}