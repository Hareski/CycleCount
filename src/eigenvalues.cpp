#include <iostream>
#include <cmath>
#include "../header/eigenvalues.h"

/**
 * Compute the Jacobi eigenvalue iteration
 * @details This function computes the eigenvalues and eigenvectors of areal symmetric matrix, using
 * Rutishauser's modifications of the classical Jacobi rotation method with threshold pivoting
 * @param n the size of the matrix
 * @param A the matrix, which must be square, real and symmetric
 * @param it_max maximum number of iterations
 * @param v matrix of eigenvectors
 * @param d eigenvalues, in descending order
 * @param it_num total number of iterations
 * @param rot_num total number of rotations
 * @author Modified by Areski Guilhem Himeur original code by John Burkardt
 * @copyright GNU LGPLv3
 */
void computeEigenvalues(int n, double *A, int it_max, double *v, double *d, int &it_num, int &rot_num) {
    double *bw;
    double c;
    double g;
    double gapq;
    double h;
    int i;
    int j;
    int k;
    int l;
    int m;
    int p;
    int q;
    double s;
    double t;
    double tau;
    double term;
    double termp;
    double termq;
    double theta;
    double thresh;
    double w;
    double *zw;

    identityDoublePrecisionMatrix(n, v);

    getDiagonalDoublePrecisionMatrix(n, A, d);

    bw = new double[n];
    zw = new double[n];
    for (i = 0; i < n; ++i) {
        bw[i] = d[i];
        zw[i] = 0.0;
    }
    it_num = 0;
    rot_num = 0;

    while (it_num < it_max) {
        it_num++;
        // The convergence threshold is based on the size of the elements in
        // the strict upper triangle of the matrix.
        thresh = 0.0;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < j; ++i) {
                thresh = thresh + A[i + j * n] * A[i + j * n];
            }
        }
        thresh = sqrt(thresh) / (double) (4 * n);
        if (thresh == 0.0) {
            break;
        }

        for (p = 0; p < n; p++) {
            for (q = p + 1; q < n; q++) {
                gapq = 10.0 * fabs(A[p + q * n]);
                termp = gapq + fabs(d[p]);
                termq = gapq + fabs(d[q]);

                if (4 < it_num &&
                    termp == fabs(d[p]) &&
                    termq == fabs(d[q])) {
                    // Annihilate tiny off diagonal elements
                    A[p + q * n] = 0.0;
                } else if (thresh <= fabs(A[p + q * n])) {
                    // Otherwise, apply A rotation
                    h = d[q] - d[p];
                    term = fabs(h) + gapq;

                    if (term == fabs(h)) {
                        t = A[p + q * n] / h;
                    } else {
                        theta = 0.5 * h / A[p + q * n];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }
                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * A[p + q * n];
                    // Accumulate corrections to diagonal elements
                    zw[p] = zw[p] - h;
                    zw[q] = zw[q] + h;
                    d[p] = d[p] - h;
                    d[q] = d[q] + h;
                    A[p + q * n] = 0.0;

                    // Rotate, using only information from the upper triangle of A
                    for (j = 0; j < p; ++j) {
                        g = A[j + p * n];
                        h = A[j + q * n];
                        A[j + p * n] = g - s * (h + g * tau);
                        A[j + q * n] = h + s * (g - h * tau);
                    }
                    for (j = p + 1; j < q; ++j) {
                        g = A[p + j * n];
                        h = A[j + q * n];
                        A[p + j * n] = g - s * (h + g * tau);
                        A[j + q * n] = h + s * (g - h * tau);
                    }
                    for (j = q + 1; j < n; ++j) {
                        g = A[p + j * n];
                        h = A[q + j * n];
                        A[p + j * n] = g - s * (h + g * tau);
                        A[q + j * n] = h + s * (g - h * tau);
                    }

                    // Accumulate information in the eigenvector matrix
                    for (j = 0; j < n; ++j) {
                        g = v[j + p * n];
                        h = v[j + q * n];
                        v[j + p * n] = g - s * (h + g * tau);
                        v[j + q * n] = h + s * (g - h * tau);
                    }
                    rot_num++;
                }
            }
        }
        for (i = 0; i < n; ++i) {
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }

    // Restore upper triangle of input matrix
    for (j = 0; j < n; ++j) {
        for (i = 0; i < j; ++i) {
            A[i + j * n] = A[j + i * n];
        }
    }

    // Ascending sort the eigenvalues and eigenvectors
    for (k = 0; k < n - 1; k++) {
        m = k;
        for (l = k + 1; l < n; l++) {
            if (d[l] < d[m]) {
                m = l;
            }
        }
        if (m != k) {
            t = d[m];
            d[m] = d[k];
            d[k] = t;
            for (i = 0; i < n; ++i) {
                w = v[i + m * n];
                v[i + m * n] = v[i + k * n];
                v[i + k * n] = w;
            }
        }
    }
    delete[] bw;
    delete[] zw;
}

/**
 * Gets the value of the diagonal of an R8MAT.
 * @details An R8MAT is A doubly dimensioned array of R8 values, stored as A vector in column-major order.
 * @param n the number of rows and columns of the matrix
 * @param A the N by N matrix
 * @param v the diagonal entries of the matrix [output]
 * @author Modified by Areski Guilhem Himeur original code by John Burkardt
 * @copyright GNU LGPLv3
 */
void getDiagonalDoublePrecisionMatrix(int n, double *A, double *v) {
    for (int i = 0; i < n; ++i) {
        v[i] = A[i + i * n];
    }
}

/**
 * Sets the square matrix A to the identity.
 * @details An R8MAT is A doubly dimensioned array of R8 values, stored as A vector in column-major order.
 * @param n the order of A
 * @param A  N by N identity matrix
 * @author Modified by Areski Guilhem Himeur original code by John Burkardt
 * @copyright GNU LGPLv3
 */
void identityDoublePrecisionMatrix(int n, double *A) {
    int k = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i == j) {
                A[k] = 1.0;
            } else {
                A[k] = 0.0;
            }
            k++;
        }
    }
}