#ifndef EIGENVALUES_H
#define EIGENVALUES_H

void computeEigenvalues(int n, double *A, int it_max, double *v, double *d, int &it_num, int &rot_num);

void getDiagonalDoublePrecisionMatrix(int n, double *A, double *v);

void identityDoublePrecisionMatrix(int n, double *A);

#endif