/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include "cblas.h"

double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");

	double* C = (double*)malloc(N * N * sizeof(double));

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans,
			N, N, N, 1.0, B, N, B, N, 0.0, C, N);

	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
			N, N, 1.0, A, N, B, N);

	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
			N, N, 1.0, A, N, B, N);

	cblas_daxpy(N*N, 1, B, 1, C, 1);
        return C;
}