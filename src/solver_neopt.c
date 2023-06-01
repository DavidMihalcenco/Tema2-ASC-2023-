/*
 *  * Tema 2 ASC
 *   * 2023 Spring
 *    */
#include "utils.h"
#include "cblas.h"


/* 
 *  * Add your BLAS implementation here
 *   */
double* my_solver(int N, double *A, double *B) {
	printf("NEOPT SOLVER\n");
	int i, j, k;
	double *BA = (double*)calloc(N * N, sizeof(double));
	double *BC_t = (double*)calloc(N * N, sizeof(double));
	double *BB_t = (double*)calloc(N * N, sizeof(double));
	double *RS = (double*)calloc(N * N, sizeof(double));

	/* BA = A * B*/
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = i; k < N; k++) {
				BA[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}

	/* BA_t = B * At */
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			for (k = j; k < N; k++)
				BB_t[i * N + j] += BA[i * N + k]
					* A[j * N + k];

	/* BB_t = Bt * Bt */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				BC_t[i * N + j] += B[k * N + i] * B[j * N + k];
			}
		}
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			RS[i * N + j] = BB_t[i * N + j] + BC_t[i * N + j];
		}
	}

	free(BA);
	free(BB_t);
	free(BC_t);

	return RS;
}
