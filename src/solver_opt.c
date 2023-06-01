/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");
	double *restrict BA = (double*)calloc(N * N, sizeof(double));
	double *restrict Bt_Bt = (double*)calloc(N * N, sizeof(double));
	double *restrict BB_t = (double*)calloc(N * N, sizeof(double));
	double *restrict RS = (double*)calloc(N * N, sizeof(double));
	double *B_t = (double*)calloc(N * N, sizeof(double));


	/* B_t = B^t */
	for (register int i = 0; i != N; ++i) {
		register double *B_t_ptr = B_t + i;
		register double *B_ptr = B + i * N;

		for (register int j = 0; j != N; ++j, 
			B_t_ptr += N, ++B_ptr) {
			*B_t_ptr = *B_ptr;
		}
	}

	/* BA = A * B */
	for (register int i = 0; i != N; ++i) {
		register double *BA_ptr = BA + i * N;
		register double *A_copy = A + i * N;

		for (register int j = 0; j != N; ++j, ++BA_ptr) {
			register double sum = 0.0;
			register double *A_ptr = A_copy + i;
			register double *B_t_ptr = B_t + j * N + i;
			for (register int k = i; k != N; ++k, ++A_ptr,
				++B_t_ptr)
				sum += *A_ptr * *B_t_ptr;

			*BA_ptr += sum;
		}
	}

	/* BA_t = BA * At */
	for (register int i = 0; i != N; ++i) {
		register double *BB_ptr = BB_t + i * N;
		register double *BA_copy = BA + i * N;

		for (register int j = 0; j != N; ++j, ++BB_ptr) {
			register double sum = 0.0;

			register double *B_ptr = BA_copy + j;
			register double *A_ptr = A + j * (N + 1);

			for (register int k = j; k < N; ++k, ++B_ptr, ++A_ptr)
				sum += *B_ptr * *A_ptr;

			*BB_ptr = sum;
		}
	}

	/*Bt * Bt*/
	for (register int i = 0; i != N; ++i) {
		for (register int j = 0; j != N; ++j) {
			register double *A_ptr = B_t + i*N;
			register double *B_ptr = B_t + j;
			register double *C_ptr = Bt_Bt + i*N + j;
			register double sum = 0.0;
			for (register int k = 0; k != N; ++k) {
				sum += (*A_ptr) * (*B_ptr);
				A_ptr++;
				B_ptr += N;
			}
			(*C_ptr) = sum;
		}
	}
	/* A * B * At + Bt* Bt */
	for (register int i = 0; i != N; ++i) {
		register double *A_ptr = BB_t + i * N;
		register double *B_ptr = Bt_Bt + i * N;
		register double *C_ptr = RS + i * N;
		for (register int j = 0; j != N; ++j) {
			*C_ptr++ = *A_ptr++ + *B_ptr++;
		}
	}

	free(B_t);
	free(BA);
	free(BB_t);
	free(Bt_Bt);
	return RS;
}
