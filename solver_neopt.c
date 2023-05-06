/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
void allocate_matrices(int N, double **AB, double **ABA_t, double **C, double **D, double **B_t)
{
	*AB = calloc(N * N, sizeof(**AB));
	*ABA_t = calloc(N * N, sizeof(**ABA_t));
	*C = calloc(N * N, sizeof(**C));
	*D = calloc(N * N, sizeof(**D));
	*B_t = calloc(N * N, sizeof(**B_t));
}

double *my_solver(int N, double *A, double *B)
{
	double *AB;	   // for A * B
	double *ABA_t; // for A * B * A_t
	double *C;	   // for B_t * B_t
	double *D;	   // for A * B * A_t + B_t * B_t
	double *B_t;   // for transposed B

	allocate_matrices(N, &AB, &ABA_t, &C, &D, &B_t);

	/*
		transpose B
	*/
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			B_t[i * N + j] = B[j * N + i];
		}
	}

	/*
		calculate A * B
	*/
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = i; k < N; k++)
			{
				AB[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}

	/*
		calculate A * B * A_t
	*/
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = j; k < N; k++)
			{
				ABA_t[i * N + j] += AB[i * N + k] * A[j * N + k];
			}
		}
	}

	/*
		calculate B_t * B_t
	*/
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				C[i * N + j] += B_t[i * N + k] * B_t[k * N + j];
			}
		}
	}

	/*
		calculate A * B * A_t + B_t * B_t
	*/
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			D[i * N + j] = C[i * N + j] + ABA_t[i * N + j];
		}
	}

	free(B_t);
	free(AB);
	free(ABA_t);
	free(C);
	return D;
}
