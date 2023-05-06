/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */

void allocate_matrices(int N, double **AB, double **ABA_t, double **C, double **D, double **B_t)
{
	register int count = N * N;
	register int size = sizeof(**AB);
	*AB = calloc(count, size);
	*ABA_t = calloc(count, size);
	*C = calloc(count, size);
	*D = calloc(count, size);
	*B_t = calloc(count, size);
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
	for (register int i = 0; i < N; i++)
	{
		register double *B_t_ptr = B_t + i;
		register double *B_ptr = B + i * N;
		for (register int j = 0; j < N; j++, B_t_ptr += N, B_ptr++)
		{
			*B_t_ptr = *B_ptr;
		}
	}

	/*
		calculate A * B
	*/
	for (register int i = 0; i < N; i++)
	{
		register double *AB_ptr = AB + i * N;
		register double *A_ptr = A + i * N;
		for (register int j = 0; j < N; j++, AB_ptr++)
		{
			register double result = 0;
			register double *B_ptr = B + j;
			for (register int k = i; k < N; k++)
			{
				result += *(A_ptr + k) * *(B_ptr + k * N);
			}
			*AB_ptr = result;
		}
	}

	/*
		calculate A * B * A_t
	*/
	for (register int i = 0; i < N; i++)
	{
		register double *ABA_t_ptr = ABA_t + i * N;
		register double *AB_ptr = AB + i * N;
		for (register int j = 0; j < N; j++, ABA_t_ptr++)
		{
			register double result = 0;
			register double *A_ptr = A + j * N;
			for (register int k = j; k < N; k++)
			{
				result += *(AB_ptr + k) * *(A_ptr + k);
			}
			*ABA_t_ptr = result;
		}
	}

	/*
		calculate B_t * B_t
	*/
	for (register int i = 0; i < N; i++)
	{
		register double *C_ptr = C + i * N;
		register double *B_t_ptr_1 = B_t + i * N;
		for (register int j = 0; j < N; j++, C_ptr++)
		{
			register double result = 0;
			register double *B_t_ptr_2 = B_t + j;
			for (register int k = 0; k < N; k++)
			{
				result += *(B_t_ptr_1 + k) * *(B_t_ptr_2 + k * N);
			}
			*C_ptr = result;
		}
	}

	/*
		calculate A * B * A_t + B_t * B_t
	*/
	for (register int i = 0; i < N; i++)
	{
		register double *C_ptr = C + i * N;
		register double *ABA_t_ptr = ABA_t + i * N;
		register double *D_ptr = D + i * N;
		for (register int j = 0; j < N; j++, C_ptr++, ABA_t_ptr++, D_ptr++)
		{
			*D_ptr = *C_ptr + *ABA_t_ptr;
		}
	}

	free(B_t);
	free(AB);
	free(ABA_t);
	free(C);
	return D;
}
