// C++ program to demonstrate working of Gaussian Elimination
// method
//#include<bits/stdc++.h>
//#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//using namespace std;

void createRandomMatrix(int N,double **mat){
	srand(time(NULL));
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N+1; j++){
			mat[i][j] = rand() % 100;
		}
	}
}

// function to reduce matrix to r.e.f. Returns a value to 
// indicate whether matrix is singular or not
int forwardElim(int N,double **mat);

int forwardElimParallel(int N,double **mat);

// function to calculate the values of the unknowns
void backSub(int N,double **mat);

// function to get matrix content
void gaussianElimination(int N,double **mat)
{
	/* reduction into r.e.f. */
	int singular_flag = forwardElim(N, mat);

	/* if matrix is singular */
	if (singular_flag != -1)
	{
		printf("Singular Matrix.\n");

		/* if the RHS of equation corresponding to
		zero row is 0, * system has infinitely
		many solutions, else inconsistent*/
		if (mat[singular_flag][N])
			printf("Inconsistent System.");
		else
			printf("May have infinitely many "
				"solutions.");

		return;
	}

	/* get solution to system and print it using
	backward substitution */
	backSub(N, mat);
}

void gaussianEliminationParallel(int N,double **mat)
{
	/* reduction into r.e.f. */
	int singular_flag = forwardElimParallel(N, mat);

	/* if matrix is singular */
	if (singular_flag != -1)
	{
		printf("Singular Matrix.\n");

		/* if the RHS of equation corresponding to
		zero row is 0, * system has infinitely
		many solutions, else inconsistent*/
		if (mat[singular_flag][N])
			printf("Inconsistent System.");
		else
			printf("May have infinitely many "
				"solutions.");

		return;
	}

	/* get solution to system and print it using
	backward substitution */
	backSub(N, mat);
}

// function for elementary operation of swapping two rows
void swap_row(int N, int i, int j, double **mat)
{
	//printf("Swapped rows %d and %d\n", i, j);

	for (int k=0; k<=N; k++)
	{
		double temp = mat[i][k];
		mat[i][k] = mat[j][k];
		mat[j][k] = temp;
	}
}

// function to print matrix content at any stage
void print(int N,double **mat)
{
	for (int i=0; i<N; i++, printf("\n"))
		for (int j=0; j<=N; j++)
			printf("%lf ", mat[i][j]);

	printf("\n");
}

// function to reduce matrix to r.e.f.
int forwardElim(int N,double **mat)
{
	for (int k=0; k<N; k++)
	{
		// Initialize maximum value and index for pivot
		int i_max = k;
		int v_max = mat[i_max][k];


		/* find greater amplitude for pivot if any */
		for (int i = k+1; i < N; i++)
			if (abs(mat[i][k]) > v_max)
				v_max = mat[i][k], i_max = i;

		/* if a principal diagonal element is zero,
		* it denotes that matrix is singular, and
		* will lead to a division-by-zero later. */
		if (!mat[k][i_max])
			return k; // Matrix is singular

		/* Swap the greatest value row with current row */
		if (i_max != k)
			swap_row(N, k, i_max, mat);


		
			
		for (int i=k+1; i<N; i++){
			/* factor f to set current row kth element to 0,
			* and subsequently remaining kth column to 0 */
			double f = mat[i][k]/mat[k][k];

			/* subtract fth multiple of corresponding kth
			row element*/
		
			for (int j=k+1; j<=N; j++)
				mat[i][j] -= mat[k][j]*f;

			/* filling lower triangular matrix with zeros*/
			mat[i][k] = 0;
		}
		

		//print(mat);	 //for matrix state
	}
	//print(mat);		 //for matrix state
	return -1;
}

// function to reduce matrix to r.e.f.
int forwardElimParallel(int N,double **mat)
{
	for (int k=0; k<N; k++)
	{
		// Initialize maximum value and index for pivot
		int i_max = k;
		int v_max = mat[i_max][k];


		/* find greater amplitude for pivot if any */
		for (int i = k+1; i < N; i++)
			if (abs(mat[i][k]) > v_max)
				v_max = mat[i][k], i_max = i;

		/* if a principal diagonal element is zero,
		* it denotes that matrix is singular, and
		* will lead to a division-by-zero later. */
		if (!mat[k][i_max])
			return k; // Matrix is singular

		/* Swap the greatest value row with current row */
		if (i_max != k)
			swap_row(N, k, i_max, mat);
			
		for (int i=k+1; i<N; i++){
			/* factor f to set current row kth element to 0,
			* and subsequently remaining kth column to 0 */
			double f = mat[i][k]/mat[k][k];

			/* subtract fth multiple of corresponding kth
			row element*/
		
			#pragma omp parallel for
			for (int j=k+1; j<=N; j++)
				mat[i][j] -= mat[k][j]*f;

			/* filling lower triangular matrix with zeros*/
			mat[i][k] = 0;
		}

	}
	return -1;
}

// function to calculate the values of the unknowns
void backSub(int N,double **mat)
{
	double x[N]; // An array to store solution

	/* Start calculating from last equation up to the
	first */
	for (int i = N-1; i >= 0; i--)
	{
		/* start with the RHS of the equation */
		x[i] = mat[i][N];

		/* Initialize j to i+1 since matrix is upper
		triangular*/
		for (int j=i+1; j<N; j++)
		{
			/* subtract all the lhs values
			* except the coefficient of the variable
			* whose value is being calculated */
			x[i] -= mat[i][j]*x[j];
		}

		/* divide the RHS by the coefficient of the
		unknown being calculated */
		x[i] = x[i]/mat[i][i];
	}

	// printf("\nSolution for the system:\n");
	// for (int i=0; i<N; i++)
	// 	printf("%lf\n", x[i]);
}

// Driver program
int main()
{
	int N;
	printf("the number of variables: ");
    scanf("%d", &N);
    printf("%d\n", N);

	/* input matrix */
	double ts, totalTime;
    double **mat;

    mat = (double **)malloc(N * sizeof(double *));
    if (mat == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    for (int i = 0; i < N; i++) {
        mat[i] = (double *)malloc((N + 1) * sizeof(double));
        if (mat[i] == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }
    }

	createRandomMatrix(N, mat);
	//show the mean of 20 iterations
	
	ts = omp_get_wtime();
	gaussianElimination(N, mat);
	totalTime = omp_get_wtime() - ts;
	
	printf("Gaussian Elimination, sequencial method time %f s\n", totalTime );

	//regenerate another unsolved matrix
	createRandomMatrix(N, mat);

	ts = omp_get_wtime();
	gaussianEliminationParallel(N, mat);
	totalTime = omp_get_wtime() - ts;

	printf("Gaussian Elimination, parallel method time %f s\n", totalTime );

	for (int i = 0; i < N; i++) {
        free(mat[i]);
    }
    free(mat);

	return 0;
}
