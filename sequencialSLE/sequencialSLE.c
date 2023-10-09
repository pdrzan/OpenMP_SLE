#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void initializeArrayMemory(double **array, int lenArray)
{
    *array = (double *)malloc(lenArray * sizeof(double));
}

void initializeMatrixMemory(double ***matrix, int lenMatrix)
{
    *matrix = (double **)malloc(lenMatrix * sizeof(double *));
    for (int i = 0; i < lenMatrix; i++)
    {
        initializeArrayMemory(&((*matrix)[i]), lenMatrix);
    }
}

void copyElementsMatrix(double **matrixOrigin, double **matrixDestination, int lenMatrix)
{
    for (int i = 0; i < lenMatrix; i++)
    {
        for (int j = 0; j < lenMatrix; j++)
        {
            matrixDestination[i][j] = matrixOrigin[i][j];
        }
    }
}

double det(double **matrix, int lenMatrix)
{
    long divide=1;
    long result=1;
    int counter=1;
    double **a;

    initializeMatrixMemory(&a, lenMatrix);
    copyElementsMatrix(matrix, a, lenMatrix);

    for (int i = 0; i < lenMatrix - 1; i++)
    {
        int temp1 = a[i][i];
        for (int j = counter++; j < lenMatrix; j++)
        {
            int temp2 = a[j][i];
            for (int k = 0; k < lenMatrix; k++)
            {
                a[j][k] = temp1 * a[j][k] - temp2 * a[i][k];
            }
        }
    }

    for (int i = 0; i < lenMatrix; i++)
    {
        result *= a[i][i];
    }

    for (int i = 0; i < lenMatrix - 1; i++)
    {
        divide *= pow(a[i][i], lenMatrix - i - 1);
    }

    result /= divide;
    free(a);
    return result;
}

void generateRandomArray(double *array, int lenArray)
{
    for (int i = 0; i < lenArray; i++)
    {
        array[i] = rand() % (lenArray * lenArray * lenArray) - rand() % (lenArray * lenArray * lenArray);
    }
}

void printArray(double *array, int n)
{
    printf("[");
    for (int i = 0; i < n - 1; i++)
    {
        printf("%.3f,\t", array[i]);
    }
    printf("%.3f]\n", array[n - 1]);
}

void printMatrix(double **matrix, int lenMatrix)
{
    for (int i = 0; i < lenMatrix; i++)
    {
        printArray(matrix[i], lenMatrix);
    }
}

void createLinearIndependentMatrix(double **matrix, int lenMatrix)
{
    double determinant = 0;
    while (determinant == 0)
    {
        for (int i = 0; i < lenMatrix; i++)
        {
            generateRandomArray(matrix[i], lenMatrix);
        }
        determinant = det(matrix, lenMatrix);
    }
    // printf("Determinant: %f\n", determinant);
}

void transposeMatrix(double **matrix, int lenMatrix)
{
    for (int i = 0; i < lenMatrix; i++)
    {
        for (int j = i + 1; j < lenMatrix; j++)
        {
            double auxInt = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = auxInt;
        }
    }
}

void resolveLinearSystem(double **matrix, double *y, int lenMatrix, double *answer)
{
    for (int i = 0; i < lenMatrix; i++)
    {
        y[i] /= matrix[i][i];
        for (int j = lenMatrix - 1; j >= 0; j--)
        {
            matrix[i][j] /= matrix[i][i];
        }
        for (int ii = i + 1; ii < lenMatrix; ii++)
        {
            for (int j = i + 1; j < lenMatrix; j++)
            {
                matrix[ii][j] -= matrix[ii][i] * matrix[i][j];
            }
            y[ii] -= matrix[ii][i] * y[i];
            matrix[ii][i] = 0;
        }
    }
    for (int i = lenMatrix - 1; i >= 0; i--)
    {
        answer[i] = y[i];
        for (int j = lenMatrix - 1; j > i; j--)
        {
            answer[i] -= matrix[i][j] * answer[j];
        }
    }
}

void printVector(double *answer, int lenAnswer, char variableName)
{
    printf("%c:% 7.3f", variableName, answer[0]);
    for (int i = 1; i < lenAnswer; i++)
    {
        printf("\n  % 7.3f", answer[i]);
    }
    printf("\n");
}

int main()
{
    double **matrixA, *y, *answer, time_execution;
    clock_t begin, end;
    int lenMatrix;
    srand(time(0));

    printf("Type how many variables you want in the System of Linear Equations: ");
    scanf("%d", &lenMatrix);

    initializeMatrixMemory(&matrixA, lenMatrix);
    initializeArrayMemory(&y, lenMatrix);
    initializeArrayMemory(&answer, lenMatrix);

    createLinearIndependentMatrix(matrixA, lenMatrix);
    generateRandomArray(y, lenMatrix);

    // Its necessary to transpose the matrix because the linear independent vetification is done by column orientation.
    // In other words, the columns of the matrix represent the linear independent equations
    transposeMatrix(matrixA, lenMatrix);

    // printMatrix(matrixA, lenMatrix);
    // printVector(y, lenMatrix, 'b');
    // printf("=========================================\n");

    begin = clock();
    resolveLinearSystem(matrixA, y, lenMatrix, answer);
    end = clock();

    // printMatrix(matrixA, lenMatrix);
    // printVector(y, lenMatrix, 'b');
    // printVector(answer, lenMatrix, 'X');

    free(matrixA);
    free(answer);
    free(y);

    printf("Time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

    return 0;
}