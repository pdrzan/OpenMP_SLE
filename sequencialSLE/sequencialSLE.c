#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double det(double **B, int m, int n)
{
    int row_size = m;
    int column_size = n;

    if (row_size != column_size)
    {
        printf("DimensionError: Operation Not Permitted \n");
        exit(1);
    }

    else if (row_size == 1)
        return (B[0][0]);

    else if (row_size == 2)
        return (B[0][0] * B[1][1] - B[1][0] * B[0][1]);

    else
    {
        double **minor;
        int row_minor, column_minor;
        int firstrow_columnindex;
        double sum = 0;
        minor = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++)
            minor[i] = (double *)malloc(n * sizeof(double));

        register int row, column;

        // exclude first row and current column
        for (firstrow_columnindex = 0; firstrow_columnindex < row_size;
             firstrow_columnindex++)
        {

            row_minor = 0;

            for (row = 1; row < row_size; row++)
            {

                column_minor = 0;

                for (column = 0; column < column_size; column++)
                {
                    if (column == firstrow_columnindex)
                        continue;
                    else
                        minor[row_minor][column_minor] = B[row][column];

                    column_minor++;
                }

                row_minor++;
            }

            m = row_minor;
            n = column_minor;

            if (firstrow_columnindex % 2 == 0)
                sum += B[0][firstrow_columnindex] * det(minor, m, n);
            else
                sum -= B[0][firstrow_columnindex] * det(minor, m, n);
        }

        return sum;
    }
}

void generateRandomArray(double *array, int lenArray)
{
    for (int i = 0; i < lenArray; i++)
    {
        array[i] = rand() % (lenArray * lenArray * lenArray);
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

void createLinearIndependentMatrix(double **matrix, int lenMatrix)
{
    double determinant = 0;
    while (determinant == 0)
    {
        for (int i = 0; i < lenMatrix; i++)
        {
            generateRandomArray(matrix[i], lenMatrix);
        }
        determinant = det(matrix, lenMatrix, lenMatrix);
    }
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
        for(int j = lenMatrix - 1; j >= 0; j--)
        {
            matrix[i][j] /= matrix[i][i];
        }
        for (int ii = i + 1; ii < lenMatrix; ii++)
        {
            for(int j = i + 1; j < lenMatrix; j++)
            {
                matrix[ii][j] -= matrix[ii][i] * matrix[i][i];
            }
            y[ii] -= matrix[ii][i] * y[i];
            matrix[ii][i] = 0;
        }
    }
    for (int i = lenMatrix - 1; i >= 0; i--)
    {
        for(int j = lenMatrix - 1; j >= i; j--)
        {
            answer[i] += matrix[i][j] * y[j];
        }
    }
}

void printVariable(char *variable, double value)
{
    printf("%s: %f", variable, value);
}

void printfAnswer(double *answer, int lenAnswer)
{
    printVariable("X", answer[0]);
    for(int i = 1; i < lenAnswer; i++)
    {
        printVariable(" X", answer[i]);
    }
    printf("\n");
}

int main()
{
    double **matrix, *y, *answer;
    int lenMatrix;
    srand(time(0));

    printf("Type how many variables you want in the System of Linear Equations: ");
    scanf("%d", &lenMatrix);

    initializeMatrixMemory(&matrix, lenMatrix);
    initializeArrayMemory(&y, lenMatrix);
    initializeArrayMemory(&answer, lenMatrix);

    createLinearIndependentMatrix(matrix, lenMatrix);
    generateRandomArray(y, lenMatrix);

    // Its necessary to transpose the matrix because the linear independent vetification is done by column orientation.
    // In other words, the columns of the matrix represent the linear independent equations
    transposeMatrix(matrix, lenMatrix);
    printMatrix(matrix, lenMatrix);
    printArray(y, lenMatrix);
    printf("=========================================\n");
    resolveLinearSystem(matrix, y, lenMatrix, answer);
    printMatrix(matrix, lenMatrix);
    printArray(y, lenMatrix);
    printfAnswer(answer, lenMatrix);

    return 0;
}