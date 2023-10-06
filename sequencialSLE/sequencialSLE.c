#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int det(int **B, int m, int n)
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
        int **minor;
        int row_minor, column_minor;
        int firstrow_columnindex;
        int sum = 0;
        minor = (int **)malloc(n * sizeof(int *));
        for (int i = 0; i < n; i++)
            minor[i] = (int *)malloc(n * sizeof(int));

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

void generateRandomArray(int *array, int lenArray)
{
    for (int i = 0; i < lenArray; i++)
    {
        array[i] = rand() % (lenArray * lenArray * lenArray);
    }
}

void printArray(int *array, int n)
{
    printf("[");
    for (int i = 0; i < n - 1; i++)
    {
        printf("%d,\t", array[i]);
    }
    printf("%d]\n", array[n - 1]);
}

void printMatrix(int **matrix, int lenMatrix)
{
    for (int i = 0; i < lenMatrix; i++)
    {
        printArray(matrix[i], lenMatrix);
    }
}

void initializeMatrixMemory(int ***matrix, int lenMatrix)
{
    *matrix = (int **)malloc(lenMatrix * sizeof(int *));
    for (int i = 0; i < lenMatrix; i++)
    {
        (*matrix)[i] = (int *)malloc(lenMatrix * sizeof(int));
    }
}

void verifyLinearIndependence(int **matrix, int lenMatrix)
{
    int determinant = 0;
    while (determinant == 0)
    {   
        for (int i = 0; i < lenMatrix; i++)
        {
            generateRandomArray(matrix[i], lenMatrix);
        }
        determinant = det(matrix, lenMatrix, lenMatrix);
    }
}

int main()
{
    int **matrix, lenMatrix, determinant = 0;
    srand(time(0));

    printf("Type how many variables you want in the System of Linear Equations: ");
    scanf("%d", &lenMatrix);

    initializeMatrixMemory(&matrix, lenMatrix);
    
    verifyLinearIndependence(matrix, lenMatrix);

    printMatrix(matrix, lenMatrix);

    return 0;
}