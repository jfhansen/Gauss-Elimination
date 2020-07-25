/*********************************************************
 * File: se_gauss_elim_pivot.cpp
 * Description: This file contains the implementation of the
 * Gauss-Elimination forward elimination and back substitution
 * algorithms in C++ using partial pivoting.
 * 
 * Author: jfhansen
 * Last modification: 23/07/2020
 *********************************************************/

#include <iostream>
#include <cstddef>
#include <math.h>

void swap_rows(double *A, double *b, size_t row_a, size_t row_b, size_t cols)
{
    if (row_a != row_b)
    {
        double temp;
        size_t i;
        for (i = 0; i < cols; i++)
        {
            temp = A[row_a*cols+i];
            A[row_a*cols+i] = A[row_b*cols+i];
            A[row_b*cols+i] = temp;
        }
        temp = b[row_a];
        b[row_a] = b[row_b];
        b[row_b] = temp;
    }
}

void forward_elimination_pivot(double *A, double *b, size_t rows, size_t cols)
{
    for (size_t pivot_row = 0; pivot_row < rows-1; pivot_row++)
    {
        // Find the pivot in pivot_row
        size_t col_max = pivot_row;
        for (size_t sub_row = pivot_row+1; sub_row < rows; sub_row++)
            if ( abs( A[col_max*cols+pivot_row] ) < abs( A[sub_row*cols+pivot_row] ) )
                col_max = sub_row;
        // If there is a pivot in pivot_row:
        if ( A[col_max*cols+pivot_row] != 0 )
        {
            swap_rows(A, b, pivot_row, col_max, cols);
            // Do for all rows beneath pivot
            for (size_t sub_row = pivot_row+1; sub_row < rows; sub_row++)
            {
                // Find scale that eliminates column
                double scale = A[sub_row*cols+pivot_row] / A[pivot_row*cols+pivot_row];
                A[sub_row * cols + pivot_row] = 0;
                // Eliminate column in sub_row
                for (size_t col = pivot_row + 1; col < cols; col++)
                    A[sub_row * cols + col] -= scale*A[pivot_row * cols + col];
                b[sub_row] -= scale * b[pivot_row];
            }
        }
    }
}

void back_substitution(double *A, double *b, size_t rows, size_t cols)
{
    for (int i = rows - 1; i >= 0; i--)
    {
        b[i] /= A[i*cols+i];
        A[i*cols+i] /= A[i*cols+i];
        for (int j = i - 1; j >= 0; j--)
        {
            //std::cout << "A[" << i << "," << i << "] = " << A[i*cols+i] << std::endl;
            //std::cout << "A[" << j << "," << i << "] = " << A[j*cols+i] << std::endl;
            double ratio = A[j*cols+i] / A[i*cols+i];
            for (int k = cols - 1; k >= j; k--)
                A[j*cols+k] = A[j*cols+k] - ratio * A[i*cols+k];
            b[j] = b[j] - ratio * b[i];
            //std::cout << "b[" << j << "] = " << b[j] << std::endl;
        }
    }
}

void serial_gauss_elimination_pivot(double *A, double *b, size_t rows, size_t cols)
{

}
