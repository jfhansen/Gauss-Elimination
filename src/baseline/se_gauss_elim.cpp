/*********************************************************
 * File: gauss-elim.cpp
 * Description: This file contains the implementation of the
 * Gauss-Elimination forward elimination and back substitution
 * algorithms in a baseline C++ implementation
 * 
 * Author: jfhansen
 * Last modification: 23/07/2020
 *********************************************************/

#include <iostream>
#include <cstddef>
#include <math.h>

// Function to swap rows, 'row_a' and 'row_b', in matrix 'A', and vector 'b'.
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

void swap_rows_inv(double *A, double *I, size_t row_a, size_t row_b, size_t cols)
{
    if (row_a != row_b)
    {
        double temp;
        for (size_t col = 0; col < cols; col++)
        {
            temp = A[row_a*cols+col];
            A[row_a*cols+col] = A[row_b*cols+col];
            A[row_b*cols+col] = temp;
            temp = I[row_a*cols+col];
            I[row_a*cols+col] = I[row_b*cols+col];
            I[row_b*cols+col] = temp;
        }
    }
}

// Forward elimination algorithm used for solving system of linear equations
void forward_elimination(double *A, double *b, size_t rows, size_t cols)
{
    for (size_t i = 0; i < rows - 1; i++)
    {
        if (A[i*cols+i] == 0)
        {
            size_t row = i;
            while( (A[row*cols+row] == 0) )
            {
                row++;
                if (row == rows)
                {
                    std::cout << "Error: No non-zero elements in column " << i << "." << std::endl;
                    break;
                }
            }
            if (row == rows)
                i++;
            else
            {
                std::cout << "Swapping rows " << i << " and " << row << std::endl;
                swap_rows(A, b, i, row, cols);
            }
        }
        for (size_t j = i + 1; j < rows; j++)
        {
            double ratio = A[j*cols+i] / A[i*cols+i];
            for (size_t k = 0; k < cols; k++)
                A[j*cols+k] = A[j*cols+k] - ratio * A[i*cols+k];
            b[j] = b[j] - ratio * b[i];
        }
    }
}

// Back substitution algorithm used in solving system of linear equations
void back_substitution(double *A, double *b, size_t rows, size_t cols)
{
    for (int i = rows - 1; i >= 0; i--)
    {
        b[i] /= A[i*cols+i];
        A[i*cols+i] /= A[i*cols+i];
        for (int j = i - 1; j >= 0; j--)
        {
            double ratio = A[j*cols+i] / A[i*cols+i];
            for (int k = cols - 1; k >= j; k--)
                A[j*cols+k] = A[j*cols+k] - ratio * A[i*cols+k];
            b[j] = b[j] - ratio * b[i];
        }
    }
}

// Function for solving system of linear equations with coefficient matrix 'A', result vector 'b'
// and unknown vector 'x'.
void serial_gauss_elimination(const double *A, const double *b, double *x, size_t rows, size_t cols)
{
    // Copy A array, so that this is not consumed by the gauss elimination
    double *A_tmp = new double[rows*cols];
    std::copy(A, A+rows*cols, A_tmp);
    // Copy b vector into x, so that result is produced in x.
    std::copy(b, b+rows, x);
    // Perform forward elimination, reducing A to upper triangular form.
    forward_elimination(A_tmp, x, rows, cols);
    // Perform back substitution, to solve system of linear equations.
    back_substitution(A_tmp, x, rows, cols);
}


double forward_elimination_inv(double *A, double *I, size_t rows, size_t cols)
{
    double determinant = 1;
    for (size_t row = 0; row < rows; row++)
    {
        if (A[row*cols+row] == 0)
        {
            size_t row_b = row;
            while (A[row_b*cols+row] == 0)
            {
                row_b++;
                if (row_b == rows)
                    break;
            }
            if (row_b == rows)
                row++;
            else
            {
                swap_rows(A, I, row, row_b, cols);
                determinant *= -1;
            }
        }

        for (size_t sub_row = row+1; sub_row < rows; sub_row++)
        {
            double ratio = A[sub_row*cols+row] / A[row*cols+row];
            if (ratio != 0)
            {
                //determinant /= ratio;
                for (size_t col = 0; col < cols; col++)
                {
                    A[sub_row*cols+col] -= ratio * A[row*cols+col];
                    I[sub_row*cols+col] -= ratio * I[row*cols+col];
                }
            }
        }
    }
    std::cout << "A after forward elimination: " << std::endl;
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
            std::cout << A[i*cols + j] << " ";
        std::cout << std::endl;
    }
    std::cout << "A inverse after forward elimination: " << std::endl;
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
            std::cout << I[i*cols + j] << " ";
        std::cout << std::endl;
    }

    for (size_t row = 0; row < rows; row++)
        determinant *= A[row*cols+row];
    std::cout << "Determinant: " << determinant << std::endl;
    return determinant;
}

void back_substitution_inv(double *A, double *I, size_t rows, size_t cols)
{
    for (int row = rows-1; row >= 0; row--)
    {
        A[row*cols+row] /= A[row*cols+row];
        for (int sub_row = row-1; sub_row >= 0; sub_row--)
        {
            double ratio = A[sub_row*cols+row];
            for (int col = cols-1; col >= 0; col--)
            {
                A[sub_row*cols+col] -= ratio * A[row*cols+col];
                I[sub_row*cols+col] -= ratio * I[row*cols+col];
            }
        }
    }
}

// Function for getting inverse, 'I', of matrix 'A'. Returns determinant of A.
double serial_matrix_inverse(const double *A, double *I, size_t rows, size_t cols)
{
    // Copy A array over in a temporary array, so it is not consumed.
    double *A_tmp = new double(rows*cols);
    std::copy(A, A+rows*cols, A_tmp);
    // Perform forward elimination on A and I, and compute the determinant.
    double det = forward_elimination_inv(A_tmp, I, rows, cols);
    // Perform back substitution algorithm on A and I
    back_substitution_inv(A_tmp, I, rows, cols);

    return det;
}
