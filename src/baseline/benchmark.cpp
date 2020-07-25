/*********************************************************
 * File: benchmark.cpp
 * Description: This file benchmarks the baseline implementation
 * of gauss elimination
 * 
 * Author: jfhansen
 * Last modification: 19/07/2020
 *********************************************************/

#include <iostream>
#include <cstddef>
#include <algorithm>
#include <random>

//#include "se_gauss_elim.cpp"

void forward_elimination(double *A, double *b, size_t rows, size_t cols);
void back_substitution(double *A, double *b, size_t rows, size_t cols);
void forward_elimination_pivot(double *A, double *b, size_t rows, size_t cols);
void serial_gauss_elimination(const double *A, const double *b, double *x, size_t rows, size_t cols);
double serial_matrix_inverse(const double *A, double *I, size_t rows, size_t cols);

const int N = 3;
const double A_ge[N*N] = {0,1,2,5,4,3,6,7,6}; //{2,1,-1,-3,-1,2,-2,1,2}; //{3,2,-4,2,3,3,5,-3,1};
const double b_vec[N] = {8,22,38}; //{3,15,14};
const double A_inv[N*N] = {2,-1,0,-1,2,-1,0,-1,2}; 

void print_matrix(const double *M, size_t rows, size_t cols)
{
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
            std::cout << M[i*cols + j] << " ";
        std::cout << std::endl;
    }
}

void benchmark_gauss_elimination()
{
    double *A, *b, *x;
    A = new double[N*N];
    b = new double[N];
    x = new double[N];

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
            A[i*N+j] = A_ge[i*N+j];
        b[i] = b_vec[i];
    }

    /*forward_elimination_pivot(A,b,N,N);
    
    std::cout << "A matrix before back substitution: " << std::endl;
    print_matrix(A,N);

    std::cout << "b vector: " << std::endl;
    for (size_t elem = 0; elem < N; elem++)
        std::cout << b[elem] << std::endl;

    back_substitution(A,b,N,N);

    std::cout << "A matrix after back substitution: " << std::endl;
    print_matrix(A,N);

    std::cout << "b vector: " << std::endl;
    for (size_t elem = 0; elem < N; elem++)
        std::cout << b[elem] << std::endl;
    */

    serial_gauss_elimination(A, b, x, N, N);

    std::cout << "A matrix:" << std::endl;
    print_matrix(A, N, N);
    std::cout << "b vector:" << std::endl;
    print_matrix(b, N, 1);
   
    std::cout << "Result: " << std::endl;
    print_matrix(x, N, 1);
}

void benchmark_matrix_inverse()
{
    double *A, *I;
    A = new double[N*N];
    I = new double[N*N];

    for (size_t i=0; i<N; i++)
        for (size_t j=0; j<N; j++)
        {
            A[i*N+j] = A_inv[i*N+j];
            I[i*N+j] = (i==j) ? 1.0 : 0.0;
        }
    
    double det = serial_matrix_inverse(A, I, N, N);
    /*
    std::cout << "Matrix A: " << std::endl;
    print_matrix(A, N, N);
    std::cout << "Determinant of A: " << det << std::endl;
    std::cout << "Inverse of A: " << std::endl;
    print_matrix(I, N, N);
    */
}

int main()
{
    std::cout << "Gauss Elimination: " << std::endl;
    benchmark_gauss_elimination();
    std::cout << "\n\nMatrix inverse:" << std::endl;
    benchmark_matrix_inverse();
    return 0;
}