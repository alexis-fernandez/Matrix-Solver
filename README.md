# Sparse Matrix Solver using Jacobi's Method

This project involves the development of a program to read a sparse matrix from an MM file, convert it to Compressed Sparse Row (CSR) format, and then solve for vector x using Jacobi's method. Additionally, the program computes the matrix multiplication A*x=b with the calculated x vector and compares the accuracy by calculating the residual between this result and the given solution vector.

## Overview

The main steps implemented in this project are as follows:

1. **Reading Matrix from MM File**: The program reads the matrix data from the MM file, and due to the CSR format's representation of only the lower triangular portion of a symmetric matrix, an algorithm is implemented to mirror the non-zero entries to the upper triangular portion.

2. **CSR Format Conversion**: The CSR format is constructed with the row pointer array, column indices array, and values array for the matrix, making it an efficient representation for sparse matrices.

3. **Jacobi's Method**: The program utilizes Jacobi's method to solve for the vector x given matrix A and solution b. Jacobi's method is chosen for its quick convergence with symmetric sparse matrices, making it suitable for the provided sets of matrices.

4. **Residual Calculation**: After solving for vector x, the program computes the residual by comparing the result of A*x=b to the provided b vector. The norm of this residual vector is then calculated to measure accuracy.

## Jacobi's Method Formula

The formula used for Jacobi's method was obtained from the following source:
https://byjus.com/maths/jacobian-method/#:~:text=Given%20an%20exact%20approximation%20x,1(k%2B1).

## Efficiency of CSR Format

CSR format is chosen for its efficiency in terms of memory size and faster operations for sparse matrices. The inherent structure of CSR format allows for efficient iteration through rows and columns without dealing with zero values in the matrix. This results in reduced storage requirements and faster matrix multiplication operations. Furthermore, CSR format facilitates easier parallelization of operations compared to full matrix format.

## Known Limitations

1. **b1 ss.mtx**: The solver cannot obtain results for this matrix as it is not symmetric, and Jacobi's method is designed for symmetric matrices.

2. **Large Matrices (e.g., 2cubes sphere.mtx, ACTIVSg70K.mtx, tmt sym.mtx, StocF-1465)**: The solver struggles to converge to a solution due to the size of these matrices. An adaptive Jacobi's method solver with relaxation factors is suggested to enhance stability and convergence speed.

## Usage

1. Clone the repository and ensure all necessary dependencies are installed.

2. Run the program, providing the MM file path as input.

3. Analyze the results, including the computed x vector, residual norm, and accuracy metrics.

## Dependencies

- Download the provided MM files for testing.

## Acknowledgments

This project was part of a foundational C programming course at McMaster University.

For detailed information, refer to the Matrix Solver Report and Results file provided.
