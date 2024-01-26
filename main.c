#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

/*
Sources:

fgets(): https://www.tutorialspoint.com/c_standard_library/c_function_fgets.htm
sscanf(): https://www.tutorialspoint.com/c_standard_library/c_function_sscanf.htm
fscanf(): https://www.tutorialspoint.com/c_standard_library/c_function_fscanf.htm
strstr(): https://www.tutorialspoint.com/c_standard_library/c_function_strstr.htm
fgetc(): https://www.tutorialspoint.com/c_standard_library/c_function_fgetc.htm
EOF: https://www.geeksforgeeks.org/eof-and-feof-in-c/
Bubble sort: https://www.geeksforgeeks.org/bubble-sort/
Jacobian Method: https://byjus.com/maths/jacobian-method/#:~:text=Given%20an%20exact%20approximation%20x,1(k%2B1).
*/

int main(int argc, char *argv[])
{
    const char *filename = argv[1];

    //Handling invalid inputs
    if (argc != 2) 
    {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    int not_on_diagonal = 0;

    CSRMatrix A;
    ReadMMtoCSR(filename, &A, &not_on_diagonal);
    
    
    // Initializing all the vector b (in Ax=b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    // Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
    }

   
    //Stop criteria for the Jacobi method
    int max_iterations = 10000;    
    double tolerance = 1e-16;

    double *x = (double *)calloc(A.num_rows, sizeof(double));
    double *y = (double *)calloc(A.num_rows, sizeof(double));
    double *r = (double *)calloc(A.num_rows, sizeof(double));
    
    printf("%d\n", not_on_diagonal);
    printf("\nThe matrix name: %s\n", filename);
    printf("The dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("Number of non-zeros (read from file): %d\n", A.num_non_zeros-not_on_diagonal);
    printf("Number of non-zeros (from symmetry): %d\n", A.num_non_zeros);


    //Starting the clock before calculating Ax=b and stoppping it after calculating the result in order to record the time taken
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    //Solving for x
    solverJacobi(&A, b, x, max_iterations, tolerance);  

    //print the CPU time taken
    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time taken to solve Ax=b: %f seconds\n", cpu_time_used);

    //Calculting Ax = b with the calculated x and computing the residual between the given b and the b calculated
    compute_residual(&A, x, y, b, r);
    
    //Calculating the norm of the residual vector
    double norm = compute_norm(&A, r);

    printf("Residual Norm: %.16e\n",norm);

    //Freeing all preallocated memory
    free(A.csr_data);
    free(A.row_ptr);
    free(A.col_ind);
    free(x);
    free(y);
    free(b);
    free(r);    
}
