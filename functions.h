#ifndef FUNCTIONS_H
#define FUNCTIONS_H


typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix, int *not_on_diagonal);
void spmv_csr(const CSRMatrix *A, const double *x, double *y);

// ###########################################################
void swapint(int *before, int *after);
void swapdouble(double *before, double *after);
void bubbleSortRows(int rows[],int columns[], double data[], int n);
void bubbleSortCols(int row_start,int row_end, int columns[], double data[]);

void solverJacobi(const CSRMatrix *A, double *b, double *x, int max_iterations, double tolerance);

void compute_residual(const CSRMatrix *A, const double *x, double *y, const double *b, double *r);
double compute_norm(const CSRMatrix *A, const double *r);



#endif
