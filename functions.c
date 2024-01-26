#include <stdio.h>
#include "functions.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

void swapint(int *before, int *after)
{
    int temp = *before;
    *before = *after;
    *after = temp;
}
void swapdouble(double *before, double *after)
{
    double temp = *before;
    *before = *after;
    *after = temp;
}

//to sort the data first by row in ascending order
void bubbleSortRows(int rows[],int columns[], double data[], int n)    //*Referenced from the source provided*
{
    int temp;
    for (int i = 0; i < n-1; i++)
    {
        for (int j = 0; j < n-i-1; j++)
        {
            if (rows[j] > rows[j+1])  //When swapping the rows, swapping all of the other data along with it
            {
                swapint(&rows[j], &rows[j+1]);
                swapint(&columns[j], &columns[j+1]);
                swapdouble(&data[j], &data[j+1]);
            }
        }
    }
}

//to sort the columns in ascending order within each of the rows
void bubbleSortCols(int row_start,int row_end, int columns[], double data[])   
{
    int temp;
    for (int i = row_start; i < row_end-1; i++)
    {
        for (int j = row_start; j < row_end-1; j++)
        {
            if (columns[j] > columns[j+1])  //When swapping the columns, swapping all of the other data along with it
            {
                swapint(&columns[j], &columns[j+1]);
                swapdouble(&data[j], &data[j+1]);
            }
        }
    }
}


void ReadMMtoCSR(const char *filename, CSRMatrix *matrix, int *not_on_diagonal)
{
    //Opening the file
    FILE *file = fopen(filename, "r");  

    if (file == NULL)    //Handling the error in case the file is not able to be open
    {
        printf("Error opening the file.\n");
    }

    int n_rows, n_cols, n_non_zeros;

    char line[1024];   //character array for the lines of the MM file (contain 1024 characters maximum)
    while(fgets(line, 1024, file) != NULL)   //doing nothing (skipping) the lines that start with a % (comments)
    { 

        if(line[0] != '%')
        {
            sscanf(line, "%d %d %d", &n_rows, &n_cols, &n_non_zeros);;
            break;
        }
        
    }

    int not_on_diagonal_count = 0;

    
    for (int i = 0; i < n_non_zeros; i++) 
    {
        int row, col;
        double value;

        fgets(line, sizeof(line), file);
        sscanf(line, "%d %d %lf", &row, &col, &value);

        // Check if the entry is not on the diagonal
        if (row != col) 
        {
            not_on_diagonal_count++;
        }
    }

    *not_on_diagonal = not_on_diagonal_count;

    rewind(file);

    while(fgets(line, 1024, file) != NULL)   //doing nothing (skipping) the lines that start with a % (comments)
    { 
        if(line[0] != '%')
        {
            break;
        }
    }

    //Reading the rows, columns and number of non-zero elements from the file and storing them in the corresponding members of the CSR matrix structure
    sscanf(line, "%d %d %d", &(*matrix).num_rows, &(*matrix).num_cols, &(*matrix).num_non_zeros);
    
    (*matrix).num_non_zeros = (*matrix).num_non_zeros + not_on_diagonal_count;  //size of nonzeros plus the values on the other side of the symmetrix matrix
    
    //Allocating memory for the arrays needed for the CSR_Matrix structure
    (*matrix).csr_data = (double *)calloc(((*matrix).num_non_zeros), sizeof(double));   
    (*matrix).row_ptr = (int *)calloc((*matrix).num_rows+1, sizeof(int));     //+1 to store the number of non-zero elements in the last index of the row pointers aray
    (*matrix).col_ind = (int *)calloc((*matrix).num_non_zeros, sizeof(int));
    

    int temp_row [(*matrix).num_non_zeros];  
    int temp_col [(*matrix).num_non_zeros];
    double temp_data [(*matrix).num_non_zeros];

    // Initialize arrays to zero
    for (int i = 0; i < (*matrix).num_non_zeros; i++) 
    {
        temp_row[i] = 0;
        temp_col[i] = 0;
        temp_data[i] = 0.0;
    }

    int ref = (*matrix).num_non_zeros - not_on_diagonal_count;   //reference to add the mirrored non-zero entries to the upper triangular part of the matrix
    int count = 0;
    //Storing all of the non-zero elements in temporary arrays to manipulate their order
    for (int i = 0; i < (*matrix).num_non_zeros - not_on_diagonal_count; i++)    //-not_on_diagonal_count to only read the non-zero values that are displayed in the file
    {
        fgets(line, 1024, file);   //Getting the line from the file
        sscanf(line, "%d %d %lf", &temp_row[i], &temp_col[i], &temp_data[i]);

        //Make matrix symmetrical beecause CSR does not store the other symmetrical side of the matrix
        if (temp_row[i] != temp_col[i])
        {
            // If the non-zero entry is not in the diagonal, mirror this entry to the upper triangular portion
            //Filling the mirrored entries with the values from the bottom triangular by copying the data to the missing indices and swapping the row and column at that index
            temp_data[count + ref] = temp_data[i];  
            temp_row [count + ref] = temp_col[i];  
            temp_col [count + ref] = temp_row[i];   
 
            count ++;
    
            //The rows and columns will be sorted below therefore the matrix will be fully symmetric due to the mirroring

            // Ensure that ref doesn't exceed the allocated space
            if ((ref + count) > (*matrix).num_non_zeros) 
            {
                fprintf(stderr, "Error: Array overflow.\n");
                exit(EXIT_FAILURE);
            }
        }        
    }

    //Sorting elements by row number
    bubbleSortRows(temp_row, temp_col, temp_data, (*matrix).num_non_zeros);

    
    //Now that rows are sorted, for every row sort by increasing column number   
    //*ChatGPT helped me in writing the logic for this for loop* 
    int current_row_start = 0;
    int counter = 1;
    (*matrix).row_ptr[0] = 0;   //Hard setting the first row to the first index of the column index
    for (int i = 0; i < (*matrix).num_non_zeros; i++)
    {
        if (i == (*matrix).num_non_zeros - 1 || temp_row[i] != temp_row[i+1])    //If the end or the next row are reached --> sort the columns within that row
        {
            bubbleSortCols(current_row_start, i+1, temp_col, temp_data);
            current_row_start = i+1;       //Updating the starting index for the next row
            (*matrix).row_ptr[counter] = i+1;    //Updating the row pointer to indicate the starting indices for the next row
            counter++;
        }
    }

    //Now that all of the data is sorted by increasing columns in every row, add it to the members of the CSR_Matrix structure
    for (int i = 0; i < (*matrix).num_non_zeros; i++)
    {
        (*matrix).col_ind[i] = temp_col[i]-1;   //changing the indexing back to 0
        (*matrix).csr_data[i] = temp_data[i];
    }
    
    printf("Number of non-zeros: %d\n", (*matrix).num_non_zeros);
    printf("Row pointer: ");
    for (int i = 0; i < (*matrix).num_rows; i++)
    {
        printf("%d ", (*matrix).row_ptr[i]);
    }
    printf("%d\n", (*matrix).num_non_zeros);

    printf("Column Index: ");
    for (int i = 0; i < (*matrix).num_non_zeros; i++)
    {
        printf("%d ", (*matrix).col_ind[i]);
    }
    printf("\n");
    printf("Values: ");
    for (int i = 0; i < (*matrix).num_non_zeros; i++)
    {
        printf("%lf ", (*matrix).csr_data[i]);
    }
    printf("\n");

    
    fclose(file);

}

void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{
    for (int i = 0; i < (*A).num_rows; i++)   //Iterating through every row in the matrix
    {
        double sum = 0.0;   //sum for the current row
        for (int j = (*A).row_ptr[i]; j < (*A).row_ptr[i+1]; j++)   //Iterating through every column in every row using the row poitner
        {
            sum += (*A).csr_data[j]*x[(*A).col_ind[j]];    //Multiplying the corresponding elements of the A CSR matrix and the x vector
        }
        
        y[i] = sum;     //Storing the sum of the products in the corresponding position in the y vector
    }
}


void solverJacobi(const CSRMatrix *A, double *b, double *x, int max_iterations, double tolerance)
{
    double *x_next = (double *)calloc((*A).num_rows, sizeof(double));    //Allocating memory for the x_next vector and initially setting it to 0 with calloc

    for (int iterations = 0; iterations < max_iterations; iterations++)    //Run the method until the max iteration runs
    {
        for (int i = 0; i < (*A).num_rows; i++)   //Approximate the solution for every row
        {
            double sum = 0.0;
            double diagonal_val = 0.0;
            for (int j = (*A).row_ptr[i]; j < (*A).row_ptr[i+1]; j++)   //Iterate through the non-zero values for every row
            {
                if ((*A).col_ind[j] != i)   //If the column is not the same number as the row   
                {
                    sum += (*A).csr_data[j]*x[(*A).col_ind[j]];    //From the formula from the provided source
                }
                else
                {
                    diagonal_val = (*A).csr_data[j];
                }
            }

            x_next[i] = (b[i] - sum)/diagonal_val;    //From the formula from the provided link
            
        }

        //Calculating the error within the values of x to compare to the tolerance after
        double error = 0.0;
        for (int i = 0; i < (*A).num_rows; i++)
        {
            error += fabs(x_next[i] - x[i]);
        }

        //Updating the solution vector with the approximations
        for (int i = 0; i < (*A).num_rows; i++)
        {
            x[i] = x_next[i];
        }

        //Check if the method should terminate due to tolerance reached
        if (error < tolerance)
        {
            printf("Converged after %d iterations \n", iterations+1);
            break;   //not continuing the solver if the tolerance is reached
        }
    }

    free(x_next);
}


void compute_residual(const CSRMatrix *A, const double *x, double *y, const double *b, double *r)
{
    //Calculting A times x with the spmv_csr function that was made before
    //The x vector was calculated with the solver method used before in the main function
    //The result to this operation is stored in the vector y
    
    spmv_csr(A, x, y);

    //Obtaining the difference between the two by subtracting the soltuion b from the calculated solution vector y
    for (int i = 0; i < (*A).num_rows; i++)
    {
        r[i] = y[i]-b[i];
    }

}

double compute_norm(const CSRMatrix *A, const double *r)
{
    double sum = 0.0;
    //The norm of a vector is the square root of the sum of the squares of the entries
    for (int i = 0; i < (*A).num_rows; i++)
    {
        sum += r[i]*r[i];
    }

    return sqrt(sum);
}


