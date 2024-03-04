#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NCELLS 5964
#define NGENES 3379

static double arithmetic_mean(double* data, int size);
static double mean_of_products(double* data1, double* data2, int size);
static double standard_deviation(double* data, int size);

//--------------------------------------------------------
// FUNCTION pearson_correlation
//--------------------------------------------------------
double pearson_correlation(double* independent, double* dependent, int size)
{
    double rho;

    // covariance
    double independent_mean = arithmetic_mean(independent, size);
    double dependent_mean = arithmetic_mean(dependent, size);
    double products_mean = mean_of_products(independent, dependent, size);
    double covariance = products_mean - (independent_mean * dependent_mean);

    // standard deviations of independent values
    double independent_standard_deviation = standard_deviation(independent, size);

    // standard deviations of dependent values
    double dependent_standard_deviation = standard_deviation(dependent, size);

    // Pearson Correlation Coefficient
    rho = covariance / (independent_standard_deviation * dependent_standard_deviation);

    return rho;
}

//--------------------------------------------------------
// FUNCTION arithmetic_mean
//--------------------------------------------------------
static double arithmetic_mean(double* data, int size)
{
    double total = 0;

    // note that incrementing total is done within the for loop
    for(int i = 0; i < size; total += data[i], i++);

    return total / size;
}

//--------------------------------------------------------
// FUNCTION mean_of_products
//--------------------------------------------------------
static double mean_of_products(double* data1, double* data2, int size)
{
    double total = 0;

    // note that incrementing total is done within the for loop
    for(int i = 0; i < size; total += (data1[i] * data2[i]), i++);

    return total / size;
}

//--------------------------------------------------------
// FUNCTION standard_deviation
//--------------------------------------------------------
static double standard_deviation(double* data, int size)
{
    double squares[size];

    for(int i = 0; i < size; i++)
    {
        squares[i] = pow(data[i], 2);
    }

    double mean_of_squares = arithmetic_mean(squares, size);
    double mean = arithmetic_mean(data, size);
    double square_of_mean = pow(mean, 2);
    double variance = mean_of_squares - square_of_mean;
    double std_dev = sqrt(variance);

    return std_dev;
}

void saveCorrelationMatrix(double** correlationMatrix, const char* filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file for writing\n");
        exit(1);
    }

    for (int i = 0; i < NGENES; i++) {
        size_t el_written = fwrite(correlationMatrix[i], sizeof(double), NGENES, file);
        if (el_written != NGENES) {
            fprintf(stderr, "Error writing data to file\n");
            exit(1);
        }
    }

    fclose(file);
}


int main() {

    // Allocate memory for corr matrix
    double** corr = (double**)malloc(NGENES * sizeof(double*));
    for (int i = 0; i < NGENES; i++) {
        corr[i] = (double*)malloc(NGENES * sizeof(double));
    }

    // Allocate memory for data matrix
    double** data = (double**)malloc(NGENES * sizeof(double*));
    for (int i = 0; i < NGENES; i++) {
        data[i] = (double*)malloc(NCELLS * sizeof(double));
    }

    // Read data from file
    FILE *file;
    file = fopen("../data/ec_filterHVG50_X.bin", "rb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    for (int i=0; i<NGENES; i++){
        size_t el_read = fread(data[i], sizeof(double), NCELLS, file);
        if (el_read != NCELLS) {
            fprintf(stderr, "Error reading data from file\n");
            exit(1);
        }
    }
    fclose(file);

    // Compute correlations
#pragma omp parallel for
    for (int i = 0; i < NGENES; i++) {
        for (int j = i; j < NGENES; j++) {
           corr[i][j] = pearson_correlation(data[i], data[j], NCELLS);
        }
    }


    //print diag corr
    printf("\n");
    for (int i = 0; i < NGENES; i++) {
        printf("%f, ", corr[i][i]);
    }
    printf("\n");
    /*
    //print data
    printf("\n");
    for (int i = 0; i < NGENES; i++) {
        for (int j = 0; j < NCELLS; j++) {
            printf("%f, ", data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    //print corr
    printf("\n");
    for (int i = 0; i < NGENES; i++) {
        for (int j = 0; j < NGENES; j++) {
            printf("%f, ", corr[i][j]);
        }
        printf("\n");
    }
    //print sum of rows
    printf("\n");
    for (int i = 0; i < NGENES; i++) {
        int sum = 0;
        for (int j = 0; j < NCELLS; j++) {
            sum += data[i][j];
        }
        printf("Row: %d, Sum: %d\n", i, sum);
    }
    */

   saveCorrelationMatrix(corr, "../data/ec_filterHVG50_corr.bin");
    
    free(corr);
    free(data);

    return EXIT_SUCCESS;
}