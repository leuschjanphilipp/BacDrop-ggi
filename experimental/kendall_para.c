#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NCELLS 25
#define NGENES 10

double kendalls_tau(int *x, int *y, int n) {
    int concordant = 0, discordant = 0;
    int i, j;

    // Ranking
    int x_ranks[NCELLS], y_ranks[NCELLS];
    for (i = 0; i < n; i++) {
        x_ranks[i] = 1;
        y_ranks[i] = 1;
        for (j = 0; j < n; j++) {
            if (x[j] < x[i] || (x[j] == x[i] && j < i)) {
                x_ranks[i]++;
            }
            if (y[j] < y[i] || (y[j] == y[i] && j < i)) {
                y_ranks[i]++;
            }
        }
    }
    printf("\n");
    for (int k = 0; k < NCELLS; k++) {
    printf("%d ", x_ranks[k]);
    }
    printf("\n");
    for (int k = 0; k < NCELLS; k++) {
    printf("%d ", y_ranks[k]);
    }
    printf("\n");

    // Counting concordant and discordant pairs
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if ((x_ranks[i] < x_ranks[j] && y_ranks[i] < y_ranks[j]) || (x_ranks[i] > x_ranks[j] && y_ranks[i] > y_ranks[j])) {
                concordant++;
            } else if ((x_ranks[i] < x_ranks[j] && y_ranks[i] > y_ranks[j]) || (x_ranks[i] > x_ranks[j] && y_ranks[i] < y_ranks[j])) {
                discordant++;
            }
        }
    }

    // Calculation of Kendall's tau coefficient
    double tau = (double)(concordant - discordant) / sqrt((double)((n * (n - 1)) / 2));

    return tau;
}


double mean(int *arr, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    return sum / n;
}


// Function to calculate the Pearson correlation coefficient between two arrays
double pearson_correlation(int* x, int* y, int n) {
    double mean_x = mean(x, n);
    double mean_y = mean(y, n);
    double numerator = 0, denominator_x = 0, denominator_y = 0;

    for (int i = 0; i < n; i++) {
        numerator += (x[i] - mean_x) * (y[i] - mean_y);
        denominator_x += pow(x[i] - mean_x, 2);
        denominator_y += pow(y[i] - mean_y, 2);
    }

    double correlation = numerator / sqrt(denominator_x * denominator_y);

    return correlation;
}

int main() {

    // Allocate memory for corr matrix
    double** corr = (double**)malloc(NGENES * sizeof(double*));
    for (int i = 0; i < NGENES; i++) {
        corr[i] = (double*)malloc(NGENES * sizeof(double));
    }

    // Allocate memory for data matrix
    int** data = (int**)malloc(NGENES * sizeof(int*));
    for (int i = 0; i < NGENES; i++) {
        data[i] = (int*)malloc(NCELLS * sizeof(int));
    }

    // Read data from file
    FILE *file;
    file = fopen("../data/toy_small.bin", "rb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    for (int i=0; i<NGENES; i++){
        size_t el_read = fread(data[i], sizeof(int), NCELLS, file);
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
    //print data
    printf("\n");
    for (int i = 0; i < NGENES; i++) {
        for (int j = 0; j < NCELLS; j++) {
            printf("%d, ", data[i][j]);
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
    double check = corr[3][0] + 1;
    printf("%f\n %f", check, corr[3][0] + 1);
    
    
    free(corr);
    free(data);

    return EXIT_SUCCESS;
}