
//  Cholesky_decomposition


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *cholesky(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
            sqrt(A[i * n + i] - s) :
            (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
    
    return L;
}

void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}

int main() {
    int n = 3;
    double m1[] = {25, 15, -5,
                   15, 18,  0,
                   -5,  0, 11};
    double *c1 = cholesky(m1, n);
    printf("\n"); 
    printf("CHOLESKY DECOMPOSITION (A= L L^*): EXAMPLE 1 of 2\n");
    printf("\n");
    printf("A=\n");
    show_matrix(m1,n);
    printf("Example 1: L=\n");
    show_matrix(c1, n);
    printf("\n");
    free(c1);
    printf("\n"); 
    printf("CHOLESKY DECOMPOSITION (A= L L^*): EXAMPLE 2 of 2\n");
    printf("\n");    
    n = 3;
    double m2[] = {4,  12, -16,
                  12,  37, -43,
                 -16, -43,  98};
    double *c2 = cholesky(m2, n);
    printf("A=\n");
    show_matrix(m1,n);
    printf("Example 2: L=\n");
    show_matrix(c2, n);
    free(c2);
    
    return 0;
}
