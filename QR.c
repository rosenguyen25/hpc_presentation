//QR Decomposition using Householder Reflections

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

//Global Variables
double **a; //Matrix whose SVD needs to be found
double *c;
double *d;

//Function
int qrdcmp(double **a, int m, int n, double *c, double *d);

int main (void)
{
    int m =3;
    int n=3;
    int i, j;
    double g[]={12,-51,4,6,167,-68,-4,24,-41};

    a = malloc(sizeof (double*)*m); //allocate M rows to dmatrix
    
    for (i =0; i < m; i++)
    {a[i] = malloc (sizeof (double)*n);} // Now for each row, allocate N actual floats

    for (i =0; i < m; i++)
    {
        for (j =0; j<n; j++)
        {
            a[i][j]=g[i*n+j];
        }
    }
    
    printf("\n");
    printf("A matrix\n");
    for (i =0; i < m; i++)
    {
        for (j =0; j< n; j++)
        {
            printf("%9.4lf\t",a[i][j]);
        }
        printf("\n");
    }
    
    c = malloc(sizeof(double)*n);
    d = malloc(sizeof(double)*n);
    
    for (i=0; i<n; i++)
    {
        c[i] = 0.0;
        d[i] = 0.0;
    }
    
    qrdcmp (a,m,n,c,d);
    double R[m][n];
    for (i =0; i < m; i++){
        for (j=0; j <n; j++){
            if (i == j){
                R[i][i] = d[i];
            }
            else if (i < j){
                R[i][j] = a[i][j];
            }
            else{
                R[i][j] = 0;
            }
        }
    }
    
    printf("\n");
    printf("R decomposition of a matrix:\n");
    for (i =0; i <m; i++){
        for (j =0; j<n; j++){
            printf("%9.4lf\t",R[i][j]);
        }
        printf("\n");
    }

}





int qrdcmp(double **a, int m, int n, double *c, double *d)
/*Constructs the QR decomposition of a[1..n][1..n]. The upper triangular matrix R is returned in the upper triangle of a, except for the diagonal elements of R which are returned in
 d[1..n]. The orthogonal matrix Q is represented as a product of n - 1 Householder matrices
 Q1
 . . . Qn-1
 , where Qj = 1 - uj ? uj /cj . The ith component of uj is zero for i = 1, . . . , j - 1
 while the nonzero components are returned in a[i][j] for i = j, . . . , n. sing returns as
 true (1) if singularity is encountered during the decomposition, but the decomposition is still
 completed in this case; otherwise it returns false (0).*/
{
    int i,j,k;
    double scale,sigma,sum,tau;
    scale = sigma = sum = tau =0.0;
    i=j=k=0;
    
    int min;
    if (m >= n) {min = n;}
    else {min = m;}
    
    for (k =0;k < min;k++) {
        //printf("K =%d \n N = %d\n", k, n);
        
        //for (i=k;i<n;i++) scale=FMAX(scale,fabs(a[i][k]));
        for (i=k;i<m;i++) scale=FMAX(scale,fabs(a[i][k]));
        //printf ("scale %lf\n", scale);
        
        if (scale == 0.0) { //Singular case.
            c[k]=d[k]=0.0;
        }
        
        else{ //Form Qk and Qk ? A.
            //for (i=k;i<n;i++)
            for (i=k;i<m;i++)
            {a[i][k] /= scale;}
            //printf ("a[%d][%d] = %lf\n",i,k, a[i][k]);}
            
            //orgfor (sum=0.0,i=k;i<n;i++)
            for (sum=0.0,i=k;i<m;i++)
            {sum += SQR(a[i][k]);}
            //printf ("Sum = %lf\n",sum);
            
            sigma=SIGN(sqrt(sum),a[k][k]);
            //printf ("sigma = %lf\n", sigma);
            
            a[k][k] += sigma;
            //printf ("a[%d][%d]= %lf\n",k ,k, a[k][k]);
            
            c[k]=sigma*a[k][k];
            //printf ("c[%d] = %lf\n", k,c[k]);
            
            d[k] = -scale*sigma;
            //printf ("d[%d] = %lf\n", k,d[k]);
            
            for (j=k+1;j < n;j++) {
                //printf("j=%d\n",j);
                
                //orgfor (sum=0.0,i=k;i<n;i++)
                for (sum=0.0,i=k;i<m;i++)
                {sum += a[i][k]*a[i][j];
                //printf("sum = %lf, a[%d][%d] = %lf, a[%d][%d]=%lf\n",sum, i, k, a[i][k], i, j, a[i][j]);
                }
                //printf("sum = %lf\n",sum);
                
                tau=sum/c[k];
                //printf("tau= %lf\n",tau);
                
                //org for (i=k;i<n;i++)
                for (i=k;i<m;i++)
                {
                    a[i][j] -= tau*a[i][k];
                    // printf("a[%d][%d] = %lf\n", i, j, a[i][j]);
                }
                
            }
        }
    }
    
    return 0;
}
