#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size, int seed)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(seed);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if ((bref[i]-b[i]>0.00000005)|| (b[i]-bref[i]> 0.0000005)) return 0;
    }
    return 1;
}

// function that set the vector [vec] to a matrix [mat]
void vecToMat(double* vec,  int n, double** mat){
    int i, j;

    for (i = 0; i < n; i++)
    {
            for (j = 0; j < n; j++)
            {
                mat[i][j]=vec[i * n + j];
            }
    }
}

//function that generate an identity matrix
void matId(double** tfloat_mat, int int_n){
    int i; /* Variable d'iteration */
    int j; /* Variable d'iteration */
    
    /*Creation de la matrice identité*/
    for( i=0; i< int_n; i++){
        for( j=0; j< int_n; j++){
            /* 1 sur la diagonale*/
            if (i==j){
                tfloat_mat[i][j]=1;
            
            }
            /*0 sinon */
            else{
                tfloat_mat[i][j]=0;
            }
        }    
    }
}

/* Apply the function term by term : Line1 <--k*Line1 */
void dilatation(double**matA, int int_taille, double float_k, int int_ligne){
    int i;/*iteration variable */

    /* Apply the function term by term : Line1 <--k*Line1 */
    for (i=0; i<int_taille; i++){
        
        matA[int_ligne][i]*= float_k;
    }
}

/* Apply the function term by term : Line1 <-- Line1 + k*Line2 */
void transvection (double** matA, int int_size, double float_k,int int_ligne1, int int_ligne2){
    int i;
    /* Apply the function term by term : Line1 <-- Line1 + k*Line2 */
    for (i=0; i<int_size; i++){
        matA[int_ligne1][i]+= float_k*matA[int_ligne2][i];
    }
    
}


// Solve the matricial equation AX=B with the Gaussian elimination
int my_dgesv(int n, double *a,  double *b) {
    
    int i; /* iteration variable */
    int j; /* iteration variable */
    int k; /* iteration variable */
    double coeff; /* transvection coefficient */
    double pivot; /* pivot value in the Gaussian elimination */
    /* matrices initialisation */ 
    double** matA; /* matrix version of the vector a */
    double** matB; /* matrix version of the vector b*/
    double** inv; /* inverse matrix*/
    double** res; /* X matrix, result of the matricial equation AX=B*/


    /* Initialization of the matrices */
    res=malloc(n*sizeof(double*));
    inv=malloc(n*sizeof(double*));
    matA=malloc(n*sizeof(double*));
    matB=malloc(n*sizeof(double*));
    for (i=0; i<n; i++){
        res[i]=malloc(n*sizeof(double));
        inv[i]=malloc(n*sizeof(double));
        matA[i]=malloc(n*sizeof(double));
        matB[i]=malloc(n*sizeof(double));
    }
    // set the initial vectors to matrices
    vecToMat(a, n, matA);
    vecToMat(b, n, matB);
    matId(inv, n);

    /* triangularistion of the matrix A */
    for (i=0; i<n; i++){
        /*Definition of the gaussian pivot */
        pivot=matA[i][i];
        for (j=0; j<n; j++){
            /* perform all the necessary transvections so that every coefficient under the pivot get set to zero*/
             if (i!=j){
                coeff= - matA[j][i]/pivot;
                transvection(matA, n,coeff, j, i);
                transvection(inv, n,coeff, j, i);
                
            }
        }
        
    }
    /* reduction of the diagonal to a diagonal containing only ones using dilatations */
    for (i=0; i<n; i++){
        pivot=matA[i][i];
        coeff= 1/pivot;
        dilatation(matA, n, coeff,i);
        dilatation(inv,n, coeff,i);
        
    }

    /* reduction of the matrix A to identity, meaning inv is set to inv(A) */
    for (i=n-1; i>0; i--){
        for (j=n-1; j>0; j--){
            /* coefficients above the pivot are set to zeros using transvections again */
            if ((i>=j)&&(i!=0)){
                coeff= - matA[j-1][i];
                transvection(matA, n,coeff, j-1, i);
                transvection(inv, n,coeff, j-1, i);
                
            }
        }
    }

    
    // calculate the solution of the equation, then set the result into the vector b

    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            for (k=0;k<n; k++){
                res[i][j]+= inv[i][k]*matB[k][j];
            }
            b[i* n + j]=res[i][j];
        }
    }

    return 0;   
}


    void main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size,1);
        aref = generate_matrix(size,1);        
        b = generate_matrix(size,2);
        bref = generate_matrix(size,2);

        
        //print_matrix("A", a, size);
        //print_matrix("B", b, size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);        
        my_dgesv(n,  a, b);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,b,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
       //print_matrix("X", b, size);
       //print_matrix("Xref", bref, size);
    }
