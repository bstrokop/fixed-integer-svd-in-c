void matrix_comparison(int m, int n, gsl_matrix* a, float A[m][n], char* mtype, float cond)
{
    int i, j;
    double d, diff, g;
    double G[m][n], D[m][n];
    printf("\n\nMatrix %s comparison. Condition number=%6.3f\n", mtype, cond);
    printf("\n               GSL double matrix                               Integer matrix (16-bit)");
    printf("                           Difference matrix\n");
    for (i=0; i<m; i++)  {
        for (j=0; j<n; j++) {
            G[i][j] = gsl_matrix_get(a, i, j);
            float sign = 1.0;
            if ( G[i][j] < 0 ) sign = -1;
            D[i][j] = fabs(G[i][j] - sign * fabs( A[i][j]));
            printf("%10.5f ", G[i][j]);
        }

        printf(" ");
        for (j=0; j<n; j++) {
            printf("%10.5f ", A[i][j]);
        }
        printf("   ");
        for (j=0; j<n; j++) {
            printf("%e ", D[i][j]);
        }
        printf("\n");
    }
}
//----------------------------------------------

int scale_float_vector(int n, const float v[], int_fixed res[n], int pow) {
    int i;
    int_fixed fixed_point_one;

    if (abs(pow) > 15) {
       printf("\n\nPossible programming error. Pow= %d in function scale_float_vector\n", pow);
       return 1; // or die
    }

    fixed_point_one = 1 << abs(pow);

    for (i=0; i< n; i++) {
//      
        res[i] = (int_fixed)(v[i] * fixed_point_one + 0.5);
    }
    return 0;
}

int scale_float_matrix(int p, int q, const float A[p][q], int_fixed RES[p][q], int pow) {
    int i, j;
    int_fixed fixed_point_one;

    fixed_point_one = 1 << abs(pow);
    //printf ("fixed point one=%hi sizeof((res[0][0])=%lu\n", fixed_point_one, sizeof(RES[0][0]));
    for (i=0; i<p; i++) {
        for (j=0; j<q; j++) {
            RES[i][j] = (int_fixed)(A[i][j] * fixed_point_one + 0.5);
            //printf ("%d %d A= %f res=%hi\n ", i,j, A[i][j], RES[i][j]); 
        }
    }
}


int scale_int_matrix(int p, int q, int_fixed A[p][q], int pow) {
    int i, j;

    for (i=0; i<p; i++) {
        for (j=0; j<q; j++) {
            if ( pow > 0 ) {
                A[i][j] = A[i][j] << pow;
            } else {
                A[i][j] = A[i][j] >> abs(pow);
            }
        }
    }

}

int unscale_int_matrix(int p, int q, int_fixed A[p][q], float res[p][q], int pow) {
    int i, j;
    int_fixed fixed_point_one = 1 << abs(pow);

    for (i = 0; i < p; i++) {
         for (j = 0; j < q; j++) {
             res[i][j] = (float)(A[i][j])/fixed_point_one;
         }
    }
}


void mat_trans(int p, int q, const float A[p][q], float result[q][p]) {
    int i, j;

    for (i = 0; i < p; i++ ) {
        for (j = 0; j < q; j++) {
            result[j][i] = A[i][j];
        }
    }
}

void int_fixed_mat_trans(int p, int q, const int_fixed A[p][q], int_fixed result[q][p]) {
    int i, j;

    for (i = 0; i < p; i++) {
        for (j=0; j < q; j++) {
             result[j][i] = A[i][j];
        }
    }

}

float vectors_dot_prod(const float *x, const float *y, size_t n)
{
    float res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

long int int_fixed_vectors_dot_prod(const int_fixed *x, const int_fixed *y, int n) {
    long int res = 0;
    int i;
    for (i=0; i<n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}
    
unsigned long int unsigned_int_fixed_vectors_dot_prod(const int_fixed *x, const int_fixed *y, int n) {
    unsigned long int res = 0;
    int i;
    for (i=0; i<n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}


float vectors_dot_prod2(const float *x, const float *y, int n)
{
    float res = 0.0;
    int i = 0;
    for (; i <= n-4; i+=4)
    {
        res += (x[i] * y[i] +
                x[i+1] * y[i+1] +
                x[i+2] * y[i+2] +
                x[i+3] * y[i+3]);
    }
    for (; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

void matrix_vector_mult(int rows, int cols, const float mat[rows][cols], const float vec[cols], float result[rows])
{ // in matrix form: result = mat * vec;
    int i,j;
    float v[cols];

    for (i = 0; i < rows; i++)
    {
        for ( j=0; j < cols; j++) { v[j] = mat[i][j]; }
        result[i] = vectors_dot_prod(v, vec, cols);
    }
    //for (j=0;j< rows;j++) { printf(" %f", result[j]);}
    //printf("\n");
}


void int_fixed_matrix_vector_mult(int rows, int cols, const int_fixed mat[rows][cols], const int_fixed vec[cols], int_fixed result[rows], int pow)
{
    // in matrix form: result = mat * vec;
    int i, j;
    int_fixed v[cols];
    long int res;

    for (i = 0; i < rows; i++)
    {
        for ( j=0; j < cols; j++) { v[j] = mat[i][j]; }
             res = int_fixed_vectors_dot_prod(v, vec, cols); // v represents one row of a matrix here. Correct!
             result[i] = res >> pow;
    }
    //for (j=0;j< rows;j++) { printf(" %hi", result[j]);}
    //printf("\n");
}


void matmul (int n, int p, int q, const float a[n][p], const float b[p][q], float result[n][q]) {

    int i, j, k;

    // Initializing all elements of result matrix to 0
    for(i=0; i<n; ++i)
        for(j=0; j<q; ++j)
        {
            result[i][j] = 0;
        }

    // Multiplying matrices a and b and
    // storing result in result matrix
    for(i=0; i<n; ++i)
        for(j=0; j<q; ++j)
            for(k=0; k<p; ++k)
            {
                result[i][j]+=a[i][k]*b[k][j];
            }

    // Displaying the result
    printf("\nOutput Matrix (matmul):\n");
    for(i=0; i<n; ++i)
        for(j=0; j<q; ++j)
        {
            printf("%f  ", result[i][j]);
            if(j == q-1)
                printf("\n\n");
        }
    return ;
}

void int_fixed_matmul (int n, int p, int q, const int_fixed  A[n][p], const int_fixed B[p][q], int_fixed result[n][q], int pow) {


    int i, j, k;
    long int res;

    // Initializing all elements of result matrix to 0
    for(i=0; i<n; ++i) {
        for(j=0; j<q; ++j)
        {
            result[i][j] = 0;
        }
    }

    // Multiplying matrices a and b and
    // storing result in result matrix
    for(i=0; i<n; ++i) {
        for(j=0; j<q; ++j) {
            res = 0;
            for(k=0; k<p; ++k)
            {
                res = res + A[i][k]*B[k][j];
                if ( A[i][k] == 0 || B[k][j] == 0) {
                    printf(" int_fixed_matmul A,B= %d %d\n", A[i][k],B[k][j]);
                }
            }

          // 4096 corresponds to unity in Q3.12 format:
          if ( abs(res >> pow) > 32767 ) {
              printf("\nPosiible programming error in INT_FIXED_MATMUL. Scaling problems...\n");
              printf("Resulting ATA matrix elem[%d][%d] has magnitude of %li while the limit is %d\n", i,j, res >> pow, 32767);
              printf("Aborting the program...\n");
              exit(1);
          }
//        scale down the result hoping the sign is preserved (compiler dependent?):
          result[i][j] = (int_fixed)(res >> pow) ;
//        printf("int_fixed_matmul> %d %d %hi\n", i,j,result[i][j] );

         }

    }


    // Displaying the result
    printf("\nOutput Matrix (int_fixed_matmul):\n");
    for(i=0; i<n; ++i)
        for(j=0; j<q; ++j)
        {
            printf("%d  ", result[i][j]);
            if(j == q-1)
                printf("\n\n");
        }
    return ;
}

float dot_product(float v[], float u[], int n){
    int i;
    float result = 0;
    for (i = 0; i < n; i++) result += v[i]*u[i];
    return result;
}

int print_matrix(FILE *f, const gsl_matrix *m)
{
    int status, n = 0;

    for (size_t i = 0; i < m->size1; i++) {
        for (size_t j = 0; j < m->size2; j++) {
//              10f --> g
                 if ((status = fprintf(f, "%10f ", gsl_matrix_get(m, i, j))) < 0)
                        return -1;
                 n += status;
        }

        if ((status = fprintf(f, "\n")) < 0)
             return -1;
             n += status;
        }

    return n;
}

int matrix_fprintf(FILE* fp, const char* format,
size_t m, size_t n, float matrix[m][n]) {

int characters = 0;

for ( size_t i = 0; i < m; ++i ) {
    for ( size_t j = 0; j < n; ++j)
        characters += fprintf(fp, format, matrix[i][j]);
        characters += fprintf(fp, "\n");
}
return characters;
}

int int_fixed_matrix_fprintf(FILE* fp, const char* format,
size_t m, size_t n, int_fixed matrix[m][n]) {

int characters = 0;

for ( size_t i = 0; i < m; ++i ) {
    for ( size_t j = 0; j < n; ++j)
        characters += fprintf(fp, format, matrix[i][j]);
        characters += fprintf(fp, "\n");
}
return characters;
}


