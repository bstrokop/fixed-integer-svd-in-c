#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <string.h>
#include "invsqrt.h"
/* Declarations for fixed point stuff */

typedef short signed int int_fixed;
#include "fixed.h"


gsl_rng * r;  /* global generator */
int powerIteration(int m, int n, int maxit, int_fixed amat[m][n], int_fixed U[m][n], int_fixed S[m][n], int_fixed V[m][n]);


int
main (int argc, char* argv[])
{
  const gsl_rng_type * T;
  int i, j;
  double x;
  int m =4, n=4;
  float cond = 2; // default condition number
  int maxit = 100;
  float amat[m][n];
  int_fixed int_fixed_A[m][n];
  float A[m][n], U[m][n], V[m][n], S[m][n];
  char* dummy;

  for (i = 0; i < argc; i++) 
  {
      printf(" argument[%d]: %s\n", i, argv[i]);

      if ( argc > i && strcmp(argv[i], "-c") == 0) {

         double tmp = strtod( argv[i+1], &dummy );

         if ( *dummy != 0 )
         {
             printf(" input was not a valid float %s\n", argv[i+1]);
             exit(2);
         }
         else
         {
             cond = tmp;
         } 

      } else if (argc > i && strcmp(argv[i], "-i") == 0) {
   
         j = strtol( argv[i+1], &dummy, 10);

         if ( *dummy != 0 )
         {
             printf(" input was not a valid integer %s\n", argv[i+1]);
             exit(2);
         }
         else
         {
             maxit = j;
         }


      }
  }

  if ( cond == 2.0 ) printf ("Failed to read in the matrix condition number. Using the default.\n");
  printf(">>>Condition number=%f max itns=%d \n", cond, maxit);

  gsl_matrix *a = gsl_matrix_alloc(m,n);
  gsl_matrix *b = gsl_matrix_alloc(m,n);
  gsl_matrix *u = gsl_matrix_alloc(m,n);
  gsl_matrix *v = gsl_matrix_alloc(m,n);
  gsl_matrix *c = gsl_matrix_alloc(m,n);
  gsl_matrix *d = gsl_matrix_alloc(m,n);
  gsl_matrix *smat = gsl_matrix_alloc(m,n);

  gsl_vector *s = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  printf ("generator type: %s\n", gsl_rng_name (r));
  printf ("seed = %lu\n", gsl_rng_default_seed);
  for ( i = 0; i < 4; i++) {
      for ( j = 0; j < 4; j++ ) {
          gsl_matrix_set(a,i,j,gsl_rng_uniform(r));
          //gsl_matrix_set(b,i,j) = gsl_matrix_get(a, i, j);
          //b[i][j] = a[i][j]; // matrix a gonna be destoryed
      }
  }
// Matrix A gonna be destroyed:
  gsl_matrix_memcpy(u,a);
  printf("Matrix A (Monte-Carlo)\n");
  print_matrix(stdout, a);
  
  gsl_linalg_SV_decomp(u, v, s, w);
  printf("Matrix U...\n");
  print_matrix(stdout, u);
  gsl_matrix_set_zero(smat);

  for ( i = 0; i < 4; i++) {
       double temp = gsl_vector_get(s, i);
       gsl_matrix_set(smat, i, i, temp);
//     printf ("s_%d = %g\n", i, gsl_vector_get (s, i));
  }

  printf("Matrix S\n");
  print_matrix(stdout,smat);
  printf("Matrix V\n");
  print_matrix(stdout, v);
// Need to multiply all three matrices:
// FIXME
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, u, smat,
                    0.0, c);
  printf("UxS matrix product:\n");
  print_matrix(stdout, c);
  

  gsl_matrix_transpose_memcpy(b,v);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, c, b,
                    0.0, d);
  printf("UxSxVTranspose matrix product:\n");
  print_matrix(stdout, d);
  
  x = ( gsl_vector_get(s,0) - 
    cond * gsl_vector_get(s,3) ) / ( cond - 1.0f);

  gsl_matrix_set_zero(smat);
  for ( i = 0; i < 4; i++ ) {
//      Maximal singular value will be unity:
       double temp = ( gsl_vector_get(s,i) + x ) / 
                     ( gsl_vector_get(s,0) + x ) ;
       gsl_matrix_set(smat, i, i, temp);
  } 

  printf("New Matrix S\n");
  print_matrix(stdout,smat);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, u, smat,
                    0.0, c);
  printf("New UxS matrix product:\n");
  print_matrix(stdout, c);

// v still untransposed:
  gsl_matrix_transpose_memcpy(b,v);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                    1.0, c, b,
                    0.0, d);

  printf("Final UxSxVTranspose matrix product:\n");
  print_matrix(stdout, d);

// make normal floating array:
  for (i=0; i<4; i++) {
      for(j=0; j<4; j++) {
          amat[i][j] = gsl_matrix_get(d, i, j);
      }
  }
  printf("\nFloating matrix:\n");
  matrix_fprintf(stdout, " %10.6f", m, n, amat);
  scale_float_matrix(m, n, amat, int_fixed_A, 12);              // to get Q3.12 format

// Interface to integer computing:
  int_fixed UU[m][n], SS[m][n], VV[m][n];

  printf("Initial integer matrix (int_fixed_A):\n");
  int_fixed_matrix_fprintf(stdout, " %6hi", m, n, int_fixed_A);
  powerIteration(m, n, maxit, int_fixed_A, UU, SS, VV);

  printf("USV matrices in integer form:\n");
  printf("\nmatrix U:\n");
  int_fixed_matrix_fprintf(stdout, " %6hi", m, n, UU);
  printf("\nmatrix S:\n");
  int_fixed_matrix_fprintf(stdout, " %6hi", m, n, SS);
  printf("\nmatrix V:\n");
  int_fixed_matrix_fprintf(stdout, " %6hi", m, n, VV);

//  exit(0);
!
  printf("\n\nRescaling U, S and V matrices back to float!\n");
  int k;
  for (k=0; k<n; k++) {
      S[k][k]= SS[k][k]/4096.0;
      for (i=0;i<n;i++) {
           V[i][k] = VV[i][k]/32768.0;
	   U[i][k] = UU[i][k]/32768.0;
	   if ( k != i ) S[k][i]=0;
      }

  }

  matrix_comparison(m, n, u, U, "U", cond); 
  matrix_comparison(m, n, smat, S, "S", cond); 
  matrix_comparison(m, n, v, V, "V", cond); 
  

  gsl_matrix_free(a); 
// free more here
  gsl_rng_free (r);
  return 0;
}

int powerIteration(int m, int n, int maxit, int_fixed A[m][n], int_fixed U[m][n], int_fixed S[m][n], int_fixed V[m][n]) {
//    float a[m][n];
    int_fixed v[n], v_old[n], av[n];
    int i, j, k;
    int p, q;
    //int_fixed A[m][n];
    int_fixed AT[m][n];
    float B[m][n];
    float AAT[m][n];
    int_fixed ATA[m][n];
    long int lambda, oldlambda;
    unsigned long int r;
    long int diff;
    float temp[n], atemp[n]; // to get U matrix
    int scale[n+1];

    //printf("\nFloating matrix confirmation:\n");
    //matrix_fprintf(stdout, " %10.6f", m, n, A_ref);

    //  printf(" sizeof A[0,0]=%lu\n", sizeof(A[0][0]));
    /*for(i=0;i<n;i++) {
       for (j=0;j<n;j++)
       printf("i=%d j=%d A=%hi\n", i,j, A[i][j]);
    }*/

    int_fixed_mat_trans(m, n, A, AT);                    // Q 3.12
    //printf("\n After my mattrans %d %d\n", m,n);

//  Need a routine for int_fixed matrices:
/*    for(i=0;i<n;i++) {
       for (j=0;j<n;j++)
           printf("i=%d j=%d AT=%hi\n", i,j, AT[i][j]);
    } */

    //exit(0);
    //matrix_fprintf(stdout, " %10.6f", m, n, AT);

// was 12 here FIXME  !!!!!
    int_fixed_matmul (n, n, n, AT, A, ATA, 12) ;    // back to Q3.12
//    matrix_fprintf(stdout, " %10.6f", m, n, ATA);
//FIXME
//    exit(0);
//  Need a routine to set matrix to const value
    for (i=0; i<n; i++) {
        for (k=0; k<n; k++) {
            S[i][k] = 0; // zero in any format
            //SS[i][k] = 0; // zero in any format
        }
    } 
    scale[0] = 1; // first matrix will not be scaled

//  Take care of ATA:
    int itns = 50;
// fixme
    for (k = 0; k < n; k++) {
        
//      Always start with equal components for the vector:
        for  (i = 0;  i < n; i++) {
           v[i] = (1 << 15) - 1;   //Q0.15
           //printf("v[%d]=%hi ", i, v[i]);
        }
        //printf("\n");


        r = unsigned_int_fixed_vectors_dot_prod(v,v,n);  //result Q0.30
        r = ref_fxrsqrt(r>>14);                   // going from Qx.30 in dot product to Q16.16 
        for ( i = 0; i < n; i++) {
          v[i] = (long int)(v[i]*r) >> 16; // initial vector was normalized here Q0.15
           //printf("v[%d]=%hi ", i, v[i]);
        }
        printf("\n");
//        exit(0);
        lambda = 1;
//      Iterate till convergence:
        for (j = 0; j < maxit; j++) {

            // v must be normalized in advance:
            int_fixed_matrix_vector_mult(m, n, ATA, v, av, 15);      // here av --> Q3.12 * Q0.15 = Q3.27 --> av Q3.12  (av comps potentially > 1)
            // for ( i = 0; i < n; i++) printf("v[%d]=%hi ", i, v[i]); printf("\n");

            oldlambda = lambda; 
            lambda = int_fixed_vectors_dot_prod(v, av, n); // VT[ATA]V                      --> Q3.12 * Q0.15 = Q3.27 --

//          Check overflow:
            if (lambda < 0 ) {
               printf("Possible programming error. Overflow in lambda. Lambda=%li\nReduce scaling... stupid.\n", lambda) ;
               exit(3);
            }

            printf("Loop Insider itns=%3d Lambda-> %f components:", j, lambda/(float)(1<<27));
            for (i=0; i<n; i++) printf("%10.6f", v[i]/(32768.0-1)); printf("\n");

            r = unsigned_int_fixed_vectors_dot_prod(av,av,n);

//          This is what is available on our architecture, i.e. 1/sqrt(x):
            if ( r > 256 ) {
                r = ref_fxrsqrt(r>>8); // for vector renormalization  av * av --> Q3.12 * Q3.12 --> Q6.24--> 16.16 1/sqrtx
            } else {
                r = ref_fxrsqrt(r);
                r = r * 16;
                printf("*** Suspicious part\n");
            } 
            if ( r == 0 ) {
                for ( i = 0; i < n; i++) printf("av[%d]=%hi ", i, av[i]); printf("\n");
                printf(" 1/sqrt=%lu %f\n",r, r/65536.0) ;
                printf(" Possible programming error.\n") ;
                exit(2);
             }

            int diffmax = -10;
            for ( i = 0; i < n; i++) {
                v_old[i] = v[i];              // keep old vector v in memory               
                v[i] = (long int)(av[i] * r) >> 13;             //  got vector v is normalized now      Q3.12 * r --> Q3.12 * Q16.16 --> Q3.28 -> Q0.15
                diff = abs(v[i] - v_old[i]); // analyze differences 
                if ( diff > diffmax ) diffmax = diff; 

            } // new normalized vector v is ready
            //for (i=0; i<n; i++) printf("%10.6f", v[i]/(32768.0-1)); printf("\n");



            if ( ( diffmax ) <  1) {

                printf (" Desired accuracy has been achieved...\n");
                float sc = 1.0/(1<<27); // division is just for pretty printing
                printf (" Lambdas, old_lambda, diff: %g %g %10.5g\n", lambda*sc, oldlambda*sc, fabs(lambda*sc - oldlambda*sc)); 
                printf("Final result: after %d itns squared singular value is %g while singular vector components: (", j+1, lambda*sc);
                for (i=0; i<n; i++) {
//                   Continue printing (dont touch this):
                     printf("%10.6f", v[i]/(32768.0));
                     //V[i][k] = v[i]/32768.0;
                     V[i][k] = v[i];

//                   Since sqrt(x) is not available at ze moment we need some sophistication:
                     r = ref_fxrsqrt(lambda>>11);  // Q3.27 -> Q16.16
                     if ( r == 0 ) {
                         printf ("Possible programming error. r= %li lambda= %li\n", r, lambda);
                         exit(4);
                     }

                     // Another possibility instead of exit(4)??
                     //r = ref_fxrsqrt(lambda>>1);
                     //r = r << 5;

                     //r = 1.0f/sqrt(lambda); // this is what we have
                     S[k][k] = (long int)((lambda >> 15) * r) >> 16; // got sqrt(lambda) as a result Q3.27 * Q16.16 --> Q3.12  (43-31=12):
                }

                printf(" )\n\n");
//                exit(0);
                break;

            }
        } // itns loop

        if ( j > itns - 1) {
            printf("*** Error. No convergence after %d iterations.\n", j);
            exit(1);
        }

//      Mudify matrix by extracting eigenvector found:
        long int vv; long int ATApq;
        long int atamax=0;
        for (p=0; p< m; p++) {
            for( q = 0; q < n; q++) {
                vv = (long int)(v[p]*v[q]) >> 18;                            // Q0.15 * Q0.15 --> Q3.12;
                ATApq = ATA[p][q] - ((long int)((lambda >> 15) * vv) >> 12); // Q3.12 * Q3.12 --> Q6.24 ---> Q3.12; 
                if ( abs(ATApq) > atamax ) atamax = abs(ATApq);
            }
        }

        if ( atamax < 4096 && k < n - 1) {
           r= ref_fxrsqrt(65536*atamax);
           scale[k+1] = (unsigned long int)(r*r) >> 20; // Q16.16 * Q16.16 --> Q32.32 ---> Q3.12
           if ( scale[k+1] < 1 ) scale[k+1] = 1;
           printf("Rescaling ATA matrix... scale[%d]=%d\n", k+1, scale[k+1]);
        }   
        for (p = 0; p < m; p++) {
            for( q = 0; q < n; q++) {
                vv = (long int)(v[p]*v[q]) >> 18;                            // Q0.15 * Q0.15 --> Q3.12;
                ATApq = ATA[p][q] - ((long int)((lambda >> 15) * vv) >> 12); // Q3.12 * Q3.12 --> Q6.24 ---> Q3.12; 
                if ( abs(ATApq) < 4096 ) {
                    ATA[p][q] = ATA[p][q] - ((long int)((lambda >> 15) * vv) >> 12); // Q3.12 * Q3.12 --> Q6.24 ---> Q3.12; 
                    if ( k < n - 1) ATA[p][q]=ATA[p][q]*scale[k+1];
                    
                } else {
                  printf("Overflow in matrix ATA element: ATA[%d][%d]=%li/n", p, q, ATApq);
                  exit(2);
                }
            }
        }
        if ( k == n-1) printf("Digital rubbish...\n");
        int_fixed_matrix_fprintf(stdout, " %hi", m, n, ATA);
       // exit(8);
        

    } // dims loop
    scale[0] =1;
    for (k=1;k<n;k++) {
       scale[k] = scale[k]*scale[k-1];
       printf(" scale[%d]=%d", k, scale[k]);
    }
    printf("\n");

//  Lets do S matrix first it is still in Q3.12 format:
    for (k=0;k<n;k++) {
       r = ref_fxrsqrt(65536*scale[k]);
       //printf("r=1/sqrt(x)=%li x=%d\n", r, scale[k]);
       S[k][k] = (unsigned int long)(S[k][k]*r) >> 16;
       //printf("1/sqrt(x)=%f where x=%f\n", 1./sqrt(scale[k]),scale[k]);
       //r = scale(k)*ref_fxrsqrt((unsigned long int)scale(k)); 
    } 

//  reconstruct U matrix (works only for m==n or min(m,n):
    printf("Reconstructing U matrix from A and V\n");
    for (k=0; k<n; k++) {
        for(i=0; i<n; i++) v[i]=V[i][k];                       // Q0.15 fmt
        int_fixed_matrix_vector_mult(m, n, A, v, av,12);       // Q3.12 * Q0.15 -> Q3.27 --> Q0.15
        r = unsigned_int_fixed_vectors_dot_prod(av,av,n);
        r = ref_fxrsqrt(r>>14);
        for(i=0; i<n; i++) {
            U[i][k] = (long int)(av[i]*r) >> 16;           // Q0.15 fmt * Q16.16
        }
    }

//  Final result printing:
/*    printf("Matrix U is ready. Please compare.\n");
    matrix_fprintf(stdout, " %10.6f", m, n, UU);
    printf("Matrix S is ready. Please compare.\n");
    matrix_fprintf(stdout, " %10.6f", m, n, SS);
    printf("Matrix V is ready. Please compare.\n");
    matrix_fprintf(stdout, " %10.6f", m, n, VV);
*/
    
    return 0;
} // function end

