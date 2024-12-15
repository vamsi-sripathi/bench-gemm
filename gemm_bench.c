#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#if defined SINGLE
#define FP_TYPE    float
#define FNAME      sgemm
#define COMPSIZE   1
#elif defined DOUBLE
#define FP_TYPE    double
#define FNAME      dgemm
#define COMPSIZE   1
#elif defined COMPLEX8
#define FP_TYPE    MKL_Complex8
#define FNAME      cgemm
#define COMPSIZE   4
#elif defined COMPLEX16
#define FP_TYPE    MKL_Complex16
#define FNAME      zgemm
#define COMPSIZE   4
#else
#error FP_TYPE is undefined
#endif

#define STRINGIFY(x)  STRINGIFY_(x)
#define STRINGIFY_(x) #x

#define  NUM_TRIALS       (30)
#define  PAGE_ALIGNMENT   (4096)
#define  AVX512_ALIGNMENT (64)

/* #define  VERBOSE */
/* #define  PAD_LD */
/* #define  FORCE_UNALIGN */

static inline int fix_ld(int ld)
{
  int bad_ld_alignment = 256*sizeof(FP_TYPE);
  int cache_line_alignment = 64;
  int ld_bytes = ld*sizeof(FP_TYPE);
  int padded_ld;

#if 0
  return ((ld_bytes + bad_ld_alignment - 1)/bad_ld_alignment*bad_ld_alignment + cache_line_alignment)/sizeof(FP_TYPE);
#endif

  if (ld_bytes%bad_ld_alignment) {
    padded_ld = ((ld_bytes + cache_line_alignment - 1)/cache_line_alignment)*cache_line_alignment;
    if (padded_ld%bad_ld_alignment) {
      return padded_ld/sizeof(FP_TYPE);
    } else {
      return (padded_ld+cache_line_alignment)/sizeof(FP_TYPE);
    }
  } else {
    return (ld_bytes + cache_line_alignment)/sizeof(FP_TYPE);
  }
}

int main (int argc, char **argv)
{
  char       transa, transb;
  int        m, n, k, lda, ldb, ldc;
  int        i, j, t, ax, bx;
  int        a_rows, a_cols, b_rows, b_cols, c_rows, c_cols;
  double     t_start, t_stop, curr_gemm_gflops, best_gemm_gflops;
  FP_TYPE    alpha, beta;
  FP_TYPE    *a = NULL, *b = NULL, *c = NULL;
  int offset = 0;

  if (argc != 8) {
    printf ("\n USAGE: %s transa<n|t> transb<n|t> m<int> n<int> k<int> alpha<float> beta<float>\n", argv[0]); fflush(0);
    exit (1);
  }

  transa      = argv[1][0];
  transb      = argv[2][0];
  m           = atoi(argv[3]);
  n           = atoi(argv[4]);
  k           = atoi(argv[5]);
#if defined (SINGLE) || defined (DOUBLE)
  alpha       = atof(argv[6]);
  beta        = atof(argv[7]);
#else
  alpha.real  = atof(argv[6]);
  alpha.imag  = 0.;
  beta.real   = atof(argv[7]);
  beta.imag   = 0.;
#endif

  ax     = ((transa == 'N' || transa == 'n') ? k : m);
  bx     = ((transb == 'N' || transb == 'n') ? n : k);
  lda    = ((transa == 'N' || transa == 'n') ? m : k);
  ldb    = ((transb == 'N' || transb == 'n') ? k : n);
  ldc    = m;
  a_rows = lda;
  a_cols = ax;
  b_rows = ldb;
  b_cols = bx;
  c_rows = ldc;
  c_cols = n;

#ifdef PAD_LD
  lda = fix_ld(lda);
  ldb = fix_ld(ldb);
  ldc = fix_ld(ldc);
#endif

#ifdef FORCE_UNALIGN
  offset++;
#endif

  FP_TYPE *a_mem = (FP_TYPE *) mkl_malloc (((1LL*lda*ax)+offset)*sizeof(FP_TYPE), AVX512_ALIGNMENT);
  FP_TYPE *b_mem = (FP_TYPE *) mkl_malloc (((1LL*ldb*bx)+offset)*sizeof(FP_TYPE), AVX512_ALIGNMENT);
  FP_TYPE *c_mem = (FP_TYPE *) mkl_malloc (((1LL*ldc*n)+offset)*sizeof(FP_TYPE),  AVX512_ALIGNMENT);

  if (a_mem == NULL || b_mem == NULL || c_mem == NULL) {
    printf ("Memory allocation failed!\n");
    exit(1);
  }

  a = a_mem+offset;
  b = b_mem+offset;
  c = c_mem+offset;

#ifdef FORCE_UNALIGN
  printf ("Alignment offset of A from %s byte boundary = %ld\n",
          STRINGIFY(AVX512_ALIGNMENT), (unsigned long)a%AVX512_ALIGNMENT);fflush(0);
  printf ("Alignment offset of B from %s byte boundary = %ld\n",
          STRINGIFY(AVX512_ALIGNMENT), (unsigned long)b%AVX512_ALIGNMENT);fflush(0);
  printf ("Alignment offset of C from %s byte boundary = %ld\n",
          STRINGIFY(AVX512_ALIGNMENT), (unsigned long)c%AVX512_ALIGNMENT);fflush(0);
#endif

  for (j=0; j<a_cols; j++) {
    for (i=0; i<a_rows; i++) {
#if defined (SINGLE) || defined (DOUBLE)
      a[(j*lda)+i] = rand()/((double) RAND_MAX - 0.5);
      /* a[(j*lda)+i] = (i>=a_rows) ? NAN : rand()/((double) RAND_MAX - 0.5); */
#else
      a[(j*lda)+i].real = a[(j*lda)+i].imag = rand()/((double) RAND_MAX - 0.5);
#endif
    }
  }

  for (j=0; j<b_cols; j++) {
    for (i=0; i<b_rows; i++) {
#if defined (SINGLE) || defined (DOUBLE)
      b[(j*ldb)+i] = rand()/((double) RAND_MAX - 0.5);
#else
      b[(j*ldb)+i].real = b[(j*ldb)+i].imag = rand()/((double) RAND_MAX - 0.5);
#endif
    }
  }

  for (j=0; j<c_cols; j++) {
    for (i=0; i<c_rows; i++) {
#if defined (SINGLE) || defined (DOUBLE)
      c[(j*ldc)+i] = 0.;
#else
      c[(j*ldc)+i].real = c[(j*ldc)+i].imag = 0.;
#endif
    }
  }

  // Warm-up call
  t_start = dsecnd();
  t_start = dsecnd();
  FNAME (&transa, &transb, &m ,&n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);

  long long n_ops = 2LL * COMPSIZE * m * n * k;
  double avg_gemm_gflops = 0.0;

  for (t=0; t<NUM_TRIALS; t++) {
    t_start = dsecnd();
    FNAME (&transa, &transb, &m ,&n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    t_stop = dsecnd();

    curr_gemm_gflops = n_ops/(t_stop-t_start) * 1.e-9;
    avg_gemm_gflops += curr_gemm_gflops;

    if (t) {
      if (curr_gemm_gflops > best_gemm_gflops) {
        best_gemm_gflops = curr_gemm_gflops;
      }
    } else {
      best_gemm_gflops = curr_gemm_gflops;
    }
#ifdef VERBOSE
    printf ("%s: transa= %c, transb= %c, m= %d, n= %d, k= %d, lda= %d, ldb= %d, ldc= %d, gflops= %.2f\n",
            STRINGIFY(FNAME), transa, transb, m, n, k, lda, ldb, ldc,
            curr_gemm_gflops); fflush(0);
#endif
  }
  printf ("%s: transa= %c, transb= %c, m= %d, n= %d, k= %d, lda= %d, ldb= %d, ldc= %d, avg-gflops= %.2f, best-gflops= %.2f\n",
          STRINGIFY(FNAME), transa, transb, m, n, k, lda, ldb, ldc,
          avg_gemm_gflops/NUM_TRIALS, best_gemm_gflops); fflush(0);

  mkl_free (a_mem);
  mkl_free (b_mem);
  mkl_free (c_mem);

  return 0;
}
