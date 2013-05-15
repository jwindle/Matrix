// Calculate LOWER Cholesky decomposition.

// The algorithm is:
// L11 = chol(A11)
// L21 L11' = A21
// L22 L22' = A22 - L21 L21'.
// Repeat.

#ifndef __CHOLESKY__

#include <stdexcept>
#include <vector>

#include "Matrix.h"

using std::vector;

////////////////////////////////////////////////////////////////////////////////
// COMPILER DEFINITIONS

// Column major indexing for sq matrix starting at 0 to array starting at 0.
#define CDX(i,j,ld)  ((i)+(j)*(ld))

// Column major indexing for sq matrix starting at 1 to array starting at 1.
#define FDX(i,j,ld)  ((i)+((j)-1)*(ld))

////////////////////////////////////////////////////////////////////////////////
// STRUCTS

// A structure for the calling parameters of syrk.
template<typename Real>
struct syrk_struct
{
  char  uplo;
  char  transa;
  int   n;
  int   k;
  Real  alpha;
  Real* a;
  int   lda;
  Real  beta;
  Real* c;
  int   ldc;

  syrk_struct(char uplo_, char transa_, int n_, int k_, Real alpha_, Real* a_, int lda_, Real beta_, Real* c_, int ldc_)
    : uplo(uplo_), transa(transa_), n(n_), k(k_), alpha(alpha_), a(a_), lda(lda_), beta(beta_), c(c_), ldc(ldc_) {}
};

template<typename Real>
struct trxx_struct
{
  char  side;
  char  uplo;
  char  transa;
  char  diag;
  int   m;
  int   n;
  Real  alpha;
  Real* a;
  int   lda;
  Real* b;
  int   ldb;
  
  trxx_struct(char side_, char transa_, char diag_, int m_, int n_, Real alpha_, Real* a_, int lda_, Real* b_, int ldb_)
    : side(side_), transa(transa_), diag(diag_), m(m_), n(n_), alpha(alpha_), a(a_), lda(lda_), b(b_), ldb(ldb_) {}
};

template<typename Real>
struct gemm_struct
{
  char  transa;
  char  transb;
  int   m;
  int   n;
  int   k;
  Real  alpha;
  Real* a;
  int   lda;
  Real* b;
  int   ldb;
  Real  beta;
  Real* c;
  int   ldc;

  gemm_struct(char transa_, char transb_, int m_, int n_, int k_, Real alpha_, Real* a_, int lda_, Real* b_, int ldb_, Real beta_, Real* c_, int ldc_)
    : transa(transa_), transb(transb_), m(m_), n(n_), k(k_), alpha(alpha_), a(a_), lda(lda_), b(b_), ldb(ldb_), beta(beta_), c(c_), ldc(ldc_) {}
};

////////////////////////////////////////////////////////////////////////////////
// CHOLESKY BY COLUMN

template<typename Real>
void chol_step_manual(Real* a, int lda, int n) {
  // a:   pointer to beginning of column major matrix.
  // lda: offset to the beginning of the next column.
  // n:   number of rows and columns of a
  // Updates LOWER triangular portion of a.

  #ifndef NDEBUG
  // if (a[0] <= 0) throw std::exception("chol_step_manual: a[0] <= 0");
  #endif

  // L11 = sqrt(A11)
  a[0] = sqrt(a[0]);

  // L21 * L11' = A21 
  for(int i=1; i<n; i++) a[i] /= a[0];

  // L22 = A22 - L21 * L21' IN LOWER TRIANGLE
  for(int j=1; j<n; j++)
    for(int i=j; i<n; i++)
      a[CDX(i,j,lda)] -= a[i]*a[j];

} // chol_step manual

template<typename Real>
void chol_step_blas(Real* a, int lda, int n)
{
  // a:   pointer to beginning of column major matrix.
  // lda: offset to the beginning of the next column.
  // n:   number of rows and columns of a
  // Updates LOWER triangular portion of a.

  #ifndef NDEBUG
  // if (a[0] <= 0) throw exception.
  #endif

  // L11 = sqrt(A11) --> A11* = sqrt(A11)
  a[0] = sqrt(a[0]); 

  if (n != 1) {

    // L21 * L11' = A21 
    rtrsm('R', 'L', 'T', 'N', n-1, 1, 1.0, a, lda, a+1, lda);
    // trxx_struct ts('R', 'L', 'N', 'N', n-1, 1, 1.0, a, lda, a+1, lda);
    // dtrsm(&ts.side, &ts.uplo, &ts.transa, &ts.diag, &ts.m, &ts.n, &ts.alpha, ts.a, &ts.lda, ts.b, &ts.ldb);

    // L22 L22' = A22 - L21 * L21' IN LOWER TRIANGLE.
    rsyrk('L', 'N', n-1, 1, -1.0, a+1, lda, 1.0, a+lda+1, lda);
    // syrk_struct ss('L', 'N', n-1, 1, -1.0, a+1, lda, a+lda+1, lda);
    // dsyrk_(&ss.uplo, &ss.trans, &ss.n, &ss.k, &ss.alpha, ss.a, &ss.lda, &ss.beta, ss.c, &ss.ldc);

  }

} // chol_step_blas

template<typename Real>
void chol(Real* a, int lda, int n)
{
  Real* a_curr = a;
  int   n_curr = n;
  for (int i=0; i<n; i++) {
    // chol_step_manual(a_curr, lda, n_curr);
    chol_step_blas(a_curr, lda, n_curr);
    a_curr += lda + 1;
    n_curr -= 1;
  }

  for(int j=1; j<n; j++)
    for(int i=0; i<j; i++)
      a[CDX(i,j,lda)] = 0.0;
} // chol

////////////////////////////////////////////////////////////////////////////////
// CHOLESKY BY BLOCK

template<typename Real>
void chol_block_manual(Real* a, int lda, int* part, int* bsize, int npart)
{
  // a:   pointer to beginning of column major matrix.
  // lda: offset to the beginning of the next column.
  // n: size of matrix
  // part: a partition of rows by where each block starts.
  // bsize: an array of block sizes corresponding to the size of each block. 
  // npart: number of elements part and bsize.

  // Updates the BLOCK lower triangular portio of a.  If a is not a symmetric
  // matrix to start then the resulting diagonal blocks will not look symmetric.

  // We could also work off of a total offset from a.  So a wouldn't change.

  // L11 = sqrt(A11)
  chol(a + CDX(part[0], part[0], lda), lda, bsize[0]);
  Real* head = a + CDX(part[0], part[0], lda);

  // L21 * L11' = A21 
  for(int i=1; i<npart; i++) {
    Real* Ai0 = a + CDX(part[i], part[0], lda); // Maybe change this.
    // Real* Ai0 = a + i; // Maybe change this.
    rtrsm('R', 'L', 'T', 'N', bsize[i], bsize[0], 1.0, a, lda, Ai0, lda);
  }

  // L22 L22' = A22 - L21 * L21' IN LOWER TRIANGLE OF BLOCKS
  for(int j=1; j<npart; j++) {
    Real* Aj0 = a + CDX(part[j], part[0], lda);
    for(int i=j; i<npart; i++) {
      Real* Ai0 = a + CDX(part[i], part[0], lda);
      Real* Aij = a + CDX(part[i], part[j], lda);
      // long int diff = Aij - head;
      // printf("diff: %li ", diff);
      printf("Ai0, Aj0, Aij: %li, %li, %li\n", Ai0-head, Aj0-head, Aij-head);
      rgemm('N', 'T', bsize[i], bsize[0], bsize[j], -1.0, Ai0, lda, Aj0, lda, 1.0, Aij, lda);
      // gemm_struct in('N', 'T', bsize[i], bsize[0], bsize[j], -1.0, Ai0, lda, Aj0, lda, 1.0, Aij, lda);
      // dgemm_(&in.transa, &in.transb, &in.m, &in.n, &in.k, &in.alpha, in.a, &in.lda, in.b, &in.ldb, &in.beta, in.c, &in.ldc);
      // Could do a rsyk if i=j.
    }
    printf("\n");
  }


} // chol_block_manual

template<typename Real>
void chol_block_blas(Real* a, int lda, int n, int bsize)
{
  // a:   pointer to beginning of column major matrix.
  // lda: offset to the beginning of the next column.
  // n:   number of rows and columns of a
  // blocks: number of blocks
  // bsize:  number of rows and columns in each block
  // Updates the lower triangular portio of a.

  // L11 = sqrt(A11)
  chol(a, lda, bsize);

  if (n != bsize) {

    // L21 * L11' = A21 
    rtrsm('R', 'L', 'T', 'N', n - bsize, bsize, 1.0, a, lda, a + bsize, lda);
    
    // L22 L22' = A22 - L21 * L21' IN LOWER TRIANGLE
    rsyrk('L', 'N', n - bsize, bsize, -1.0, a+bsize, lda, 1.0, a+CDX(bsize, bsize, lda), lda);

  }

} // chol_block_blas

template<typename Real>
void chol(Real* a, int lda, int n, int* bsize, int npart)
{
  Real* a_curr = a;
  int   n_curr = n;

  vector<int> part(npart); part[0] = 0;
  for(int i=1; i<npart; i++) part[i] = part[i-1] + bsize[i-1];

  for (int i=0; i<npart; i++) {
    // printf("i, part[i], n_curr, a_curr: %i, %i, %i, %p\n", i, part[i], n_curr, a_curr);
    chol_block_manual(a, lda, &part[i], bsize+i, npart-i);

    // chol_block_blas(a_curr, lda, n_curr, bsize[i]);

    // chol_block_manual(a_curr, lda, &part[i], bsize+i, npart-i);
    // int remove = part[i+1];
    // for (int j=i+1; j<npart; j++) {
    //   part[j] -= remove; // Could be part of a CUDA step.
    // }

    a_curr += CDX(bsize[i], bsize[i], lda);
    n_curr -= bsize[i];

  }
} // chol

template<typename Real>
void chol(Real* a, int lda, int n, int bsize)
{
  int remain = n % bsize;
  int npart  = n / bsize + (remain ? 1 : 0);

  printf("n, bsize, remain, npart: %i, %i, %i, %i\n", n, bsize, remain, npart);

  vector<int> size(npart);

  for (int i=0; i<npart-1; i++)
    size[i] = bsize;
  size[npart-1] = remain ? remain : bsize;

  chol(a, lda, n, &size[0], npart);

  for(int j=1; j<n; j++)
    for(int i=0; i<j; i++)
      a[CDX(i,j,lda)] = 0.0;

} // chol

////////////////////////////////////////////////////////////////////////////////

#undef CDX
#undef FDX

#endif
