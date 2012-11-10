/*
  This file contains the functions that make use of MatrixFrame.
  These include Blas wrappers for matrix arithmetic, LAPACK
  wrappers, and variouis matrix manipulations.

  You probably want to view this file with at least a 160 character
  wide screen.

  The BLAStoC.perl file will help you write BLAS wrappers.
 */

#include "MatrixFrame.h"

//////////////////////////////////////////////////////////////////////
			 // BLAS / LAPACK //
//////////////////////////////////////////////////////////////////////

/*
  See BLAS / LAPACK documentation at netlib.org.  The Fortran source
  code is very regular.  The perl script BLAStoC.perl will take a BLAS
  subroutine/function file and convert it to the necessary C code.
 */

//////////////////////////////////////////////////////////////////////
		    // MatrixFrame BLAS WRAPPER //
//////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------//
// y = alpha x + y.

void daxpy(int n, double da, double* dx, int incx, double* dy, int incy)
{ daxpy_(&n, &da, dx, &incx, dy, &incy); }

void axpy(double alpha, MF x, MF y)
{
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  daxpy((int)x.rows(), alpha, &x(0), 1, &y(0), 1);
}

//------------------------------------------------------------------//
// x'y

double ddot(int n, double* dx, int incx, double* dy, int incy)
{ return ddot_(&n, dx, &incx, dy, &incy); }

double dot(MF x, MF y)
{
  #ifndef NDEBUG
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  #endif
  return ddot(x.rows(), x.getp(), 1, y.getp(), 1);
}

//------------------------------------------------------------------//
// c = alpha op(a) * op(b) + beta c.

void dgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{ dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

// void gemm(MF c, MF a, MF b, char ta='N', char tb='N', double alpha=1.0, double beta=0.0)
void gemm(MF c, MF a, MF b, char ta, char tb, double alpha, double beta)
{
  #ifndef NDEBUG
  memcheck(!overlap(c,a) && !overlap(c,b));
  #endif
  // Get the dimensionality information we need.
  int cnr = (int)c.rows(); int cnc = (int)c.cols();
  int anr = (int)a.rows(); int bnr = (int)b.rows();
  int k   = (int)pconform(c, a, b, ta, tb);
  // Make sure things conform.
  #ifndef NDEBUG
  sizecheck(k!=0);
  #endif
  dgemm(ta, tb, cnr, cnc, k, alpha, &a(0), anr, &b(0), bnr, beta, &c(0), cnr);
} // gemm

//------------------------------------------------------------------//
// b = alpha op(a) * b  OR  b = alpha b * op(a) where a is triangular.

void dtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

// void trmm(MF a, MF b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
void trmm(MF a, MF b, char uplo, char side, char ta, char diag, double alpha)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='L' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  dtrmm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trmm

//------------------------------------------------------------------//
// Solve x:  op(a) x = alpha b  OR  x op(a) = alpha b, a triangular.
// i.e: x = alpha inv(op(a)) b  OR  x = alpha b inv(op(a)).
// The solution is overwriten into b.

void dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

// void trsm(MF a, MF b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
void trsm(MF a, MF b, char uplo, char side, char ta, char diag, double alpha)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='L' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  dtrsm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trsm

//------------------------------------------------------------------------------
// C := alpha*A*A**T + beta*C, ta='N'
// C := alpha*A**T*A + beta*C. ta='T'

void dsyrk(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double beta, double* c, int ldc)
{ return dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }

void syrk(MF c, MF a, char ta, double alpha, double beta)
{
  memcheck(!overlap(c,a));
  char tb = ta=='N' ? 'T' : 'N';
  pconform(c, a, a, ta, tb);
  int k = ta=='N' ? a.cols() : a.rows();
  
  dsyrk('U', ta, c.rows(), k, alpha, &a(0), a.rows(), beta, &c(0), c.rows());

  // Better way?
  for(int j=0; j<c.cols()-1; j++)
    for(int i=j+1; i<c.rows(); i++)
      c(i,j) = c(j,i);
}

//////////////////////////////////////////////////////////////////////
		  // MATRIX FRAME LAPACK WRAPPER //
//////////////////////////////////////////////////////////////////////

// Solve a general linear system, ax = b for x.

void dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int& info)
{ dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); }

int gesv(MF a, MF b)
{
  memcheck(!overlap(a, b));       // No overlap in memory.
  sizecheck(pconform(b, a, b)!=0); // a is square and b conforms.
  int info;
  std::vector<int> ipiv(a.rows());
  dgesv(a.rows(), b.cols(), &a(0), a.rows(), &ipiv[0], &b(0), b.rows(), info);
  return info;
}

// Shorthand.
int solve(MF a, MF b)
{
  return gesv(a, b);
}

//------------------------------------------------------------------//
// Solves ax = b for x where a is sym. pos. def.  Note: the lower (or
// upper) portion of A is overwritten with the Cholesky decomposition.

void dposv(char uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info)
{ dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info); }

int posv(MF a, MF b, char uplo)
{
  memcheck(!overlap(a,b));
  sizecheck(pconform(b, a, b)!=0);
  int info;
  dposv(uplo, a.rows(), b.cols(), &a(0), a.rows(), &b(0), b.rows(), info);

  if (info != 0) {
    printf("Error in posv: info = %i\n", info);
    throw std::runtime_error("aborted in posv\n");
  }

  return info;
}

//------------------------------------------------------------------//
// Cholesky Decomposition

void dpotrf(char uplo, int n, double* a, int lda, int& info)
{ dpotrf_(&uplo, &n, a, &lda, &info); }

int potrf(MF a, char uplo)
{
  sizecheck(a.rows()==a.cols());
  int info = 0;
  dpotrf(uplo, a.rows(), &a(0), a.rows(), info);
  return info;
}

// int chol(MF a, char uplo='L')
int chol(MF a, char uplo)
{
  return potrf(a, uplo);
}

//------------------------------------------------------------------//

