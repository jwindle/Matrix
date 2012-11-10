#include "Matrix.h"

//////////////////////////////////////////////////////////////////////
		   // SIMPLE CONJUGATE GRADIENTS //
//////////////////////////////////////////////////////////////////////

// http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// by Jonathan Richard Shewchuk

// Assumes A is positive definite.

int cg(MF& x, const MF& A, const MF& b, double tol, int max_iter)
{
  uint P = b.rows();

  int iter = 0;

  Matrix r(b);
  gemm(r, A, x, 'N', 'N', -1.0, 1.0);

  Matrix d(r);

  double delta_new = dot(r, r);
  // See note below...
  // double delta_0   = delta_new;
  double delta_old = 0;

  Matrix q(P);
  double alpha;
  double beta;

  //while(iter < max_iter && delta_new > tol*delta_0){
  //The above condition is suggested by Shewchuk, but we want an absolute tolerance.
  while(iter < max_iter && delta_new > tol){
    gemm(q, A, d);
    alpha = delta_new / dot(d, q);
    axpy(alpha, d, x);

    // You might want to discard this step.
    if (iter % 50 == 0) {
      r.clone(b);
      gemm(r, A, x, 'N', 'N', -1.0, 1.0);
    }
    else{
      axpy(-1.0 * alpha, q, r);
    }

    delta_old = delta_new;
    delta_new = dot(r, r);

    beta = delta_new / delta_old;
    hsum(d, d, r, beta, 0.0);
    //for(uint j = 0; j < P; j++)
    //  d(j) = beta * d(j) + r(j);

    iter++;
  }

  return iter;
}

//////////////////////////////////////////////////////////////////////
		     // BLAS / LAPACK WRAPPERS //
//////////////////////////////////////////////////////////////////////

// Solve a symmetric, positive definite system of equations.
// ax = b -> b := x;
// int symsolve(MF a, Matrix& b, char uplo='L')
int symsolve(MF a, Matrix& b, char uplo)
{
  // b.fill(0.0);
  // for(uint i = 0; i < b.cols(); ++i)
  //   b(i,i) = 1.0;

  Matrix temp(a);
  int info = posv(temp, b, uplo);

  if (info) {
    fprintf(stderr, "Problem with symsolve: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);

    throw std::runtime_error("potrf failed\n");
  }

  return info;
}

//------------------------------------------------------------------------------
// int syminv(MF a, Matrix& ainv, char uplo='L')
int syminv(MF a, Matrix& ainv, char uplo)
{
  ainv.resize(a.rows(), a.cols());
  ainv.fill(0.0);
  for(uint i = 0; i < ainv.cols(); ++i)
     ainv(i,i) = 1.0;

  return symsolve(a, ainv, uplo);
}

//------------------------------------------------------------------------------
// Get the Cholesky decomposition of a matrix a.
// int chol(Matrix& c, MF a, char uplo='L')
int chol(Matrix& c, MF a, char uplo)
{
  c.clone(a);
  int info = chol(c, uplo);

  // FIX FIX Set other entries to zero.
  // Do I want to have an option for letting a person not do this?

  if (uplo=='L') {
    for (uint j = 1; j < c.cols(); ++j)
      for (uint i = 0; i < j; ++i)
	c(i,j) = 0.0;
  }
  else {
    for (uint j = 0; j < c.cols(); ++j)
      for (uint i = j+1; i < c.cols(); ++i)
	c(i,j) = 0.0;
  }

  if (info) {
    fprintf(stderr, "Problem with chol: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);

    throw std::runtime_error("potrf failed\n");
  }

  return info;
}

//------------------------------------------------------------------------------
void dsyevd(char jobz, char uplo, int n, double* a, int lda, double* w, double* work, int lwork, int* iwork, int liwork, int* info)
{ return dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }

int symeigen(Matrix& evec, Matrix& eval, MF symmat)
{
  sizecheck(symmat.rows()==symmat.cols());
  int info = 0;
  int N = symmat.rows();
  
  evec.clone(symmat);
  eval.resize(N);

  // int lwork  = 1 + 6 * N + 2*N*N;
  // int liwork = 3 + 5*N;

  int lwork = -1;
  int liwork = -1;

  std::vector<double> work(1);
  std::vector<int>    iwork(1);

  dsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);

  lwork = (int) work[0];
  liwork = (int) iwork[0];

  work.resize(lwork);
  iwork.resize(liwork);

  dsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);
  
  if (info != 0) {
    fprintf(stderr, "problem in symeigen; info=%i.\n", info);
    throw std::runtime_error("symeigen failed\n");
  }

  return info;
}

//------------------------------------------------------------------------------
int symsqrt(Matrix& rt, MF symmat)
{
  int N = symmat.rows();

  Matrix evec;
  Matrix eval;
  int info = symeigen(evec, eval, symmat);

  Matrix d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

//------------------------------------------------------------------------------
int syminvsqrt(Matrix& rt, MF symmat)
{
  int N = symmat.rows();

  Matrix evec;
  Matrix eval;
  int info = symeigen(evec, eval, symmat);

  Matrix d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(1.0 / eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

//--------------------------------------------------------------------

void dgesvd(char jobu, char jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double*work, int lwork, int* info)
{ return dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }

int svd(Matrix& U, Matrix& S, Matrix& tV, Matrix& X)
{
  // Always calculate a thinned U and a non-thinned V.  Pad out the
  // singular values by 0.0 when needed.
  Matrix A(X);

  int m = A.rows();
  int n = A.cols();

  tV.resize(n,n);
  S.resize(n); S.fill(0.0); // The singular values.
  // Should be of size min(m,n), but I want to pad things out.
  
  if (m >= n) U.resize(m,n);
  else U.resize(m,m);
  int ldu = U.rows();

  char jobu  = 'S';
  char jobvt = 'A';

  int maxmn = m > n ? m : n;
  int minmn = m < n ? m : n;

  vector<double> work(1);
  int lwork = -1;

  int info;

  // Workspace query.
  dgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);
  printf("lwork: %g\n", work[0]);
  
  // SVD.
  lwork = (int)work[0];
  lwork = lwork > 1 ? lwork : 5 * minmn + maxmn;

  work.resize(lwork);
  dgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);

  if (info != 0) {
    fprintf(stderr, "problem in svd; info=%i.\n", info);
    throw std::runtime_error("svd failed\n");
  }

  return info;
}

//------------------------------------------------------------------------------
// this = op(a) * op(b) + alpha * this.
void mult(Matrix& c, const Frame<double>& a, const Frame<double>& b, char ta, char tb, double alpha, double beta)
{
  uint opa_rows = ta=='T' ? a.cols() : a.rows();
  uint opb_cols = tb=='T' ? b.rows() : b.cols();
  c.resize(opa_rows, opb_cols, 1);
  
  gemm(c, a, b, ta, tb, alpha, beta);
}
