// -*- mode: c++; -*-

#include "Matrix.h"
#include <iostream>
#include <sstream>
#include <string>

using std::cout;
using std::string;
using std::stringstream;
using std::vector;

//using namespace Mine;

// typedef float TReal;
typedef double TReal;

//------------------------------------------------------------------------------
void test_svd2(Block<TReal>& A) {

  // SVD
  // A << " 4 3 0 "
  //      " 0 4 1 "
  //      " 0 1 4 ";

  Block<TReal> U;
  Block<TReal> S;
  Block<TReal> tV;

  svd2(U, S, tV, A);
  cout << "A:\n" << A;
  cout << "U:\n" << U;
  cout << "S:\n" << S;
  cout << "tV:\n" << tV;

}

//------------------------------------------------------------------------------
void test_svd(Block<TReal>& A, char jobz='A', bool padS = false) {

  // SVD
  // A << " 4 3 0 "
  //      " 0 4 1 "
  //      " 0 1 4 ";

  Block<TReal> U;
  Block<TReal> S;
  Block<TReal> tV;

  svd(U, S, tV, A, jobz, padS);
  cout << "A:\n" << A;
  cout << "U:\n" << U;
  cout << "S:\n" << S;
  cout << "tV:\n" << tV;

}

//------------------------------------------------------------------------------
void test_cg(Block<TReal>& A)
{

  // // CONJUGATE GRADIENT
  Block<TReal> b("N", 3);

  Block<TReal> x(b);
  symsolve(A, x);

  cout << "x:\n" << x;

  Block<TReal> y(b);
  try {
  cout << "iter: " << cg(y, A, b, (TReal)10e-8, 100) << "\n";
  }
  catch (std::exception& e) {
    cout << e.what();
  }

  cout << "y:\n" << y;

}

//------------------------------------------------------------------------------
void test_block_constructor(){

  Block<int> C("W", 3);
  cout << "C:\n" << C << "\n";

}

//------------------------------------------------------------------------------
void test_read_write() {

  Matrix A("I", 3);
  cout << A.write("A.mat", true) << "\n";

  A.resize(1);
  A.read("A.mat", true);
  cout << A;

}

//------------------------------------------------------------------------------
void test_readNatural() {
  Matrix A;

  A.readNatural("Nat1.mat");
  A.out(cout, true) << "\n";

  A.readNatural("Nat2.mat");
  A.out(cout, true) << "\n";
}

//------------------------------------------------------------------------------
void test_out(Block<TReal>& A) {
  A.out(cout, true);

  cout << "\n";

  A.out(cout, false);

  cout << "\n";
}

//------------------------------------------------------------------------------
void test_syrk() {

  Block<TReal> A(3,3);

  A << " 4 0 0 "
       " 0 4 1 "
       " 0 1 4 ";

  cout << "A:\n" << A << "\n";

  Block<TReal> AA(3, 3); AA.fill(1.0);
  syrk(AA, A, 'N', (TReal)1.0, (TReal)0.0);

  cout << "AA:\n" << AA << "\n";

  gemm(AA, A, A, 'N', 'T');
  cout << "AA:\n" << AA << "\n";

}

//------------------------------------------------------------------------------
void test_symeigen(Block<TReal>& A) {

  // EIGENVECTOR DECOMP OF SYM MATRIX
  Block<TReal> evec;
  Block<TReal> eval;
  Block<TReal> rt;
  Block<TReal> irt;

  symeigen(evec, eval, A);
  symsqrt(rt, A);
  syminvsqrt(irt, A);

  cout << "A:\n" << A;
  cout << "U:\n" << evec;
  cout << "D:\n" << eval;

  cout << "rt:\n" << rt;
  cout << "irt:\n" << irt;

}

//------------------------------------------------------------------------------
void test_scan() {
  // SCAN
  Matrix b(2,2);
  b.scan("b.dat");
  cout << "b:\n" << b;
}

//------------------------------------------------------------------------------
void test_max_and_min(Block<TReal>& A) {

  cout << "Max: " << maxAll(A) << "\n";
  cout << "Min: " << minAll(A) << "\n";

}

//------------------------------------------------------------------------------
void test_op() {

  Block<TReal> A(3, 3);

  A << " 1 2 3 "
       " 4 5 6 "
       " 7 8 9 " ;

  cout << A;

  Block<TReal> Id("I", 3);
  Block<TReal> B = 5.0 * Id;

  Block<TReal> C(0.0, B);

  cout << B ;

  cout << hmin(B, A);

  Frame<TReal> AMF(&A(0), 3, 1);
  B = Block<TReal>("I", 2);
  B = (A + A) + (TReal)1.0 + AMF;
  cout << B;

}

//------------------------------------------------------------------------------
void test_copy() {

  Block<TReal> A(3, 3);

  A << " 1 2 3 "
       " 4 5 6 "
       " 7 8 9 " ;

  Block<TReal> rows(2); rows << "0 1";
  Block<TReal> cols(2); cols << "1 2";
  Block<TReal> B(2,2);
  B.copy(A, rows, cols);
  cout << "A:\n" << A << "\n";
  cout << "B:\n" << B << "\n";

}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

  // Test SVD;
  Block<TReal> A(3,4);

  // Read in transposed.
  A << " 4 3 0 "
       " 0 4 1 "
       " 0 1 4 "
       " 3 1 1 ";

  // test_svd(A, 'A');
  // test_svd(A, 'S');
  // test_svd2(A);

  // A.resize(3,4);

  // test_svd(A, 'A', true);
  // test_svd(A, 'S', true);
  // test_svd2(A);

  test_readNatural();

  // test_out(A);

  return 0;

}
