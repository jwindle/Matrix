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
  cout << "A:\n" << A << "\n";
  cout << "U:\n" << U << "\n";
  cout << "S:\n" << S << "\n";
  cout << "tV:\n" << tV << "\n";

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
  cout << "A:\n" << A << "\n";
  cout << "U:\n" << U << "\n";
  cout << "S:\n" << S << "\n"; cout << "tV:\n" << tV << "\n";

}

//------------------------------------------------------------------------------
void test_cg(Block<TReal>& A)
{

  // // CONJUGATE GRADIENT
  Block<TReal> b("N", 3);

  Block<TReal> x(b);
  symsolve(A, x);

  cout << "x:\n" << x << "\n";

  Block<TReal> y(b);
  try {
    cout << "iter: " << cg(y, A, b, (TReal)10e-8, 100) << "\n";
  }
  catch (std::exception& e) {
    cout << e.what();
  }

  cout << "y:\n" << y << "\n";

}

//------------------------------------------------------------------------------
void test_block_constructor(){

  printf("Test block constructor:\n");

  Block<int> C("W", 3);
  cout << "C:\n" << C << "\n";

}

//------------------------------------------------------------------------------
void test_dump_load() {

  printf("Test dump/load:\n");

  Matrix A("I", 3);
  cout << A.dump("A.mat", true) << "\n";

  A.resize(1);
  A.load("A.mat", true);
  cout << A << "\n";

}

//------------------------------------------------------------------------------
void test_read() {
  Matrix A;

  printf("Test read:\n");

  printf("Read natural.\n");

  A.read("Nat1.mat");
  A.out(cout, true) << "\n";

  A.read("Nat2.mat");
  A.out(cout, true) << "\n";

  printf("Read transpose.\n");

  A.read("Nat1.mat", false);
  A.out(cout, true) << "\n";

  A.read("Nat2.mat", false);
  A.out(cout, true) << "\n";
}

//------------------------------------------------------------------------------
void test_out(Block<TReal>& A) {
  printf("Test out:\n");

  A.out(cout, true) << "\n";

  A.out(cout, false) << "\n";
}

//------------------------------------------------------------------------------
void test_syrk() {

  printf("Test syrk:\n");

  Block<TReal> A(3,3);

  A.scanString( " 4 0 0 "
                " 0 4 1 "
                " 0 1 4 " );

  cout << "A:\n" << A << "\n";

  Block<TReal> AA(3, 3); AA.fill(1.0);
  syrk(AA, A, 'N', (TReal)1.0, (TReal)0.0);

  cout << "AA:\n" << AA << "\n";

  gemm(AA, A, A, 'N', 'T');
  cout << "AA:\n" << AA << "\n";

}

//------------------------------------------------------------------------------
void test_symeigen(Block<TReal>& A) {

  printf("Test symeigen:\n");

  // EIGENVECTOR DECOMP OF SYM MATRIX
  Block<TReal> evec;
  Block<TReal> eval;
  Block<TReal> rt;
  Block<TReal> irt;

  symeigen(evec, eval, A);
  symsqrt(rt, A);
  syminvsqrt(irt, A);

  cout << "A:\n" << A << "\n";
  cout << "U:\n" << evec << "\n";
  cout << "D:\n" << eval << "\n";

  cout << "rt:\n" << rt << "\n";
  cout << "irt:\n" << irt << "\n";

}

//------------------------------------------------------------------------------
void test_scan() {

  printf("Test scan:\n");

  Matrix b(2,2);
  b.scan("b.dat");
  cout << "b:\n" << b;
}

//------------------------------------------------------------------------------
void test_max_and_min(Block<TReal>& A) {

  printf("Test max and min:\n");

  cout << "Max: " << maxAll(A) << "\n";
  cout << "Min: " << minAll(A) << "\n";

}

//------------------------------------------------------------------------------
void test_op() {

  printf("Test op:\n");

  Block<TReal> A(3, 3);

  A.scanString( " 1 2 3 "
                " 4 5 6 "
                " 7 8 9 " );

  cout << A << "\n";

  Block<TReal> Id("I", 3);
  Block<TReal> B = 5.0 * Id;

  Block<TReal> C(0.0, B);

  cout << B << "\n";

  cout << hmin(B, A);

  Frame<TReal> AMF(&A(0), 3, 1);
  B = Block<TReal>("I", 2);
  B = (A + A) + (TReal)1.0 + AMF;
  cout << B;

}

//------------------------------------------------------------------------------
void test_copy() {

  printf("Test copy:\n");

  Block<TReal> A(3, 3);

  A.scanString( " 1 2 3 "
                " 4 5 6 "
                " 7 8 9 " );

  Block<TReal> rows(2); rows.scanString("0 1");
  Block<TReal> cols(2); cols.scanString("1 2");
  Block<TReal> B(2,2);
  B.copy(A, rows, cols);
  cout << "A:\n" << A << "\n";
  cout << "B:\n" << B << "\n";

  Block<TReal> C(A);
  cout << "A:\n" << A << "\n";
  cout << "C:\n" << C << "\n";

  try {
    Frame<TReal> D(A);
    D.copy(A);
  }
  catch(std::exception& e) {
    std::cerr << e.what();
  }
}

//------------------------------------------------------------------------------
void test_reshape() {

  printf("Test reshape:\n");

  Block<TReal> A(3, 3);

  A.scanString( " 1 2 3 "
                " 4 5 6 "
                " 7 8 9 " );

  A.reshape(1,9);

  cout << A << "\n";

}

//------------------------------------------------------------------------------
void test_cbind() {

  printf("Test cbind:\n");

  Block<TReal> A(3, 3);

  A.scanString( " 1 2 3 "
                " 4 5 6 "
                " 7 8 9 " );

  Block<TReal> B(3,1);
  B.scanString( " 4 5 6 " );

  A.cbind(A);
  cout << A << "\n";

  A.cbind(B);
  cout << A << "\n";

  Block<TReal> C(2,1);
  C.scanString( " 4 5 " );

  A.cbind(C);
  cout << A << "\n";

  B.cbind(B).cbind(B);
  cout << B << "\n";

}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

  // Test SVD;
  Block<TReal> A(3,4);

  // Read in transposed.
  A.scanString( " 4 3 0 "
                " 0 4 1 "
                " 0 1 4 "
                " 3 1 1 " );

  // test_svd(A, 'A');
  // test_svd(A, 'S');

  // test_svd2(A);

  // A.resize(3,4);

  // test_svd(A, 'A', true);
  // test_svd(A, 'S', true);

  // test_svd2(A);

  // test_syrk();

  // test_read();

  // test_out(A);

  // test_copy();

  // test_reshape();

  test_cbind();

  return 0;

}
