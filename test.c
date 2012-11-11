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

// typedef float Real;
typedef double Real;

int main(int argc, char** argv)
{

  // vector<double> a(100);
  // Frame<double> aframe(&a[0], 10, 10);
  // vector<double> b(100);
  // Frame<double> bframe(&b[0], 10, 10);

  // operator+=(aframe, bframe);

  // Matrix A("I", 3);
  // cout << A.write("A.mat", true) << "\n";

  // Matrix A;
  // A.read("A.mat", true);
  // cout << A;

  Block<Real> A(3, 3);
  // string av = 
  //   " 1 2 3 " 
  //   " 4 5 6 " 
  //   " 7 8 9 " ;
  
  // stringstream ss(av);

  // A << " 1 2 3 "
  //      " 4 5 6 "
  //      " 7 8 9 " ;

  // cout << A;

  // Matrix B = 5 * Matrix("I", 3);

  // Matrix C(0.0, B);

  // cout << B ;

  // cout << hmin(B, A);

  // MatrixFrame C(&A(0), 3, 1);
  // Matrix B("I", 2);
  // B = (A + A) + 1 + C;
  // cout << B;

  // Matrix B(3, 1);
  // string b = " 1 2 3 ";
  // B.readstring(b);
  // A -= B;
  // cout << A + A;

  //   Matrix rows(2); rows << "0 1";
  //   Matrix cols(2); cols << "1 2";
  //   Matrix B(2,2);
  //   B.copy(A, rows, cols);
  //   cout << "A:\n" << A << "\n";
  //   cout << "B:\n" << B << "\n";

  A << " 4 0 0 "
       " 0 4 1 "
       " 0 1 4 ";

  cout << "A:\n" << A << "\n";
  
  Block<Real> AA(3, 3); AA.fill(1.0);
  syrk(AA, A, 'N', (Real)1.0, (Real)0.0);

  cout << "AA:\n" << AA << "\n";

  // CONJUGATE GRADIENT
  Block<Real> b("N", 3);

  Block<Real> x(b);
  symsolve(A, x);

  cout << "x:\n" << x;

  Block<Real> y(b);
  try {
  cout << "iter: " << cg(y, A, b, (Real)10e-8, 100) << "\n";
  }
  catch (std::exception& e) {
    cout << e.what();
  }

  cout << "y:\n" << y;

  // SCAN
  // Matrix b(2,2);
  // b.scan("b.dat");
  // cout << "b:\n" << b;

  // // EIGENVECTOR DECOMP OF SYM MATRIX
  // Matrix evec;
  // Matrix eval;
  // Matrix rt;
  // Matrix irt;

  // symeigen(evec, eval, A);
  // symsqrt(rt, A);
  // syminvsqrt(irt, A);

  // cout << "A:\n" << A;
  // cout << "U:\n" << evec;
  // cout << "D:\n" << eval;

  // cout << "rt:\n" << rt;
  // cout << "irt:\n" << irt;

  // cout << "Max: " << maxAll(A) << "\n";
  // cout << "Min: " << minAll(A) << "\n";

  // // SVD
  // A << " 4 3 0 "
  //      " 0 4 1 "
  //      " 0 1 4 ";

  // Matrix U;
  // Matrix S;
  // Matrix tV;

  // svd(U, S, tV, A);
  // cout << "A:\n" << A;
  // cout << "U:\n" << U;
  // cout << "S:\n" << S;
  // cout << "tV:\n" << tV;

  //Block<int> C("W", 3);
  //cout << "C:\n" << C << "\n";

  return 0;

}
