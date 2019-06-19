#include <iostream>
#include <string>
#include <float.h>
#include <limits>
#include "MPRGPSolver.h"
using namespace std;
using namespace MATH;

typedef FixedSparseMatrix<double> MAT;

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

int main(int argc, char *argv[]){

  // init QP problem
  SparseMatrix<double> A, J;
  VectorXd B, c, x;
  const string qp_file = "./data/qp33.b";
  const bool load_qp_succ = loadQP(A, B, J, c, x, qp_file);
  assert(load_qp_succ);

  cout << "dimension: " << A.rows() << endl;
  cout << "constraints: " << J.rows() << endl;

  // solve
  MAT FA(A);
  const int rlst_code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, x);

  // check results
  assert(rlst_code == 0);
  const VectorXd cx = J*x-c;
  for (int i = 0; i < cx.size(); ++i){
	assert(cx[i] >= 0);
  }
  
  return 0;
}
