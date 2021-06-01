//
// Created by kazem on 2020-05-18.
//

#ifndef NASOQ_LBL_EIGEN_H
#define NASOQ_LBL_EIGEN_H

#include <cholesky_solver.h>
//#include <parsy/common/def.h>
#include <Eigen/Core>
#include <Eigen/Sparse>


 namespace sympiler
 {
  //  Solving Ax=b
  // Inputs:
  //   H  n by n sparse Hessian matrix **lower triangle only** (see
  //     .triangularView<Eigen::Lower>() )
  //   q  n by 1 vector
  // Outputs:
  //   x  n by 1 solution vector
  // Returns sympiler exit flag
  //
  int sympiler_spd_linsolve(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A,
    Eigen::Matrix<double,Eigen::Dynamic,1> b,
    Eigen::Matrix<double,Eigen::Dynamic,1> & x);

  sym_lib::parsy::CSC* sympiler_csc_format(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A);

  sym_lib::parsy::SolverSettings* sympiler_chol_symbolic(
    // Pass inputs by copy so we get non-const and casted data
    sym_lib::parsy::CSC *A);


  sym_lib::parsy::SolverSettings* sympiler_chol_symbolic(
    // Pass inputs by copy so we get non-const and casted data
    sym_lib::parsy::CSC *A, int num_threads, int mode = 4);

 }



#endif //NASOQ_LBL_EIGEN_H
