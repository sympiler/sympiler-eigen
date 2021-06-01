//
// Created by kazem on 2020-05-18.
//


#include "includes/sympiler_cholesky.h"
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
    Eigen::Matrix<double,Eigen::Dynamic,1> & x){
   assert(A.isApprox(A.triangularView<Eigen::Lower>(),0) &&
          "P should be lower triangular");
   assert(A.isCompressed());
   assert(A.rows()==A.cols());
   assert(A.rows()==b.rows());

   auto *H = new sym_lib::parsy::CSC;
   H->nzmax = A.nonZeros();
   H->ncol= H->nrow=A.rows();
   H->p = A.outerIndexPtr();
   H->i = A.innerIndexPtr();
   H->x = A.valuePtr();
   H->stype=-1;
   H->xtype=SYMPILER_REAL;
   H->packed=TRUE;
   H->nz = NULL;
   H->sorted = TRUE;

   auto *sym_chol = new sym_lib::parsy::SolverSettings(H,b.data());
   sym_chol->ldl_variant = 1;
   sym_chol->solver_mode = 0;
   sym_chol->symbolic_analysis();
   sym_chol->numerical_factorization();
   double *sol = sym_chol->solve_only();
   //lbl->compute_norms();
   //std::cout<<"residual: "<<lbl->res_l1<<"\n";

   x = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
     sol,A.rows(),1);

   // Exitflag TODO
   int exitflag = 0;
   delete sym_chol;
   delete H;
   return exitflag;
  }


  sym_lib::parsy::CSC* sympiler_csc_format(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A) {
   assert(A.isApprox(A.triangularView<Eigen::Lower>(), 0) &&
          "P should be lower triangular");
   assert(A.isCompressed());
   assert(A.rows() == A.cols());

   auto *H = new sym_lib::parsy::CSC;
   H->nzmax = A.nonZeros();
   H->ncol = H->nrow = A.rows();
   H->p = A.outerIndexPtr();
   H->i = A.innerIndexPtr();
   H->x = A.valuePtr();
   H->stype = -1;
   H->xtype = SYMPILER_REAL;
   H->packed = TRUE;
   H->nz = NULL;
   H->sorted = TRUE;
   return H;
  }

  sym_lib::parsy::SolverSettings* sympiler_chol_symbolic(
    // Pass inputs by copy so we get non-const and casted data
    sym_lib::parsy::CSC *A){

   auto *sym_chol = new sym_lib::parsy::SolverSettings(A);
   sym_chol->ldl_variant = 1;
   sym_chol->solver_mode = 0;
   sym_chol->symbolic_analysis();
   return sym_chol;
  }


  sym_lib::parsy::SolverSettings* sympiler_chol_symbolic(
    // Pass inputs by copy so we get non-const and casted data
    sym_lib::parsy::CSC *A, int num_threads, int mode){

   auto *sym_chol = new sym_lib::parsy::SolverSettings(A);
   sym_chol->ldl_variant = mode;
   sym_chol->num_thread = num_threads;
   sym_chol->solver_mode = 0;
   sym_chol->symbolic_analysis();
   return sym_chol;
  }

 }



