//
// Created by kazem on 2020-05-18.
//

#include <unsupported/Eigen/SparseExtra>
#include <sympiler_cholesky.h>
#include "../sympiler/utils/includes/def.h"

int main(int argc, const char *argv[]){

 std::string h_path = argv[1];
 /// Reading input matrices.
 Eigen::SparseMatrix<double,Eigen::ColMajor,int> H;
 std::string message = "Could not load ";
 if( !Eigen::loadMarket( H, h_path ) ){ std::cout<<message<<"H"; return 1; }
 //if( !Eigen::loadMarketVector( q, rhs_path ) ){ std::cout<<message<<"q"; return 1; }
 Eigen::VectorXd q = Eigen::VectorXd::Ones(H.cols());

  //std::cout<< q;
 /// output vectors
 Eigen::VectorXd x, x1;

 ///A1: call the linear solver wrapper.
 int ret=0;
 ret = sympiler::sympiler_spd_linsolve(H,q,x);
 //std::cout<<x;

 ///A2: Call each step separately.
 Eigen::VectorXd q1 = Eigen::VectorXd::Ones(H.cols());
 double *rhs_q = q1.data();

 auto *A = new sym_lib::parsy::CSC;
 A->nzmax = H.nonZeros();
 A->ncol = A->nrow = H.rows();
 A->p = H.outerIndexPtr();
 A->i = H.innerIndexPtr();
 A->x = H.valuePtr();
 A->stype = -1;
 A->xtype = SYMPILER_REAL;
 A->packed = TRUE;
 A->nz = NULL;
 A->sorted = TRUE;
 timing_measurement t;
 t.start_timer();
 auto *sym_chol1 = sympiler::sympiler_chol_symbolic(A, rhs_q);
 t.measure_elapsed_time();
 std::cout<<"t: "<<t.elapsed_time<<"\n";
 sym_chol1->numerical_factorization();
 auto *sol = sym_chol1->solve_only();
 sym_chol1->compute_norms();
 x1 = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> >(
   sol,H.rows(),1);
 //std::cout<<"\n\n "<<x1;
 std::cout<<"res norm" << sym_chol1->res_l1<<"\n";

 delete sym_chol1;
 delete A;
 delete []sol;
 return ret;
}
