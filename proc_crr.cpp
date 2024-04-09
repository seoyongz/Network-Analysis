// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
using namespace arma;
#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// The function of Procrustes matching
// [[Rcpp::export]]
arma::mat proc_crr(arma::mat Z, arma::mat Z0){
  
  int i, j, r;
  int nrow = Z.n_rows;
  int ncol = Z.n_cols;
  
  
  // Compute A = t(Z)%*%(  Z0%*%t(Z0)  )%*%Z
  //// Z0tZ0 = Z0%*%t(Z0)
  arma::mat Z0tZ0(nrow, nrow, fill::zeros);
  for(i=0; i<nrow; i++){
    for(j=0; j<nrow; j++){
      for(r=0; r<ncol; r++){
        Z0tZ0(i, j) += Z0(i, r) * Z0(j, r);
      }
    }
  }
  
  //// ZZ0tZ0 = t(Z)%*%(  Z0%*%t(Z0)  )
  arma::mat ZZ0tZ0(ncol, nrow, fill::zeros);
  for(i=0; i<ncol; i++){
    for(j=0; j<nrow; j++){
      for(r=0; r<nrow; r++){
        ZZ0tZ0(i, j) += Z(r, i) * Z0tZ0(r, j);
      }
    }
  }
  
  //// A = t(Z)%*%(  Z0%*%t(Z0)  )%*%Z
  arma::mat A(ncol, ncol, fill::zeros);
  for(i=0; i<ncol; i++){
    for(j=0; j<ncol; j++){
      for(r=0; r<nrow; r++){
        A(i, j) = ZZ0tZ0(i, r) * Z(r, j);
      }
    }
  }
  
  // arma::mat eigvec;
  // arma::vec eigval;
  // arma::eig_sym(eigval, eigvec, A);
  
  // // Compute Ahalf = eigvec%*%diag(sqrt(eigval))%*%t(eigvec)
  // 
  // // eigvec%*%diag(sqrt(eigval))
  // arma::mat mat_tmp1(ncol, ncol, fill::zeros);
  // for(i=0; i<ncol; i++){
  //   for(j=0; j<ncol; j++){
  //     mat_tmp1(i, j) = eigvec(i, j) * sqrt(eigval(j));
  //   }
  // }
  // 
  // arma::mat Ahalf(ncol, ncol, fill::zeros);
  // for(i=0; i<ncol; i++){
  //   for(j=0; j<ncol; j++){
  //     for(r=0; r<ncol; r++){
  //       Ahalf(i, j) += mat_tmp1(i, r) * eigvec(j, r);
  //     }
  //   }
  // }
  // arma::mat Ahalf_inv = arma::inv(Ahalf_inv);
  // 
  // 
  // // Compute Zproc = t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z))
  // // Zproc = Z %*% t(solve(Ahalf))%*% t(Z) %*% Z0
  // 
  // // mat_tmp2 = Z %*% t(solve(Ahalf))
  // arma::mat mat_tmp2(nrow, ncol, fill::zeros);
  // for(i=0; i<nrow; i++){
  //   for(j=0; j<ncol; j++){
  //     for(r=0; r<ncol; r++){
  //       mat_tmp2(i, j) += Z(i, r) * Ahalf_inv(j, r);
  //     }
  //   }
  // }
  // 
  // // mat_tmp3 = Z %*% t(solve(Ahalf))%*% t(Z) = mat_tmp2 %*% t(Z)
  // arma::mat mat_tmp3(nrow, nrow, fill::zeros);
  // for(i=0; i<nrow; i++){
  //   for(j=0; j<nrow; j++){
  //     for(r=0; r<ncol; r++){
  //       mat_tmp3(i, j) += mat_tmp2(i, r) * Z(j, r);
  //     }
  //   }
  // }
  // 
  // // Zproc = Z %*% t(solve(Ahalf))%*% t(Z) %*% Z0 = mat_tmp3 %*% Z0
  // arma::mat Zproc(nrow, ncol, fill::zeros);
  // for(i=0; i<nrow; i++){
  //   for(j=0; j<ncol; j++){
  //     for(r=0; r<nrow; r++){
  //       Zproc(i, j) += mat_tmp3(i, r) * Z0(r, j);
  //     }
  //   }
  // }
  // 
  // return(Zproc);
  return(A);
  
}

