// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppArmadillo.h>
using namespace arma;
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

#include <RcppEigen.h>
#include <complex>
#include <vector>
#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>


///////////////////////////////////////////////////////////////////////
////////////////////////// Basic functions ////////////////////////////
///////////////////////////////////////////////////////////////////////

arma::vec resample(int N) {
  // vector permutation function
  
  int i;
  arma::vec out(N);
  arma::uvec tmp = arma::randperm(N);
  for(i=0; i<N; i++){
    out(i) = tmp(i);
  }
  return (out);
} 

arma::mat Null(arma::mat M){
  // Find Null space
  
  int i, j, r;
  mat Q;
  mat R;
  qr(Q, R, M);  // QR decomposition -> M = QR
  int Qcol;
  
  r = arma::rank(M, 1e-14);
  if(r == 0) {
    Qcol = M.n_cols; 
  }
  else {
    Qcol = M.n_rows - r;
  }
  
  arma::mat Q2(M.n_rows, Qcol);   // Q = [Q1 Q2] : concatenate 
  
  for(j=0; j<Qcol; j++){
    for(i=0; i<M.n_rows; i++){
      Q2(i, j) = Q(i, j+r);
    }
  }
  
  return(Q2);
}


// Sampling a sample from startN:Nmax (Natural number) with specific probability(no need normalizing)
int sample_N(int startN, int Nmax, arma::vec prob){
  int i; 
  double N;
  arma::vec Nseq(Nmax-startN+1, fill::zeros);
  
  for(i=startN; i<=Nmax; i++){ 
    Nseq(i-startN) = i; 
  }
  N = Rcpp::RcppArmadillo::sample(Nseq, 1, TRUE, prob)(0);
  
  return((int)N);
}


///////////////////////////////////////////////////////////////////////
///////////// Find roots for cubic and quartic equation. //////////////
///////////////////////////////////////////////////////////////////////
// 3차 방정식 해 구하는 함수
typedef std::complex<double> dcomplex;
cx_vec roots_cubic(const dcomplex a2, const dcomplex a1, const dcomplex a0) {
  dcomplex Q, R, D, S, T;
  
  Q = (dcomplex(3.0)*a1 - a2*a2)/dcomplex(9.0);
  R = (dcomplex(9.0)*a2*a1 - dcomplex(27.0)*a0 - dcomplex(2.0)*std::pow(a2, 3.0))/dcomplex(54.0);
  D = std::pow(Q, 3) + std::pow(R, 2);
  
  S = std::pow( R + std::pow(D, 0.5), 0.3333333333333333); // 왜 때문인지 1/3로하면 계산이 안됨...
  T = std::pow( R - std::pow(D, 0.5), 0.3333333333333333);
  
  arma::cx_vec roots(3);
  roots[0] = -a2/dcomplex(3.0) + S + T;
  roots[1] = -a2/dcomplex(3.0) - dcomplex(0.5)*(S + T) + dcomplex(0, 1)*dcomplex(0.5)*dcomplex(std::pow(3, 0.5), 0)*(S - T);
  roots[2] = -a2/dcomplex(3.0) - dcomplex(0.5)*(S + T) - dcomplex(0, 1)*dcomplex(0.5)*dcomplex(std::pow(3, 0.5), 0)*(S - T);
  
  return(roots);
}

// 4차 방정식 해 구하는 함수
arma::cx_vec roots_quartic(dcomplex a3, dcomplex a2, dcomplex a1, dcomplex a0) {
  
  arma::cx_vec z3 = roots_cubic(-a2, a1*a3-dcomplex(4.0)*a0, dcomplex(4.0)*a2*a0 -a1*a1-a3*a3*a0);
  arma::cx_vec roots(4);        
  dcomplex y, R, D, E;
  int i;
  double tol;
  tol = std::pow(10, -10);
  
  for(i=0; i<z3.n_elem; i++){
    if(std::abs(imag(z3[i])) < tol){
      y = real(z3[i]);
      break;
    }
  }
  
  R = std::pow(dcomplex(0.25)*a3*a3 - a2 + y, 0.5);
  D = std::pow(dcomplex(0.25*3.0)*a3*a3 - R*R - dcomplex(2.0)*a2 + dcomplex(0.25)*(dcomplex(4.0)*a3*a2 - dcomplex(8.0)*a1 - std::pow(a3, 3))/R, 0.5);
  E = std::pow(dcomplex(0.25*3)*a3*a3 - R*R - dcomplex(2.0)*a2 - dcomplex(0.25)*(dcomplex(4.0)*a3*a2 - dcomplex(8)*a1 - std::pow(a3, 3))/R, 0.5);
  
  roots[0] = -a3/dcomplex(4.0) + R/dcomplex(2.0) + D/dcomplex(2.0);
  roots[1] = -a3/dcomplex(4.0) + R/dcomplex(2.0) - D/dcomplex(2.0);
  roots[2] = -a3/dcomplex(4.0) - R/dcomplex(2.0) + E/dcomplex(2.0);
  roots[3] = -a3/dcomplex(4.0) - R/dcomplex(2.0) - E/dcomplex(2.0);
  
  return(roots);
}



///////////////////////////////////////////////////////////////////////
////////////////////////// Sample from vMF ////////////////////////////
///////////////////////////////////////////////////////////////////////

Rcpp::NumericVector rvmf(arma::dvec kmu){
  // Sampling from vmf
  
  const int m = kmu.n_elem;
  Rcpp::NumericVector u(m);
  int i, j;
  
  double b, x0, c, W, Z, V_norm, kap=0.0;
  arma::dmat mu(m, 1, fill::zeros);
  arma::dvec V(m, fill::zeros);
  arma::dmat R(m, m, fill::zeros);
  
  
  if(m == 1) {
    u = std::pow(-1, R::rbinom(1, 1/(1 + std::exp(2.0*kmu(0)))));
  }
  else{
    kap = arma::norm(kmu);
    for(i=0; i<m; i++){
      mu(i, 0) = kmu(i)/kap;
    }
    arma::mat Nullmu = Null(mu);
    
    b = (-2.0*kap + std::sqrt((4.0*kap*kap) + (m-1)*(m-1)*1.0))/((m*(1.0)-1.0));
    x0 = (1.0-b)/(1.0+b);
    c = kap*x0 + ((m-1)*1.0) * std::log(1.0 - x0*x0);
    
    
    // Sample W using Rejection sampling 
    bool done = FALSE;
    while(!done){
      Z = R::rbeta((m-1.0)/2.0, (m-1.0)/2.0);
      W = (1.0-(1.0+b)*Z) / (1.0-(1.0-b)*Z);
      done = (kap*W + (m-1)*log(1.0 - x0*W) - c > log(R::runif(0, 1)) );
    }
    
    // Sample V from the Uniform distribution on the sphere S^(m-1)
    V_norm = 0.0;
    for(i=0; i<(m-1); i++) {
      V(i) = R::rnorm(0, 1);
    }
    for(i=0; i<(m-1); i++) {
      V_norm += V(i)*V(i);
    }
    V_norm = std::sqrt(V_norm);
    for(i=0; i<(m-1); i++) {
      V(i) = std::sqrt(1.0 - W*W)*V(i)/V_norm; 
    }
    V(m-1) = W;       
    
    for(i=0; i<(m-1); i++) {
      for(j=0; j<m; j++) {
        R(j, i) = Nullmu(j, i);
      }
    }
    for(j=0; j<m; j++) {
      R(j, m-1) = mu(j);   
    }  
    
    for(i=0; i<m; i++) {
      u(i) = 0;
    }
    for(i=0; i<m; i++){
      for(j=0; j<m; j++) {
        u(i) += R(i, j)*V(j);     // R%*%V
      }
    }
  }
  return(u);
}

///////////////////////////////////////////////////////////////////////
//////////////// Computing infinite sum part(al, bl) /////////////////
///////////////////////////////////////////////////////////////////////

arma::vec ln2moment(double mu, double s2, int mmax){
  // bl : Computing even moments of Normal distribution(log)
  
  int i;
  arma::dvec lmom(mmax*2+1, fill::zeros);
  arma::dvec l2mom(mmax+1, fill::zeros);
  lmom(1) = std::log(mu);
  
  for(i=2; i<=(mmax)*2; i++){
    lmom(i) = lmom(i-1) + std::log((i-1)*s2*std::exp(lmom(i-2)-lmom(i-1) + mu));
    if(i%2 == 0){
      l2mom(i/2) = lmom(i);
    }
  }
  
  return(l2mom);
}


arma::vec lcr(arma::vec theta, const int lmax){
  // al : Computing for E[lambda*q]
  
  arma::vec lc(lmax+1, fill::zeros);
  arma::vec ls(lmax+1, fill::zeros);
  arma::vec lr(lmax+1, fill::zeros);
  int n = theta.n_elem;
  double tmp;
  int l, k;

  for(l=1; l<lmax+1; l++){
    tmp = 0.0;
    for(k=0; k<n; k++) {
      tmp += std::pow(theta(k)/theta(0), l);
    }
    
    ls(l) = l*std::log(theta(0)) + log(tmp);
    tmp = 0.0;
    for(k=0; k<l; k++){
      tmp += std::exp(lc(k) + ls(l-k) - ls(l) - std::log(2.0*l));
    }
    
    lc(l) = ls(l) + std::log(tmp);
    lr(l) = lc(l) + std::lgamma(l+1) + std::lgamma(n/2.0) - std::lgamma(n/2.0 + l);
  }
  
  return(lr);
}

///////////////////////////////////////////////////////////////////////
////////////////// Sampling d from infinite mixture ///////////////////
///////////////////////////////////////////////////////////////////////
double ldxl(const double x, const double mu, const double sigma, const int l){
  double ln2m = ln2moment(mu,sigma,l)(l);
  
  return((l*1.0)*std::log(x*x) - std::log(sigma) - 0.5*std::log(2.0*M_PI) - 0.5*std::pow((x - mu)/sigma, 2) - ln2m);
}

double rxl(const double mu, const double sigma, const int l, double nu){
  nu = 1.0;
  int i, count=0;
  
  double theta, tau, a, b, c, d, x, lM, lrratio;
  double tol = std::pow(10, -10);
  
  theta = 0.5*mu*(1 + std::sqrt(1.0 + 8.0*(l*1.0)*sigma*sigma/(mu*mu)));
  tau = 1/std::sqrt( 1/sigma*sigma + 2*(l*1.0)/(theta*theta));
  
  
  a = -2.0 * theta - mu;
  b = (theta*theta + 2*mu*theta + nu*tau*tau +sigma*sigma*(-nu - 2*(l*1.0) - 1.0));
  c = (-mu*theta*theta + (nu + 4.0*(l*1.0) + 1.0)*sigma*sigma*theta - mu*nu*tau*tau  );
  d = (-2.0*(l*1.0)*sigma*sigma*theta*theta -2*(l*1.0)*nu*sigma*sigma*tau*tau);
  
  
  cx_vec z4 = roots_quartic(a, b, c, d);
  
  // 실근의 갯수 -> xc(root vector) 크기 지정
  for(i=0; i < z4.n_elem; i++){
    if(std::abs(imag(z4[i])) < tol){
      count += 1;
    }
  }
  arma::vec xc(count);
  
  for(i=0; i < count; i++){
    if(std::abs(imag(z4[i])) < tol){
      xc(i) = real(z4[i]);
    }
  }
  
  arma::vec lM_temp(count, fill::zeros);
  for(i=0; i < count; i++){
    lM_temp(i) = ldxl(xc(i), mu, sigma, l) - ( -std::log(tau) +  R::dt( (xc(i) - theta)/tau, nu, 1));
  }
  lM = max(lM_temp);
  
  bool samp = TRUE;
  while(samp==TRUE){
    x = theta + R::rt(nu)*tau;
    lrratio = ldxl(x, mu, sigma, l) - (-log(tau) + R::dt((x - theta)/tau, nu, 1)) - lM;
    samp = (std::log(R::runif(0,1)) > lrratio);
  }
  
  return(x);
}


///////////////////////////////////////////////////////////////////////
////////////////////////// Step A에서 u, v 뽑기 //////////////////////
///////////////////////////////////////////////////////////////////////
Rcpp::List ruv_A(arma::mat A, const int nsamp=25) {
  int i, j, t;
  int m = A.n_rows;
  int n = A.n_cols;
  
  arma::vec tmp_d;
  arma::mat tmp_u;
  arma::mat tmp_v;
  arma::svd(tmp_u, tmp_d, tmp_v, A);
  
  int s;
  double rb;
  arma::vec v(n, fill::zeros);
  arma::vec u(m, fill::zeros);
  
  // sample a mode
  s = sample_N(1, n, exp( tmp_d - max(tmp_d)) ); 
  rb = R::rbinom(1, 0.5);
  
  for(i=0; i<n; i++){
    v(i) = std::pow(-1, rb)*tmp_v(i, s-1);    
  }
  
  // sample around a mode
  for(t=0; t<nsamp; t++){
    arma::vec Av(m, fill::zeros);
    
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
        Av(i) += A(i, j)*v(j);
      }
    }
    u = rvmf(Av);   
    
    arma::vec Au(n, fill::zeros);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
        Au(j) += A(i, j)*u(i);
      }
    }
    v = rvmf(Au);
  }

  Rcpp::List output;
  output["u"] = u;
  output["v"] = v;  
  return(output);
}

///////////////////////////////////////////////////////////////////////
////                          Main Function                        ////
///////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List BSVD_var(arma::mat Y, arma::mat U, arma::mat V, arma::mat D, const int llb, const int lub,
                   const int nsamp, const int niter, const int nburn, const int nthin, const int nprint,
                   const double mu_pr_precision, const double mu_pr_mean, const double nu0, const double t20, 
                   const double eta0, const double s20){
  
  // Define variables
  int i, l, j, k, t, count, l_samp, col_j;
  int m0, n0; // the number of rows and columns of E_tilde
  count = 0;
  const int Yrow = Y.n_rows;
  const int Ycol = Y.n_cols;
  int nrank = 0;
  int nmcmc = nmcmc;
  double mle;
  
  arma::dmat Umat = U;
  arma::dmat Vmat = V;
  arma::dmat Dmat = D;
  
  arma::dmat UDV(Yrow, Ycol, fill::zeros);
  arma::dmat UDV_j(Yrow, Ycol, fill::zeros);
  
  // Starting values
  double phi=1/s20;
  double mu=mu_pr_mean; 
  double psi=1/t20;
  

  // U, V, E matrix except jth column
  arma::dmat U_j(Yrow, Ycol-1, fill::zeros);
  arma::dmat V_j(Ycol, Ycol-1, fill::zeros);
  arma::dmat E_j(Yrow, Ycol, fill::zeros);
  arma::dmat D_j(Ycol-1, Ycol-1, fill::zeros);
  arma::vec col_ind;

  double d_mu, d_sd, mu_post_mean, mu_post_sd;
  double phi_post_a, phi_post_b, psi_post_a, psi_post_b;

  
  arma::dcube samp_U(nmcmc, Yrow, Ycol, fill::zeros);
  arma::dcube samp_V(nmcmc, Ycol, Ycol, fill::zeros);
  arma::dmat samp_d(nmcmc, Ycol, fill::zeros);
  arma::dvec samp_phi(nmcmc, fill::zeros);
  arma::dvec samp_mu(nmcmc, fill::zeros);
  arma::dvec samp_psi(nmcmc, fill::zeros);
  arma::dvec samp_phi_beta(nmcmc, fill::zeros);
  arma::dvec samp_psi_post_b(nmcmc, fill::zeros);
  arma::dvec samp_mu_post_mean(nmcmc, fill::zeros);
  arma::dvec samp_mu_post_sd(nmcmc, fill::zeros);
  arma::dvec samp_rank(nmcmc, fill::zeros);
  arma::dvec samp_like(nmcmc, fill::zeros);
  
  double tmp, lor, newd_tmp, del, lcrit, odds_dj, lpe_r10, amax, E02;
  int lmax;

  
  for(t=0; t<niter; t++){
    
    ////////////////////////////////////////////////////////
    ///////////// Variable dimension sampler ///////////////
    ////////////////////////////////////////////////////////
    
    col_ind = resample(Ycol);
    for(j=0; j<nsamp; j++){    // j-th column gibbs update
      col_j = col_ind(j);
      
      // // Setting U[, -j] and V[, -j], D[-j, -j]
      for(k=0; k<col_j; k++){             // column index
        for(i=0; i<Yrow; i++) U_j(i, k) = Umat(i, k);
        for(i=0; i<Ycol; i++) V_j(i, k) = Vmat(i, k);
        D_j(k, k) = Dmat(k, k);
      }
      for(k=col_j+1; k<Ycol; k++){
        for(i=0; i<Yrow; i++) U_j(i, k-1) = Umat(i, k);
        for(i=0; i<Ycol; i++) V_j(i, k-1) = Vmat(i, k);
        D_j(k-1, k-1) = Dmat(k, k);
      }
      
      
      // // Construct E_{-j} = Y - U[,-j]D[-j,-j]V[,-j]^T matrix
      arma::mat U_jD_j(Yrow, Ycol-1, fill::zeros);
      for(l=0; l<Ycol-1; l++){
        for(i=0; i<Yrow; i++){
          U_jD_j(i, l) = U_j(i, l)*D_j(l, l);
        }
      }
      
      UDV_j.fill(0.0);
      for(i=0; i<Yrow; i++){
        for(l=0; l<Ycol; l++){
          for(k=0; k<Ycol-1; k++){ 
            UDV_j(i, l) += U_jD_j(i, k)*V_j(l, k);
          }
        }
      }
      
      for(i=0; i<Yrow; i++){
        for(l=0; l<Ycol; l++){
          E_j(i, l) = Y(i, l) - UDV_j(i, l);
        }
      }
    
  
      
      // // Compute Null(-j) matrix
      arma::mat Nu = Null(U_j);
      arma::mat Nv = Null(V_j);
      n0 = Nv.n_cols;
      m0 = Nu.n_cols;


      
      // // Compute E tilde
      arma::mat NuE_j(m0, Ycol, fill::zeros);
      for(i=0; i<m0; i++){
        for(l=0; l<Ycol; l++){
          for(k=0; k<Yrow; k++){ 
            NuE_j(i, l) += Nu(k, i)*E_j(k, l); 
          }
        }
      }
      arma::mat E_tilde(m0, n0, fill::zeros);
      for(i=0; i<m0; i++){
        for(l=0; l<n0; l++){
          for(k=0; k<Ycol; k++){ 
            E_tilde(i, l) += NuE_j(i, k)*Nv(k, l); 
          }
        }
      }

      
      ////////////////////////////////////////////////////////
      ////////////////////// Step A-(1) //////////////////////
      ////////////////////////////////////////////////////////
      
      // Infinite sum part
      // lmax 정하기
      arma::vec z = arma::svd(E_tilde);
      z = arma::pow(z, 2);
      E02 = arma::sum(z);
      arma::vec theta = z/E02;
      amax = z(0) * phi*phi *E02 * max(theta);
      

      
      lcrit = 0.5*(std::sqrt( std::pow((n0*1.0)/2.0 + 1.0, 2.0) + 4*(amax - (n0*1.0)/2) ) - (n0/2 + 1) );
      lcrit = 1.25*lcrit;
      if(lcrit < llb) lcrit = llb; 
      if(lcrit > lub) lcrit = lub;
      lmax = 2*std::floor(lcrit/2);
    
      
      
      // Compute al, bl
      arma::dvec la = lcr(theta, lmax); 
      for(l=0; l<lmax+1; l++){
        la(l) += lgamma((m0*1.0)/2.0) - std::lgamma((m0*1.0)/2.0 + (l*1.0)) - std::lgamma((l*1.0)+1.0) - (l*1.0)*std::log(4.0);
      }
      
      arma::dvec lb(lmax+1, fill::zeros);
      if(mu != 0){
        lb = ln2moment(mu*psi/(phi+psi), 1.0/std::sqrt(psi+phi), lmax) + 0.5*std::log(psi/(phi+psi)) - 0.5*mu*mu*psi*phi/(phi+psi);
      }
      else{
        for(l=0; l<lmax+1; l++){
          lb(l) = 0.5*std::log(psi/(phi+psi)) + (l*1.0)*std::log(1.0/(phi+psi)) + ( std::lgamma(2.0*l+1.0) - std::lgamma((l*1.0)+1.0) - (l*1.0)*std::log(2.0) );
        }
      }
      // Rprintf("error in Step A 3 - %d\n", j);
      
      arma::dvec lc(lmax+1, fill::zeros);
      arma::dvec lt(lmax+1, fill::zeros);

      // ||E_tilde||^(2l)
      for(l=0; l<lmax+1; l++){
        lc(l) = (l*1.0)*std::log(E02*phi*phi);
      }
      for(l=0; l<lmax+1; l++){
        lt(l) = la(l) + lb(l) + lc(l);
      }
      
      // A-(1) 에서 summation part 계산 : lpe_r10
      lpe_r10 = max(lt) + log(arma::sum(exp(lt - max(lt))));
      
      
      
      // Sample dj from ({dj!=0}, {dj=0})

      double nonzero_d = 0.0;
      for(i=0; i<Ycol-1; i++){
        if(D_j(i,i) != 0){
          nonzero_d += 1.0;
        }
      }
      lor = log( (nonzero_d + 1.0)/(Ycol*1.0 - nonzero_d) );
      
      odds_dj = 1.0/(1.0 + std::exp(-(lor + lpe_r10)));
      newd_tmp = R::rbinom(1, odds_dj);
      
      
      
      // {dj!=0} 이면 Sample dj from infinite mixture, Sample(u[,j], v[,j]) from joint dist
      arma::vec pdl(lmax+1, fill::zeros);
      arma::vec lvec(lmax+1, fill::zeros);
      
      arma::dvec uv_samp_u(m0, fill::zeros);
      arma::dvec uv_samp_v(n0, fill::zeros);
      del = 0;
      
      if(newd_tmp == 1){
        // sample del (d element)
        for(l=0; l<lmax+1; l++) {pdl(l) = std::exp(lt(l) - max(lt)); }
        l_samp = sample_N(0, lmax, pdl);
        if(l_samp > 0 & mu!=0) del = rxl( mu*psi/(psi+phi), std::sqrt(1.0/(phi+psi)), l_samp, 1 );
        if(l_samp > 0 & mu==0) del = std::sqrt( R::rchisq(2.0*l_samp + 1.0)/(phi + psi) );
        if(l_samp == 0) del = R::rnorm(mu*psi/(psi + phi), std::sqrt( 1/(phi + psi)));
        // Rprintf("error in Step A 4 - %d\n", j);
        
        // Sample (U[,j], V[,j])
        Rcpp::List uv_samp = ruv_A(E_tilde*phi*del);
        arma::dvec uv_samp_u = uv_samp["u"];
        arma::dvec uv_samp_v = uv_samp["v"];
        
        for(i=0; i<Yrow; i++){
          tmp = 0.0;
          for(k=0; k<m0; k++){ tmp += Nu(i, k)*uv_samp_u(k);}
          Umat(i, col_j) = tmp;
        }
        for(l=0; l<Ycol; l++){
          tmp = 0.0;
          for(k=0; k<n0; k++){ tmp += Nv(l, k)*uv_samp_v(k); }
          Vmat(l, col_j) = tmp;
        }
        Dmat(col_j, col_j) = del;
      }
      else{
        for(i=0; i<Yrow; i++){Umat(i, col_j) = 0.0;}
        for(l=0; l<Ycol; l++){Vmat(l, col_j) = 0.0;}
        Dmat(col_j, col_j) = 0.0;
      }
      // Rprintf("Check del : %f\n", del);

    } // end Step A(iteration for j (column))
  
    

    
    
    
    
    ////////////////////////////////////////////////////////
    //////////////////////// Step B ////////////////////////
    ////////////////////////////////////////////////////////
    // Step B : Fixed dimension sampler
    
    
    // Update U[, j], V[, j], D[j, j]
    col_ind = resample(Ycol);
    nrank = 0;
    for(j=0; j<Ycol; j++){    // j-th column gibbs update
      col_j = col_ind(j);
      
      if(Dmat(col_j, col_j) != 0){          // d_j 가 0이 아닐때
        nrank += 1;
        // // Construct U[, -j] and V[, -j]
        for(k=0; k<col_j; k++){             // column index
          for(i=0; i<Yrow; i++) U_j(i, k) = Umat(i, k);
          for(i=0; i<Ycol; i++) V_j(i, k) = Vmat(i, k);
          D_j(k, k) = Dmat(k, k);
        }
        for(k=col_j+1; k<Ycol; k++){
          for(i=0; i<Yrow; i++) U_j(i, k-1) = Umat(i, k);
          for(i=0; i<Ycol; i++) V_j(i, k-1) = Vmat(i, k);
          D_j(k-1, k-1) = Dmat(k, k);
        }
        
        
        
        // // Construct E(-j) matrix
        arma::mat U_jD_j(Yrow, Ycol-1, fill::zeros);
        for(l=0; l<Ycol-1; l++){
          for(i=0; i<Yrow; i++){
            U_jD_j(i, l) = U_j(i, l)*D_j(l, l);
          }
        }
        
        UDV_j.fill(0.0);
        for(i=0; i<Yrow; i++){
          for(l=0; l<Ycol; l++){
            for(k=0; k<Ycol-1; k++){ 
              UDV_j(i, l) += U_jD_j(i, k)*V_j(l, k);
            }
          }
        }
        
        for(i=0; i<Yrow; i++){
          for(l=0; l<Ycol; l++){
            E_j(i, l) = Y(i, l) - UDV_j(i, l);
          }
        }
        
        
        
        // // Compute Null(-j) matrix
        arma::mat Nu = Null(U_j); 
        arma::mat Nv = Null(V_j); 
        n0 = Nv.n_cols;
        m0 = Nu.n_cols;
        // Rprintf("error in Step B 2 - %d\n", j);
        
        
       
        // // Update u[,j] and v[,j]
   
        // Sample vj from vMF
        arma::mat NvE_j(n0, Yrow, fill::zeros);
        for(i=0; i<n0; i++){
          for(l=0; l<Yrow; l++){
            for(k=0; k<Ycol; k++){ 
              NvE_j(i, l) += Nv(k, i)*E_j(l, k); 
            }
          }
        }
        arma::dvec v_param(n0, fill::zeros);
        for(i=0; i<n0; i++){
          for(k=0; k<Yrow; k++){ 
            v_param(i) += NvE_j(i, k)*Umat(k, col_j); 
          }
        }
        v_param = phi*Dmat(col_j, col_j)*v_param;
        arma::vec vj = rvmf(v_param);
        
        // Update V[,j]
        for(l=0; l<Ycol; l++){
          Vmat(l, col_j) = 0.0;
          for(k=0; k<n0; k++){ 
            Vmat(l, col_j) += Nv(l, k)*vj(k); 
          }
        }
        
        
        // Sample uj from vMF
        arma::mat NuE_j(m0, Ycol, fill::zeros);
        for(i=0; i<m0; i++){
          for(l=0; l<Ycol; l++){
            for(k=0; k<Yrow; k++){ 
              NuE_j(i, l) += Nu(k, i)*E_j(k, l); 
            }
          }
        }
        arma::dvec u_param(m0, fill::zeros);
        for(i=0; i<m0; i++){
          for(k=0; k<Ycol; k++){ 
            u_param(i) += NuE_j(i, k)*Vmat(k, col_j); 
          }
        }
        u_param = phi*Dmat(col_j, col_j)*u_param;
        arma::vec newu = rvmf(u_param);
        
        // Update U[,j]
        for(i=0; i<Yrow; i++){
          Umat(i, col_j) = 0.0;
          for(k=0; k<m0; k++){ 
            Umat(i, col_j) += Nu(i, k)*newu(k);
          }
        }
        
        // Rprintf("error in Step B 3 - %d\n", j);
        
        
        
        
        // Update D[j, j]
        arma::vec UjE_j(Ycol, fill::zeros);
        for(i=0; i<Ycol; i++){
          for(k=0; k<Yrow; k++){ 
            UjE_j(i) += Umat(k, col_j)*E_j(k, i); 
          }
        }
        double UjE_jVj = 0.0;
        for(k=0; k<Ycol; k++){ 
          UjE_jVj += UjE_j(k)*Vmat(k, col_j); 
        }
        
        d_mu = (UjE_jVj*phi + mu*psi)/(phi + psi);
        d_sd = std::sqrt(1/(phi + psi));
        
        // Sampling d(j)
        Dmat(col_j, col_j) = R::rnorm(d_mu, d_sd);
        
      } // End iteration for update columns of  U, D, V when D[j,j]!=0
    }
      
      
      
      
      
      
      

    
    ////////////////////////////////////////////////////////
    //////////////////////// Step C ////////////////////////
    ////////////////////////////////////////////////////////
    // Step C : Other terms
   
    // Update phi
    
    // Compute UDV matrix
    arma::mat UD(Yrow, Ycol, fill::zeros);
    for(l=0; l<Ycol; l++){
      for(i=0; i<Yrow; i++){ 
        UD(i, l) = Umat(i, l)*Dmat(l, l); 
      }
    }
    
    UDV.fill(0.0);
    for(i=0; i<Yrow; i++){
      for(l=0; l<Ycol; l++){
        for(k=0; k<Ycol; k++){ 
          UDV(i, l) += UD(i, k)*Vmat(l, k); 
        }
      }
    }
    
    double Y_UDV_norm = 0.0;
    for(i=0; i<Yrow; i++){
      for(l=0; l<Ycol; l++){
        Y_UDV_norm += pow(Y(i ,l) - UDV(i, l), 2.0);
      }
    }
    
    phi_post_b = (nu0*s20 + Y_UDV_norm)/2.0;
    phi_post_a = (nu0 + Yrow*Ycol)/2.0;
    
    // Sampling
    phi = R::rgamma(phi_post_a, 1.0/phi_post_b);
    
    
    
    // Update mu
    double sum_dj = 0.0;
    for(j=0; j<Ycol; j++) {
      sum_dj += Dmat(j, j);
    }
    
    mu_post_mean = (psi*tmp + mu_pr_mean*mu_pr_precision)/(psi*(nrank*1.0) + mu_pr_precision);
    mu_post_sd = std::sqrt(1/(psi*(nrank*1.0) + mu_pr_precision));
    
    // Sampling
    mu = R::rnorm(mu_post_mean, mu_post_sd);
    
    
    
    
    // Update psi
    double dj_mu_diff2 = 0.0;
    for(j=0; j<Ycol; j++) {
      dj_mu_diff2 += std::pow(Dmat(j, j) - mu, 2.0);   
    }
    psi_post_b = (eta0*t20 + dj_mu_diff2)/2.0;
    psi_post_a = (eta0 + nrank*1.0)/2.0;

    // Sampling
    psi = R::rgamma(psi_post_a, 1.0/psi_post_b);
    
    
    
    
    
    
    
    
    
    
    // Store variables
    // Store samples with burning and thinning
    if(t >= nburn) {
      if(t % nthin == 0){
        
        for(i=0; i<Yrow; i++){
          for(j=0; j<Ycol; j++){
            samp_U(count, i, j) = Umat(i, j);
          }
        }
        for(l=0; l<Ycol; l++){
          for(j=0; j<Ycol; j++){
            samp_V(count, l, j) = Vmat(l, j);
          }
        }
        for(j=0; j < Ycol; j++) {
          samp_d(count, j) = Dmat(j, j);
        }
        samp_phi(count) = phi;
        samp_mu(count) = mu;
        samp_psi(count) = psi;
        samp_phi_beta(count) = phi_post_b;
        samp_psi_post_b(count) = psi_post_b;
        samp_mu_post_mean(count) = mu_post_mean;
        samp_mu_post_sd(count) = mu_post_sd;
        samp_rank(count) = nrank;
        
        
        mle = 0.0;
        // likelihood
        mle += (Yrow*Ycol)*1.0 * std::log(phi) - phi/2.0*Y_UDV_norm;
        
        // prior : d, phi, mu, psi
        for(k=0; k<nrank; k++){
          mle += R::dnorm4(Dmat(k, k), mu, std::sqrt(1/psi), 1);
        }
        mle += R::dgamma(phi, nu0/2.0, nu0*s20/2.0, 1);
        mle += R::dgamma(psi, eta0/2.0, eta0*s20/2.0, 1);
        mle += R::dnorm4(mu, mu_pr_mean, std::sqrt(1/mu_pr_precision), 1);
        
        samp_like(count) = mle;

        
        count++;
      }
      
    }// end burning and thinning
    if(t % nprint == 0) {
      Rprintf("| iter - %d ", t);
      Rprintf("niter %d phi %f mu %f psi %f phi_post_b %f psi_post_b %f nrank %d \n",t, phi, mu, psi, phi_post_b, psi_post_b, nrank);
    }
  
  }
  Rcpp::List output;    
  output["U"] = samp_U;
  output["D"] = samp_d;
  output["V"] = samp_V;
  output["phi"] = samp_phi;
  output["mu"] = samp_mu;
  output["psi"] = samp_psi;
  output["phi_post_b"] = samp_phi_beta;
  output["psi_post_b"] = samp_psi_post_b;
  output["mu_post_mean"] = samp_mu_post_mean;
  output["mu_post_sd"] = samp_mu_post_sd;
  output["rank"] = samp_rank;
  output["map"] = samp_like;
  
  return(output);
  
  
} // end function
  
  