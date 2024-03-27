// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;

arma::mat l2_dist(arma::mat z){
  
  const int nresp = z.n_rows;
  const int ndim = z.n_cols;
  int i, j, k;
  double dist_tmp;
  arma::mat dist_mat(nresp, nresp);
  
  for(i=0; i < nresp; i++){
    for(j=i; j < nresp; j++){
      dist_tmp = 0.0;
      for(k=0; k < ndim; k++){
        dist_tmp += std::pow(z(i, k) - z(j, k), 2.0);
      }
      dist_tmp = sqrt(dist_tmp);
      dist_mat(i, j) = -dist_tmp;
      dist_mat(j, i) = -dist_tmp;
    }
  }
  
  return(dist_mat);
}

arma::mat l1_dist(arma::mat z){
  
  const int nresp = z.n_rows;
  const int ndim = z.n_cols;
  int i, j, k;
  double dist_tmp;
  arma::mat dist_mat(nresp, nresp);
  
  for(i=0; i < nresp; i++){
    for(j=i; j < nresp; j++){
      dist_tmp = 0.0;
      for(k=0; k < ndim; k++){
        dist_tmp += std::abs(z(i, k)-z(j, k));
      }
      dist_mat(i, j) = -dist_tmp;
      dist_mat(j, i) = -dist_tmp;
    }
  }
  
  return(dist_mat);
}

arma::mat l_projection(arma::mat z){
  
  const int nresp = z.n_rows;
  const int ndim = z.n_cols;
  int i, j, k;
  double dist_tmp;
  double norm_tmp;
  arma::mat dist_mat(nresp, nresp);
  
  for(i=0; i < nresp; i++){
    for(j=0; j < nresp; j++){
      dist_tmp = 0.0;
      norm_tmp = 0.0;
      for(k=0; k < ndim; k++){
        dist_tmp += z(i, k) * z(j, k);
        norm_tmp += std::pow(z(j, k), 2.0);
      }
      dist_mat(i, j) = dist_tmp/norm_tmp;
    }
  }
  
  return(dist_mat);
}


// [[Rcpp::export]]
Rcpp::List LSM_my(arma::mat data, const int dist_type, const int ndim, 
                  const int niter, const int nburn, const int nthin, const int nprint,
                  const double jump_alpha, const double jump_z, 
                  const double alpha_pr_mean, const double alpha_pr_sd, const double z_pr_mean, const double z_pr_sd){
  
  // dist_type = 1 : l1 norm
  // dist_type = 2 : l2 norm
  // dist_type = 3 : projection model
  
  // Define the type and size of variables
  int t, i, j, k, accept, count;
  accept = count = 0;
  const int nresp = data.n_rows;  // the number of people in the network data
  const int nmcmc = (niter-nburn)/nthin;
  
  double oldalpha, newalpha, mle;
  double old_like, new_like;
  
  double u, ratio, accept_alpha=0.0, accept_z=0.0;
  
  arma::dvec samp_alpha(nmcmc, fill::zeros);
  arma::dcube samp_z(nmcmc, nresp, ndim, fill::zeros);
  arma::dvec sample_mle(nmcmc, fill::zeros);
  
  // Initial value of parameters
  oldalpha = R::rnorm(alpha_pr_mean, alpha_pr_sd);
  arma::dmat oldz(nresp, ndim, fill::randu);
  // arma::dmat oldz = Z0;
  arma::dmat newz = oldz;
  arma::dmat olddist(nresp, nresp);
  
  if(dist_type == 1){
    olddist = l1_dist(oldz);
  }
  else if(dist_type == 2){
    olddist = l2_dist(oldz);
  }
  else if(dist_type == 3){
    olddist = l_projection(oldz);
  }

  arma::dmat newdist = olddist;
  
  
  for(t=0; t < niter; t++){

    ////// Update Z with MH algorithm //////
    for(i=0; i < nresp; i++){
      for(k=0; k < ndim; k++) {
        // Propose z* from symmetric proposal distribution
        newz(i, k) = R::rnorm(oldz(i, k), jump_z);    
      }
    }
    
    if(dist_type == 1){
      newdist = l1_dist(newz);
    }
    else if(dist_type == 2){
      newdist = l2_dist(newz);
    }
    else if(dist_type == 3){
      newdist = l_projection(newz);
    }
    
    
    old_like = new_like = 0.0;
    for(i=0; i<nresp; i++){
      for(j=0; j < nresp; j++){
        if(j != i){
          // Calculate likelihood
          old_like += data(i, j) * (oldalpha + olddist(i, j)) - log(1 + exp(oldalpha + olddist(i, j)));
          new_like += data(i, j) * (oldalpha + newdist(i, j)) - log(1 + exp(oldalpha + newdist(i, j)));
        }
      }
    }
      
    //Rprintf("Error03-%d\n", t);
    ratio = new_like - old_like;
    for(i=0; i<nresp; i++){
      for(k=0; k < ndim; k++){
        ratio += R::dnorm4(newz(i, k), z_pr_mean, z_pr_sd, 1) - R::dnorm4(oldz(i, k), z_pr_mean, z_pr_sd, 1);
      }
    }
    

    // Accept or Reject
    if(ratio > 0.0) accept = 1;
    else{
      u = std::log(R::runif(0, 1));
      if(ratio > u) accept = 1;
      else accept = 0;
    }
    
    if(accept==1){
      oldz = newz;
      olddist = newdist;
      accept_z += 1.0 / (1.0 * niter);
    }
    //Rprintf("Error04-%d\n", t);
      
    ////// Update alpha with MH algorithm //////
    newalpha = R::rnorm(oldalpha, jump_alpha);    // Propose alpha* from symmetric proposal distribution
    // Calculate likelihood
    old_like = new_like = 0.0;
    for(i=0; i < nresp; i++){
      for(j=0; j < nresp; j++){
        if(i!=j){
          old_like += data(i, j) * (oldalpha + olddist(i, j)) - log(1 + exp(oldalpha + olddist(i, j)));
          new_like += data(i, j) * (newalpha + olddist(i, j)) - log(1 + exp(newalpha + olddist(i, j)));
        }
      }
    }
    //Rprintf("Error01-%d\n", t);
    
    // Accept or Reject
    ratio = new_like - old_like;
    ratio += R::dnorm4(newalpha, alpha_pr_mean, alpha_pr_sd, 1) - R::dnorm4(oldalpha, alpha_pr_mean, alpha_pr_sd, 1);
    
    if(ratio > 0) accept = 1;
    else{
      u = std::log(R::runif(0, 1));
      if(ratio > u) accept = 1;
      else accept = 0;
    }
    
    if(accept==1){
      oldalpha = newalpha; 
      accept_alpha += 1.0 / (1.0 * niter);
    }
    //Rprintf("Error02-%d\n", t);
    
    // Burning, Thinning and Store samples
    if(t >= nburn && t % nthin == 0){
      // Store samples
      samp_alpha(count) = oldalpha;
      for(i = 0; i < nresp; i++){
        for(k = 0; k < ndim; k++){
          samp_z(count, i, k) = oldz(i, k);
        }
      }
      
      // Calculate the posterior value
      mle = R::dnorm4(oldalpha, alpha_pr_mean, alpha_pr_sd, 1);
      for(i=0; i < nresp; i++){
        for(k=0; k < ndim; k++){
          mle += R::dnorm4(oldz(i, k), z_pr_mean, z_pr_sd, 1);
        }
      }
      for(i = 0; i < nresp; i++){
        for(j =0; j < nresp; j++){
          if(j != i){
            mle += data(i, j) * (oldalpha + olddist(i, j)) - log(1 + exp(oldalpha + olddist(i, j)));
          }
        }
      }
      sample_mle(count) = mle;
      count++;
    } // end burning and thinning and storing
    //Rprintf("Error05-%d\n", t);
    
    if(t % nprint == 0) {
      Rprintf("************* niter %d alpha %f z(0,0) %f likelihood %f ***************\n",t, oldalpha, oldz(0,0), mle);
    }
    
  } // end iteration
  
  // Output
  Rcpp::List output;
  output["alpha"] = samp_alpha;
  output["z"] = samp_z;
  output["map"] = sample_mle;
  output["accept_alpha"] = accept_alpha;
  output["accept_z"] = accept_z;
  
  return(output);
  
} // end function
  
  
  
