// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace arma;


arma::mat l2_dist(arma::mat a, arma::mat b){
  
  const int nresp = a.n_rows;
  const int nitem = b.n_rows;
  const int ndim = a.n_cols;
  int i, j, k;
  double dist_tmp;
  arma::mat dist_mat(nresp, nitem);
  
  for(i=0; i < nresp; i++){
    for(j=0; j < nitem; j++){
      dist_tmp = 0.0;
      for(k=0; k < ndim; k++){
        dist_tmp += std::pow(a(i, k) - b(j, k), 2.0);
      }
      dist_tmp = sqrt(dist_tmp);
      dist_mat(i, j) = -dist_tmp;
    }
  }
  
  return(dist_mat);
}

// [[Rcpp::export]]
Rcpp::List LSIRM_my(arma::mat data, const int niter, const int nburn, const int nthin, const int ndim, const int nprint,
                     const double jump_alpha, const double jump_beta, const double jump_gamma, const double jump_a, const double jump_b,
                     const double gamma_pr_mean, const double gamma_pr_sd,
                     const double a_pr_sd, const double b_pr_sd,
                     const double sigma_pr_a, const double sigma_pr_b, const double beta_pr_sd){
  
  // Define the type and size of variables
  int t, i, j, k, accept, count;
  accept = count = 0;
  double alpha_pr_mean = 0.0, beta_pr_mean = 0.0, a_pr_mean = 0.0, b_pr_mean = 0.0;
  double old_like, new_like;
  double u, ratio, mle;
  double sigma_post_a, sigma_post_b;
  
  const int nresp = data.n_rows;
  const int nitem = data.n_cols;
  const int nmcmc = (niter-nburn)/nthin;
  
  arma::dmat samp_alpha(nmcmc, nresp, fill::zeros);
  arma::dmat samp_beta(nmcmc, nitem, fill::zeros);
  arma::dvec samp_gamma(nmcmc, fill::zeros);
  arma::dcube samp_a(nmcmc, nresp, ndim, fill::zeros);
  arma::dcube samp_b(nmcmc, nitem, ndim, fill::zeros);
  arma::dvec samp_sigma2(nmcmc, fill::zeros);
  arma::vec samp_mle(nmcmc, fill::zeros);
  
  double sigma2 = 1.0;
  
  arma::dvec oldalpha(nresp, fill::randu);
  oldalpha = oldalpha * 4.0 - 2.0;
  arma::dvec newalpha = oldalpha;
  
  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;
  
  double newgamma = 1.0, oldgamma = 1.0;
 
  arma::dmat olda(nresp, ndim, fill::randu);
  olda = olda * 2.0 - 1.0;
  arma::dmat newa = olda;
  
  arma::dmat oldb(nitem, ndim, fill::randu);
  oldb = oldb * 2.0 - 1.0;
  arma::dmat newb = oldb;
  
  arma::dvec accept_alpha(nresp, fill::zeros);
  arma::dvec accept_beta(nitem, fill::zeros);
  double accept_gamma; 
  accept_gamma = 0.0;
  arma::dvec accept_a(nresp, fill::zeros);
  arma::dvec accept_b(nitem, fill::zeros);
  
  arma::dmat olddist = l2_dist(olda, oldb);
  arma::dmat newdist = olddist;
  
  
  // iteration for niter
  // Rprintf("Error00-%d\n", 0);
  for(t=0; t < niter; t++){
    
    // Update beta : items main effects
     for(j=0; j < nitem; j++){
       newbeta(j) = R::rnorm(oldbeta(j), jump_beta);
    
       // Calculate Likelihood
       new_like = old_like = 0.0;
       for(i=0; i < nresp; i++){
         new_like += data(i, j) * (oldalpha(i) + newbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + newbeta(j) + oldgamma*olddist(i, j)));
         old_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
       }
       // Rprintf("Error01-%d\n", t);
    
       // Calculate acceptance Ratio
       ratio = new_like - old_like;
       ratio += R::dnorm4(newbeta(j), beta_pr_mean, beta_pr_sd, 1) - R::dnorm4(oldbeta(j), beta_pr_mean, beta_pr_sd, 1);
    
       // Accept or Reject
       if(ratio > 0.0) {
         accept=1;
       }
       else{
         u = std::log(R::runif(0, 1));
         if(u < ratio) {
           accept = 1;
         }
         else {
           accept = 0;
         }
       }
       if(accept==1){
         oldbeta(j) = newbeta(j);
         accept_beta(j) += 1.0 / (niter * 1.0);
       }
     }
     // Rprintf("Error02-%d\n", t);

    
    
    // Update alpha : respondent main effects
    for(i=0; i < nresp; i++){
      newalpha(i) = R::rnorm(oldalpha(i), jump_alpha);

      // Calculate Likelihood
      new_like = old_like = 0.0;
      for(j=0; j < nitem; j++){
        new_like += data(i, j) * (newalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(newalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
        old_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
      }
      // Calculate acceptance Ratio
      ratio = new_like - old_like;
      ratio += R::dnorm4(newalpha(i), alpha_pr_mean, sqrt(sigma2), 1) - R::dnorm4(oldalpha(i), alpha_pr_mean, sqrt(sigma2), 1);

      // Accept or Reject
      if(ratio > 0.0) accept=1;
      else{
        u = std::log(R::runif(0, 1));
        if(u < ratio) accept = 1;
        else accept = 0;
      }

      if(accept==1){
        oldalpha(i) = newalpha(i);
        accept_alpha(i) += 1.0 / (niter * 1.0);
      }
    }
    // Rprintf("Error03-%d\n", t);


    

    // Update gamma with Metropolis algorithm (not symmetric proposal distribution)
    newgamma = R::rlnorm(std::log(oldgamma), jump_gamma);

    // Calculate Likelihood
    new_like = old_like = 0.0;
    for(i=0; i < nresp; i++){
      for(j=0; j < nitem; j++){
        new_like += data(i, j) * (oldalpha(i) + oldbeta(j) + newgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + newgamma*olddist(i, j)));
        old_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
      }
    }

    // Calculate acceptance Ratio
    ratio = new_like - old_like;
    ratio += R::dlnorm(oldgamma, std::log(newgamma), jump_gamma, 1) - R::dlnorm(newgamma, std::log(oldgamma), jump_gamma, 1) ;
    ratio += R::dlnorm(newgamma, gamma_pr_mean, gamma_pr_sd, 1) - R::dlnorm(oldgamma, gamma_pr_mean, gamma_pr_sd, 1);
    
                           
    if(ratio > 0.0) {
      accept=1;
    }
    else{
      u = std::log(R::runif(0, 1));
      if(u < ratio) {
        accept=1;
      }
      else {
        accept=0;
      }
    }

    if(accept==1){
      oldgamma = newgamma;
      accept_gamma += 1.0 / (niter * 1.0);
    }
    // Rprintf("Error04-%d\n", t);
    
    

    // Update latent position of respondent : ai
    for(i = 0; i < nresp; i++){
      for(k = 0; k < ndim; k++) {
        // Propose new a_i
        newa(i, k) = R::rnorm(olda(i, k), jump_a);
      }
    }
    newdist = l2_dist(newa, oldb);
    
    
    for(i = 0; i < nresp; i++){
      // Calculate likelihood for each i
      old_like = new_like = 0.0;
      for(j=0; j<nitem; j++){
        new_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*newdist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*newdist(i, j)));
        old_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
      }
      ratio = new_like - old_like;
      for(k=0; k<ndim; k++){
        ratio += R::dnorm4(newa(i, k), a_pr_mean, a_pr_sd, 1) - R::dnorm4(olda(i, k), a_pr_mean, a_pr_sd, 1);
      }
      
      if(ratio > 0.0) {
        accept = 1;
      }
      else{
        u = R::runif(0,1);
        if(std::log(u) < ratio) {
          accept = 1;
        }
        else {
          accept = 0;
        }
      }
      
      if(accept == 1){
        for(k = 0; k < ndim; k++) {
          olda(i, k) = newa(i, k);
        }
        accept_a(i) += 1.0 / (niter * 1.0);
        for(j = 0; j < nitem; j++){
          olddist(i, j) = newdist(i, j);
        }
      }
    }
    // Rprintf("Error05-%d\n", t);


    
    // bj update
    for(j = 0; j < nitem; j++){
      for(k = 0; k < ndim; k++) {
        // Propose new b_j
        newb(j, k) = R::rnorm(oldb(j, k), jump_b);
      }
    }
    newdist = l2_dist(olda, newb);
      
    for(j=0; j<nitem; j++){
      //calculate likelihood
      old_like = new_like = 0.0;
      for(i = 0; i < nresp; i++){
        new_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*newdist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*newdist(i, j)));
        old_like += data(i, j) * (oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + oldgamma*olddist(i, j)));
      }
      ratio = new_like - old_like;
      for(k=0; k<ndim; k++){
        ratio += R::dnorm4(newb(j, k), b_pr_mean, b_pr_sd, 1) - R::dnorm4(oldb(j, k), b_pr_mean, b_pr_sd, 1);
      }
      
      if(ratio > 0.0) {
        accept = 1;
      }
      else{
        u = R::runif(0,1);
        if(std::log(u) < ratio) {
          accept = 1;
        }
        else {
          accept = 0;
        }
      }
      
      if(accept == 1){
        for(k = 0; k < ndim; k++) {
          oldb(j, k) = newb(j, k);
        }
        accept_b(j) += 1.0 / (niter * 1.0);
        for(i=0; i<nresp; i++){
          olddist(i, j) = newdist(i, j);
        }
      }
    }
    // Rprintf("Error06-%d\n", t);
    

    // Update sigma with Gibbs sampling
    sigma_post_a = sigma_pr_a  + nresp/2.0;
    sigma_post_b = sigma_pr_b;
    for(i = 0; i < nresp; i++) {
      sigma_post_b += std::pow((oldalpha(i) - alpha_pr_mean), 2.0) / 2.0;
    }
    sigma2 = 1.0/R::rgamma(sigma_post_a, 1.0/sigma_post_b);
    // Rprintf("Error07-%d\n", t);
    
    // Store samples with burning and thinning
    if(t >= nburn && t % nthin == 0) {
      for(i=0; i < nresp; i++) {
        samp_alpha(count, i) = oldalpha(i);
      }
      
      for(j=0; j < nitem; j++) {
        samp_beta(count, j) = oldbeta(j);
      }
      
      samp_gamma(count) = oldgamma;
      
      for(i=0; i < nresp; i++) {
        for(k=0; k < ndim; k++){
          samp_a(count, i, k) = olda(i, k);
        }
      }
      
      for(j=0; j < nitem; j++){
        for(k=0; k < ndim; k++){
          samp_b(count, j, k) = oldb(j, k);
        }
      }
      
      samp_sigma2(count) = sigma2;
      // Rprintf("Error08-%d\n", t);
      
      // Compute log-likelihood 
      mle = 0.0;
      for(i=0; i<nresp; i++){
        mle += R::dnorm4(oldalpha(i), alpha_pr_mean, sqrt(sigma2), 1);
      }
      
      for(j=0; j<nitem; j++){
        mle += R::dnorm4(oldbeta(j), beta_pr_mean, beta_pr_sd, 1);
      }
      
      mle += R::dgamma(1.0/sigma2, sigma_pr_a, 1.0/sigma_pr_b, 1);
      
      mle += R::dlnorm(oldgamma, gamma_pr_mean, gamma_pr_sd, 1);
      
      for(i=0; i < nresp; i++){
        for(k=0; k < ndim; k++){
          mle += R::dnorm4(olda(i, k), a_pr_mean, a_pr_sd, 1);
        }
      }
      
      for(j=0; j < nitem; j++){
        for(k=0; k < ndim; k++){
          mle += R::dnorm4(oldb(j, k), b_pr_mean, b_pr_sd, 1);
        }
      }
      
      for(i = 0; i < nresp; i++){
        for(j =0; j < nitem; j++){
          mle += data(i, j) * (oldalpha(i) + oldbeta(j) + olddist(i, j)) - log(1 + exp(oldalpha(i) + oldbeta(j) + olddist(i, j)));
        }
      }
      samp_mle(count) = mle;
      
      count++;
    }
    
    if(t % nprint == 0) {
      Rprintf("************* niter %d alpha(0) %f a(0,0) %f likelihood %f ***************\n",t, oldalpha(0), olda(0,0), mle);
    }
    

  } // end iteration

  
  // Define output
  Rcpp::List output;    
  output["alpha"] = samp_alpha;
  output["beta"] = samp_beta;
  output["gamma"] = samp_gamma;
  output["a"] = samp_a;
  output["b"] = samp_b;
  output["sigma_alpha"] = samp_sigma2;
  output["accept_alpha"] = accept_alpha;
  output["accept_beta"] = accept_beta;
  output["accept_gamma"] = accept_gamma;
  output["accept_a"] = accept_a;
  output["accept_b"] = accept_b;
  output["map"] = samp_mle;
  
  return(output);
  
} //  end function

  
  


