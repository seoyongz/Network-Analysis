library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(ggrepel)

setwd("/Users/seoyoung/Desktop/yonsei_study/papers/LSIRM")
sourceCpp("LSIRM_my.cpp")

## functions
z_postmean = function(Z, max.address){
  ##  Procrustes transform with MAP and estimate with posterior mean
  nmcmc = dim(Z)[1]
  nresp = dim(Z)[2]
  ndim = dim(Z)[3]
  
  z_star = Z[max.address,,]
  z_proc = array(0,dim=c(nmcmc,nresp,ndim))
  
  for(iter in 1:nmcmc){
    z_iter = Z[iter,,]
    
    if(iter != max.address) {
      for( i in 1:ndim) { 
        z_iter[, i] = z_iter[, i] - mean(z_iter[, i]) + mean(z_star[,i]) 
      } #translation
      
      A = t(z_iter)%*%(  z_star%*%t(z_star)  )%*%z_iter
      eA = eigen(A, symmetric=T)
      Ahalf = eA$vec[,1:ndim]%*%diag(sqrt(eA$val[1:ndim]))%*%t(eA$vec[,1:ndim])
      
      z_proc[iter,,] = t(t(z_star)%*%z_iter%*%solve(Ahalf)%*%t(z_iter))
    }
    else {
      z_proc[iter,,] = z_iter
    }
  }
  
  z_est = apply(z_proc, c(2,3), mean)
  return(z_est)
}

## data
elem_data = read.csv("/Users/seoyoung/Desktop/Meeting/My_work/HLSIRM_Incheon_fitting/data/elem_data1207.csv")
head(elem_data)
data = as.matrix(elem_data[, -1])
head(data)

## Parameters Setting
ndim = 2
niter = 10000
nburn = 3000
nthin = 5
nprint = 50
jump_alpha = 0.8
jump_beta = 0.3
jump_gamma = 0.01
jump_a = 3.0
jump_b = 1.5

gamma_pr_mean = 1.0
gamma_pr_sd = 10
sigma_pr_a = 0.001
sigma_pr_b = 0.001
beta_pr_sd = 100
a_pr_sd = 10
b_pr_sd = 10

## Fitting
LSIRM_fit2 = LSIRM_my(data=data, niter=niter, nburn=nburn, nthin=nthin, ndim=ndim, nprint=nprint,
                      jump_alpha=jump_alpha, jump_beta=jump_beta, jump_gamma=jump_gamma, jump_a=jump_a, jump_b=jump_b,
                      gamma_pr_mean=gamma_pr_mean, gamma_pr_sd=gamma_pr_sd,
                      a_pr_sd=a_pr_sd, b_pr_sd=b_pr_sd,
                      sigma_pr_a=sigma_pr_a, sigma_pr_b=sigma_pr_b, beta_pr_sd=beta_pr_sd)
  


## Convergence Check
output = LSIRM_fit2
output$accept_alpha
output$accept_beta
output$accept_a
output$accept_b
output$accept_gamma
ts.plot(output$alpha[,1])
ts.plot(output$beta[,1])
ts.plot(output$a[,1,1])
ts.plot(output$b[,1,1])
ts.plot(output$gamma)
ts.plot(output$map)

## procrustes matching
max.address = which.max(output$map)
a_est = as.data.frame(z_postmean(output$a, max.address))
b_est = as.data.frame(z_postmean(output$b, max.address))
colnames(a_est) = c("z1", "z2")
colnames(b_est) = c("z1", "z2")

## Plot the latent positions
data_labels = sub("^.{3}([^0-9]*)", "\\1", colnames(elem_data[,-1]))

## latent position map (distance model)
ggplot(b_est) +
  geom_point(aes(x = z1, y=z2),) +
  geom_text_repel(aes(x = z1, y=z2, label = 1:82), fontface=2, size=5)+
  theme_minimal()
ggsave("/Users/seoyoung/Desktop/yonsei_study/papers/LSIRM/Incheon_elem_item_map.jpg", width=8, height=8)


## latent position map (projection model)
ggplot(z_est) +
  geom_vline(aes(xintercept=0), color="red", size=0.5, linetype="dashed")+
  geom_hline(aes(yintercept=0), color="red", size=0.5, linetype="dashed")+
  geom_segment(aes(x=rep(0, length(z1)), y=rep(0, length(z2)), xend=z1, yend=z2), 
               arrow = arrow(type = "open", length = unit(0.05, "inches")), size = 0.2) +
  geom_text_repel(aes(x = z1, y=z2, label = data_labels), fontface=2, size=5)+
  theme_minimal()+
  labs(x="z1", y="z2")
ggsave("/Users/seoyoung/Desktop/yonsei_study/papers/LSIRM/Incheon_elem_item_map.jpg", width=8, height=8)




