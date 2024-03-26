library(Rcpp)
library(RcppArmadillo)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/Null.cpp")
Null(U)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/lthmom.cpp")
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/functions/eigval.cpp")
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/lcr_my.cpp")
ln2_test = n2mom(1,1,3)
lthmom(theta=c(1/2,1/2,1/2,1/2), lambda=c(0.4,0.3,0.2,0.1), lmax=20)
lcr(c(0.4,0.3,0.2,0.1), 20)

lambda=c(0.33,0.27,0.25,0.15)
lthmom(theta=c(1/2,1/2,1/2,1/2), lambda, 20)
exp(lcr(lambda, 20))

sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/rgamma_test.cpp")
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/rbinon_test.cpp")

binom_samp = c()
for(i in 1:1000){
  binom_samp[i] = rbinom_test(1, 0.95)
}

par(mfrow=c(1,3))
hist(binom_samp)
hist(gamma_samp2)
hist(gamma_samp_r)

samp1 = c(); samp2 = c()
for(i in 1:1000){
  samp1[i] = rxl(10, 1, 500, 1)
}

for(i in 1:1000){
  samp2[i] = rxl(10, 1, 500, 1)
}

samp11 = rvmf(c(0,0,1))
samp55 = rvmf(c(0,0,5))

Null<-function (M){
  tmp <- qr(M)
  set <- if (tmp$rank == 0) 
    1:ncol(M)
  else -(1:tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
j=1
Nu<-Null(U[,-j])
Nv<-Null(V[,-j])
E_j<- data$Y-U[,-j]%*%D[-j,-j]%*%t(V[,-j])


E_tilde<- t(Nu)%*%E_j%*%Nv
E_tilde_n<-dim(E_tilde)[2]
E_tilde_m<-dim(E_tilde)[1] # 51
svd_Etilde<-svd(E_tilde)
z<-svd_Etilde$d^2
E02<-sum(z)
theta<-z/sum(z)
norm(E_tilde)^4

E02*phi^2
eigval(E_tilde)

nonzero_d = 0
newd = diag(D)
for(i in 1:46){
  if(i!=j){
    if(newd[i] != 0){
      nonzero_d = nonzero_d + 1;
    }
  }
}


# lmax 구하기..
a = svd(E_tilda)$d[1]^2 * phi^2 * sum(svd(E_tilda)$d^2) * max(svd(E_tilda)$d^2/sum(svd(E_tilda)$d^2))


# eigen value test
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/eigen_test.cpp")
eigen_test(X)


### functions 
# Moment 구하는 함수 test
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/functions/ln2moment.cpp")
ln2moment(1,2,10)  # log(1) log(4) (성공)

sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/lcr.cpp")
lcr(c(0.4,0.3,0.2,0.1),20)
lcr_my()


# 3,4차 방정식 test (성공)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/functions/roots_cubic.cpp")
roots_cubic(-3,2,0)  # roots : 0,1,2
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/functions/roots_quartic.cpp")
roots_quartic(-6,11,-6,0)  # roots : 0,1,2,3


# d 뽑는 rxl 함수 test (성공)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/functions/rxl.cpp")
rxl(mu = 1, sigma=1, 50, nu=1)

sampd = c()
sampd_my = c()
for(i in 1:1000){
  sampd[i] = rxl(mu = 1, sigma=1, 50, nu=1)
}

for(i in 1:1000){
  sampd_my[i] = rxl(mu = 1, sigma=1, 50, nu=1)
}
par(mfrow=c(1,2))
hist(sampd, breaks=20)
hist(sampd_my, breaks=20)
table


mu=1
sigma=1
l=50
nu=1


# Sample 함수 test(성공)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/sample_my.cpp")
sample_samp = sample_N(50, 1000, 1:50)

hist(sample_samp)
par(mfrow=c(1,1))


# (u, v) joint에서 뽑는 ruv 함수(아마도 성공)
set.seed(529)
A = matrix(sample(1:20, 20, 1), nrow = 5)
sourceCpp("/Users/seoyoung/Desktop/Meeting/bayesian_svd/bsvd_varrank/ruv_A.cpp")
ruv_A(A)


