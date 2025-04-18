library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(latex2exp)
library(glmnet)
library(egg)
library(cowplot)
library(tictoc)
library(here)
library(clusterGeneration)
library(pracma)
library(matrixcalc)
library(CVXR)
library(ggplot2)
library(ggpubr)
library(Rfast)
library(matrixcalc)

soft_shrink_LU = function(Y,lambda,L,U){
  d=nrow(Y)
  return(Pmin(Y-lambda*L,matrix(0,d,d))+Pmax(Y-lambda*U,matrix(0,d,d)))
}

objective = function(C, X, xi, S, eigL, v, nu){
  return(sum(sum(C)))
}

laplsadmm_fit = function(C, lambda1_list, lambda2_list, Lpen, Upen){
  d = nrow(C)
  beta = 1.0/d
  para <- list(
    alpha = 1,
    eta = 4,
    tau = 1.001*(2 + 1)/2, # the 1 is alpha
    TOL = 1e-5,
    tol1 = 1e-6,
    v = 0, #lambda1
    nu=0, #lambda2
    MAXITER = 500,
    continuation = 1,
    num_continuation=10,
    muf = 1e+6
  )
  i=1
  results = list()
  for(lambda1 in lambda1_list){
    for(lambda2 in lambda2_list){
      para$v = lambda1
      para$nu = lambda1 * lambda2
      res = laplsadmm(C, beta, para, Lpen, Upen)
      results[[i]] = res
      i = i + 1
    }
  }
  return(results)
}

laplsadmm = function(C, beta, para, Lpen, Upen){
  n = nrow(C)
  alpha = para$alpha
  TOL = para$TOL
  tol1 = para$tol1
  r_L = 1.001*beta
  r_S = r_L
  tau = para$tau
  
  v = para$v
  nu = para$nu
  
  MAXITER = para$MAXITER
  
  # TODO STARTING CONDITIONS
  EY = diag(n)
  X = diag(n)
  S = diag(n)
  L = matrix(0, nrow = n, ncol = n)
  P = eigen(diag(n)-matrix(1/n,n,n))$vectors[,-n]
  C = t(P)%*%C%*%P
  lambda = matrix(0, nrow = n, ncol = n)
  
  for (k in 2:MAXITER){
    X_old = X
    S_old = S
    L_old = L
    
    ev = eigen(C + beta*t(P)%*%(L-S)%*%P - t(P)%*%lambda%*%P)
    Q = ev$vectors
    xi = (-ev$values + sqrt(ev$values^2 + 4*beta))/(2*beta)
    Xi = Q%*%diag(xi)%*%t(Q)
    Xi = (Xi+t(Xi))/2
    X = P%*%Xi%*%t(P)
    
    lambda_mid = lambda - alpha*beta*(X-S+L)
    
    S_tem = S - 1/(tau*r_S)*(lambda_mid)
    v_tem = v/(tau*r_S)
    S = soft_shrink_LU(S_tem, v_tem, Lpen, Upen)
    diag(S) = pmax(diag(S),0)
    #S = soft_shrink(S_tem,v_tem)
    # TODO CHECK
    S = (S+t(S))/2
    
    L_tem = L +1/(tau*r_L)*(lambda_mid - nu*EY)
    ev = eigen(L_tem)
    eigL = pmax(ev$values,rep(0,n))
    L = ev$vectors%*%diag(eigL)%*%t(ev$vectors)
    # TODO CHECK
    L = (L+t(L))/2
    
    LL_OLD = L - L_old
    SS_OLD = S - S_old
    lambda = lambda_mid - beta*(LL_OLD - SS_OLD)
    # TODO CHECK
    lambda = (lambda + t(lambda))/2
    
    # TODO ADD CHECKS
    equ = X - S + L
    equnorm = norm(equ,"F")
    
    S_error=norm(SS_OLD,"F")/(1+norm(S_old,"F"))
    X_error=norm(X- X_old,"F")/(1+norm(X_old,"F"))
    L_error=norm(LL_OLD,"F")/(1+norm(L_old,"F"))
    error = max(S_error,max(L_error,X_error))
    if ((equnorm <= tol1) && error <= TOL){
      break
    }
  }
  return(list(X = X, S = S, L = L))
}

generate_latent_model_cycle <- function(p, h) {
  W <- matrix(1, p + h, p + h)
  S <- matrix(0, p, p)
  for (i in 1:(p - 1)) {
    S[i, i + 1] <- 1
  }
  S[1, p] <- 1
  S <- S + t(S)
  diag(S) <- 0
  
  W[1:p, 1:p] <- S
  diag(W) <- 0
  L <- matrix(runif((p + h)^2, 2, 2), nrow = p + h)
  if(h!=0){
    L[(p + 1):(p + h), (p + 1):(p + h)] <- diag(h)
    
    
    z <- 50 / (sqrt(as.integer(p / h))) * matrix(runif(as.integer(p / h), 1, 1.5), nrow = as.integer(p / h))
    
    for (i in 1:h) {
      L[1:p, (p + i):(p + i)] <- 0
      L[(p + i):(p + i), 1:p] <- 0
      L[seq(i, p, h), (p + i):(p + i)] <- z
      L[(p + i):(p + i), seq(i, p, h)] <- z
    }
  }
  W <- W * L
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  Theta <- diag(rowSums(W)) - W
  G <- Theta2Gamma(Theta)
  return(list(Gamma = G, graph = Gamma2graph(G), Lst = NULL))
}


twocycles_lap = function(p){
  p=p/2
  aux = matrix(0,p,p)
  aux[1,p] = aux[p,1] = -2
  for (k in 1:(p-1)){
    aux[k,k+1] = aux[k+1,k] = -2
  }
  K = matrix(0,2*p+1,2*p+1)
  K[1:p,1:p] = aux
  K[(p+1):(2*p),(p+1):(2*p)] = aux
  for(i in 1:(2*p)){
    if(i%%2==0) {K[i,2*p+1] = K[2*p+1,i] = -0.4}
    else{ K[i,2*p+1] = K[2*p+1,i] = 0.2 }
  }
  for(i in 1:(2*p+1)){
    aux = sum(K[i,])
    K[i,i] = -aux
  }
  return(K)
}
