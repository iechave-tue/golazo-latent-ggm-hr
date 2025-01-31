soft_shrink = function(z,tau){
  return(sign(z)*pmax(abs(z)-tau,0))
}

soft_shrink_LU = function(Y,lambda,L,U){
  d=nrow(Y)
  return(Pmin(Y-lambda*L,matrix(0,d,d))+Pmax(Y-lambda*U,matrix(0,d,d)))
}


LU_constraints = function(d){
  aux = list()
  constraints = list()
  
  L=matrix(-1,d,d)
  U=matrix(1,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[1]] = aux
  
  L=matrix(-1,d,d)
  U=matrix(Inf,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[2]] = aux
  
  L=matrix(0,d,d)
  U=matrix(Inf,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[3]] = aux
  
  L=matrix(0,d,d)
  U=matrix(1,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[4]] = aux
  
  return(constraints)
}

GGM_constraints = function(d){
  aux = list()
  constraints = list()
  
  L=matrix(-1,d,d)
  U=matrix(1,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[1]] = aux
  
  L[1:(d/2),(d/2+1):d] = L[(d/2+1):d,1:(d/2)] = -Inf
  U[1:(d/2),(d/2+1):d] = U[(d/2+1):d,1:(d/2)] = Inf
  aux$L = L
  aux$U = U
  constraints[[2]] = aux
  
  L=matrix(0,d,d)
  U=matrix(Inf,d,d)
  diag(U)=diag(L)=rep(0,d)
  aux$L = L
  aux$U = U
  constraints[[3]] = aux
  
  L[1:(d/2),(d/2+1):d] = L[(d/2+1):d,1:(d/2)] = -Inf
  U[1:(d/2),(d/2+1):d] = U[(d/2+1):d,1:(d/2)] = Inf
  aux$L = L
  aux$U = U
  constraints[[4]] = aux
  
  return(constraints)
}

objective = function(C, X, xi, S, eigL, v, nu){
  return(sum(sum(C)))
}

lsadmm_parameters = function(d,lambda,gamma=2){
  beta = 1.0/d
  para <- list(
    alpha = 1,
    eta = 4,
    tau = 1.001*(2 + 1)/2,
    TOL = 1e-5,
    tol1 = 1e-6,
    v = lambda,
    nu=lambda*gamma,
    MAXITER = 500,
    continuation = 1,
    num_continuation=10,
    muf = 1e+6
  )
  aux = list()
  aux$beta = beta
  aux$para = para
  return(aux)
}


lsadmm = function(S, beta, para, Lpen, Upen){
  n = nrow(S)
  alpha = para$alpha
  TOL = para$TOL
  tol1 = para$tol1
  r_B = 1.001*beta
  r_A = r_B
  tau = para$tau
  
  v = para$v
  nu = para$nu
  
  MAXITER = para$MAXITER
  
  EY = diag(n)
  M = diag(n)
  A = diag(n)
  B = matrix(0, nrow = n, ncol = n)
  lambda = matrix(0, nrow = n, ncol = n)
  
  for (k in 2:MAXITER){
    M_old = M
    A_old = A
    B_old = B
    
    ev = eigen(S + beta*(B-A) - lambda)
    Q = ev$vectors
    xi = (-ev$values + sqrt(ev$values^2 + 4*beta))/(2*beta)
    M = Q%*%diag(xi)%*%t(Q)
    M = (M+t(M))/2
    
    lambda_mid = lambda - alpha*beta*(M-A+B)
    
    A_tem = A - 1/(tau*r_A)*(lambda_mid)
    v_tem = v/(tau*r_A)
    A = soft_shrink_LU(A_tem, v_tem, Lpen, Upen)
    A = (A+t(A))/2
    
    B_tem = B +1/(tau*r_B)*(lambda_mid - nu*EY)
    ev = eigen(B_tem)
    eigB = pmax(ev$values,rep(0,n))
    B = ev$vectors%*%diag(eigB)%*%t(ev$vectors)
    B = (B+t(B))/2
    
    BB_OLD = B - B_old
    AA_OLD = A - A_old
    lambda = lambda_mid - beta*(BB_OLD - AA_OLD)
    lambda = (lambda + t(lambda))/2
    
    equ = M - A + B
    equnorm = norm(equ,"F")
    
    A_error=norm(AA_OLD,"F")/(1+norm(A_old,"F"))
    M_error=norm(M- M_old,"F")/(1+norm(M_old,"F"))
    B_error=norm(BB_OLD,"F")/(1+norm(B_old,"F"))
    error = max(A_error,max(B_error,M_error))
    if ((equnorm <= tol1) && error <= TOL){
      break
    }
  }
  return(list(M = M, A = A, B = B))
}


twocycles = function(p){
  p=p/2
  aux = diag(5,p,p)
  aux[1,p] = aux[p,1] = -2
  for (k in 1:(p-1)){
    aux[k,k+1] = aux[k+1,k] = -2
  }
  K = matrix(0,2*p+1,2*p+1)
  K[1:p,1:p] = aux
  K[(p+1):(2*p),(p+1):(2*p)] = aux
  for(i in 1:(2*p)){
    K[i,2*p+1] = K[2*p+1,i] = 5/p
  }
  K[2*p+1,2*p+1]=5
  return(K)
}
