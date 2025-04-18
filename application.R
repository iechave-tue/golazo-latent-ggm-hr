library(Rfast)
library(matrixcalc)
library(MASS)
library(golazo)
library(ggplot2)
#devtools::install_github('stefano-meschiari/latex2exp') needed for some symbols
library(latex2exp)
source("functions.R")

# Two cycles simulation


listp = c(50)
listh = c(1)
listn = c(100)
listlambda = seq(1e-8,1,length.out = 50)

gamma = 0.5

trials = 20

set.seed(0)

loglikelihood = array(data = 0,dim = c(length(listp),length(listh),length(listn),length(listlambda),6))
edges = array(data = 0,dim = c(length(listp),length(listh),length(listn),length(listlambda),6))
rank = array(data = 0,dim = c(length(listp),length(listh),length(listn),length(listlambda),6))

estimators = array(data = 0,dim = c(length(listp),length(listh),length(listn),length(listlambda),6,trials,3,listp[1],listp[1]))

for (i in 1:length(listp)){
  for (j in 1:length(listh)){
    p = listp[i]
    h = listh[j]
    K = twocycles(p)
    for (k in 1:length(listn)){
      n = listn[k]
      constraints = GGM_constraints(p)
      for(aux in 1:trials){
        print(aux)
        data = mvrnorm(n=n,mu=rep(0,p+h), Sigma = solve(K))[,1:p]
        test = mvrnorm(n=n,mu=rep(0,p+h), Sigma = solve(K))[,1:p]
        for (l in 1:length(constraints)){
          constraint = constraints[[l]]
          L = constraint$L
          U = constraint$U
          for (m in 1:length(listlambda)){
            lambda = listlambda[[m]]
            d <- p
            beta = 1.0/d
            para <- list(
              alpha = 1,
              eta = 4,
              tau = 1.001*(2 + 1)/2,
              TOL = 1e-5,
              tol1 = 1e-6,
              v = gamma*lambda,
              nu = lambda,
              MAXITER = 500,
              continuation = 1,
              num_continuation=10,
              muf = 1e+6
            )
            S = cov(data)
            res <- lsadmm(S,beta,para,Lpen=L,Upen=U)
            true = K[1:p,1:p]
            loglikelihood[i,j,k,m,l] = loglikelihood[i,j,k,m,l] + (sum(log(eigen(res$M)$values))-sum(diag((res$M)%*%cov(test))))/trials
            edges[i,j,k,m,l] = edges[i,j,k,m,l] + ((sum(abs(res$A)>10e-5) - p)/2)/trials
            rank[i,j,k,m,l] = rank[i,j,k,m,l] + sum(eigen(res$B)$values>10e-5)/trials
            
            estimators[i,j,k,m,l,aux,1,,] = res$M
            estimators[i,j,k,m,l,aux,2,,] = res$A
            estimators[i,j,k,m,l,aux,3,,] = res$B
          }
        }
      }
    }
  }
}

df = data.frame(listlambda, loglikelihood[1,1,1,,1],  loglikelihood[1,1,1,,2], loglikelihood[1,1,1,,3], loglikelihood[1,1,1,,4])
colnames(df) = c("lambda","lasso","lassoggm","positive","positiveggm")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso,color="A"),size=0.8)+
  geom_line(aes(y=lassoggm,color="B"),size=0.8)+
  geom_line(aes(y=positive,color="C"),size=0.8)+
  geom_line(aes(y=positiveggm,color="D"),size=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.3))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with partial GGM constraints"),TeX("${MTP}_2$ constraint"),TeX("${MTP}_2$ constraint with partial GGM constraints")),values = c("black","red","blue","green"))

df = data.frame(listlambda, edges[1,1,1,,1],  edges[1,1,1,,2], edges[1,1,1,,3], edges[1,1,1,,4])
colnames(df) = c("lambda","lasso","lassoggm","positive","positiveggm")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso,color="A"),size=0.8)+
  geom_line(aes(y=lassoggm,color="B"),size=0.8)+
  geom_line(aes(y=positive,color="C"),size=0.8)+
  geom_line(aes(y=positiveggm,color="D"),size=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.5))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with partial GGM constraints"),TeX("${MTP}_2$ constraint"),TeX("${MTP}_2$ constraint with partial GGM constraints")),values = c("black","red","blue","green"))

df = data.frame(listlambda, rank[1,1,1,,1],  rank[1,1,1,,2], rank[1,1,1,,3], rank[1,1,1,,4])
colnames(df) = c("lambda","lasso","lassoggm","positive","positiveggm")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso,color="A"),size=0.8)+
  geom_line(aes(y=lassoggm,color="B"),size=0.8)+
  geom_line(aes(y=positive,color="C"),size=0.8)+
  geom_line(aes(y=positiveggm,color="D"),size=0.8)+
  geom_hline(yintercept = 1, color = "black", linewidth = 0.3, linetype = "solid") +
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.5))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with partial GGM constraints"),TeX("${MTP}_2$ constraint"),TeX("${MTP}_2$ constraint with partial GGM constraints")),values = c("black","red","blue","green"))


# Gene data (Rosetta dataset)

gene = read.csv("Rosetta.csv",header = FALSE)
gene = data.matrix(gene)
variance = apply(gene,1,var)
indices = rev(order(variance))

d=25
data = t(gene[indices[1:d],])

n = nrow(data)
folds = 5
nsplit = floor(n/folds)
lambdalist = seq(1e-8,0.4,length.out = 30)
loglikeli = matrix(0,length(lambdalist),folds)
edges = matrix(0,length(lambdalist),folds)

constraints = LU_constraints(d)
results = list()
edge_results = list()

gamma = 0.1

for(j in 1:length(constraints)){
  loglikeli = matrix(0,length(lambdalist),folds)
  for(i in 1:folds){
    print(i)
    lambdatimes = 1
    valInd = c(((nsplit*(i-1))+1):(nsplit*i))
    validation = data[valInd,]
    train = data[setdiff(1:nrow(data),valInd),]
    for(lambda in lambdalist){
      d <- ncol(train)
      beta = 1.0/d
      para <- list(
        alpha = 1,
        eta = 4,
        tau = 1.001*(2 + 1)/2,
        TOL = 1e-5,
        tol1 = 1e-6,
        v = gamma*lambda,
        nu = lambda,
        MAXITER = 500,
        continuation = 1,
        num_continuation=10,
        muf = 1e+6
      )
      L=constraints[[j]]$L
      U=constraints[[j]]$U
      S=cov(train)
      res <- lsadmm(S,beta,para,Lpen=L,Upen=U)
      loglikeli[lambdatimes,i] = sum(log(eigen(res$M)$values))-sum(diag((res$M)%*%cov(validation)))
      edges[lambdatimes,i] = (sum(abs(res$A)>10e-5) - d)/2
      lambdatimes = lambdatimes + 1
    }
    results[[j]] = loglikeli
    edge_results[[j]] = edges
  }
}

loglikeli = lapply(results,rowSums)
for(i in 1:4){
  loglikeli[[i]] = loglikeli[[i]]/folds
}

edges = lapply(edge_results,rowSums)
for(i in 1:4){
  edges[[i]] = edges[[i]]/folds
}

df = data.frame(lambdalist, loglikeli[[1]],  loglikeli[[2]], loglikeli[[3]], loglikeli[[4]])
colnames(df) = c("lambda","lasso","lassopositive","positive","positivelasso")
colors = c("Graphical lasso"="black","Modified positivity"="red","Positivity"="green","Positive lasso"="blue")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso,color="A"),size=0.8)+
  geom_line(aes(y=lassopositive,color="B"),size=0.8)+
  geom_line(aes(y=positive,color="C"),size=0.8)+
  geom_line(aes(y=positivelasso,color="D"),size=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.50))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("Modified $\\ell_1$-penalty"),TeX("${MTP}_2$ constraint"),TeX("Positive lasso")),values = c("black","red","green","blue"))

df = data.frame(lambdalist, edges[[1]],  edges[[2]], edges[[3]], edges[[4]])
colnames(df) = c("lambda","lasso","lassopositive","positive","positivelasso")
colors = c("Graphical lasso"="black","Modified positivity"="red","Positivity"="green","Positive lasso"="blue")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso,color="A"),size=0.8)+
  geom_line(aes(y=lassopositive,color="B"),size=0.8)+
  geom_line(aes(y=positive,color="C"),size=0.8)+
  geom_line(aes(y=positivelasso,color="D"),size=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.70))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("Modified $\\ell_1$-penalty"),TeX("${MTP}_2$ constraint"),TeX("Positive lasso")),values = c("black","red","green","blue"))




