#install.packages("devtools")
#devtools::install_github("sebastian-engelke/graphicalExtremes")

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


source("functions_extreme.R")

IATAs <- getFlightDelayData("IATAs", "tcCluster")

MIN_N_FLIGHTS <- 20
flight_graph <- getFlightGraph(IATAs, minNFlights = MIN_N_FLIGHTS)

# Plot airports + connections (indicating number of flights by thickness)

gg0 <- plotFlights(
  IATAs,
  graph = flight_graph,
  useAirportNFlights = TRUE,
  useConnectionNFlights = FALSE,
  returnGGPlot = TRUE,
  clipMap = 1.3,
  xyRatio = 1
)

gg0

# Check whether all dates from the train-test-split are available
# (due to size restrictions, the CRAN version of this package does not contain all dates)

allDatesAvailable <- tryCatch(
  {
    getFlightDelayData("dates", dateFilter = c("tcAll"))
    TRUE
  },
  error = function(...) FALSE
)
cat("All dates available:", allDatesAvailable, "\n")

# Full data set 3603 x 29 (if allDatesAvailable == TRUE)
mat <- rbind(drop(getFlightDelayData("delays", "tcCluster", "tcTrain")), drop(getFlightDelayData("delays", "tcCluster", "tcTest")))

mixed_test = 0
lasso_test = 0
positivity_test = 0

lasso_edges = 0
mixed_edges = 0
positivity_edges = 0

lasso_rank = 0
mixed_rank = 0
positivity_rank = 0


# Probability threshold to define extremes
n <- nrow(mat)
p <- 0.95 #1-n^0.65/n #0.85 
k <- (1-p) * n
d <- ncol(mat)

# Regularization parameter to enforce sparsity
lambda1_range <- seq(0.0000000001, .1, by = 0.01) # seq(0, 0.1, by = 0.002)
ll <- length(lambda1_range)
# Regularization parameter to enforce low-rank
lambda2_range <- c(4) # 2 # seq(0, 0.3, by = 0.08) #seq(0, 0.3, by = 0.02)
# for eglearn
rholist <- lambda1_range



## cross-validation for model assessment
set.seed(4134999)



splitInd <- floor(nrow(mat) * 1 / 5)

for (cvsplit in 1:5) {
  
  print(cvsplit)
  
  # Validation set for this iteration
  ValInd <- c((splitInd * (cvsplit - 1) + 1):(splitInd * cvsplit))
  matVal <- mat[ValInd, ]
  # Training set for this iteration
  matEst <- mat[setdiff(1:nrow(mat), ValInd), ]
  
  # Normalize data to multivariate Pareto scale
  train_data <- data2mpareto(data = matEst, p = p)
  test_data <- data2mpareto(data = matVal, p = p)
  
  # empirical variogram on training and test sets
  Gamma_train <- emp_vario(data = train_data)
  Gamma_test <- emp_vario(data = test_data)
  
  L = matrix(-1,d,d)
  U = matrix(1,d,d)
  diag(L) = diag(U) = rep(0,d)
  fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                             Lpen = L, Upen = U)
  lasso_test <- lasso_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
  })
  lasso_edges <- lasso_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
  })
  lasso_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
  })
  
  L = matrix(-1,d,d)
  U = matrix(Inf,d,d)
  diag(L) = diag(U) = rep(0,d)
  fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                             Lpen = L, Upen = U)
  mixed_test <- mixed_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
  })
  mixed_edges <- mixed_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
  })
  mixed_rank <- mixed_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
  })
  
  L = matrix(0,d,d)
  U = matrix(Inf,d,d)
  diag(L) = diag(U) = rep(0,d)
  fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                             Lpen = L, Upen = U)
  positivity_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
  })
  positivity_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
  })
  positivity_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
  })
}


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_test,  mixed_test, positivity_test)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_test/5,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_test/5,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_test/5,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.4,.2))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_edges,  mixed_edges, positivity_edges)
colnames(df) = c("lambda","lasso","mixed","positivity")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_edges/5,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_edges/5,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_edges/5,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.7))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))

df = data.frame(lambda1_range*lambda2_range[[1]], lasso_rank,  mixed_rank, positivity_rank)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_rank/5,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_rank/5,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_rank/5,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.6))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))

## Fiting the models on the whole data set
set.seed(0)

lambda2_range <- c(4)
p=0.95
tol = 0.00001
d=ncol(mat)
train_data <- data2mpareto(data = mat, p = p)
Gamma_train <- emp_vario(data = train_data)

lambda1_range = c(0.03)

L = matrix(-1,d,d)
U = matrix(1,d,d)
diag(L) = diag(U) = rep(0,d)
fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                           Lpen = L, Upen = U)
A <- 1 * (abs(fit_golazo[[1]]$S) > tol)
graph = igraph::graph_from_adjacency_matrix(A, mode = "undirected", 
                                            diag = FALSE)
gg <- plotFlights(
  IATAs,
  graph = graph,
  xyRatio = 1,
  clipMap = 1.3,
  returnGGPlot = TRUE,
  useAirportNFlights = TRUE,
)
gg
graph_lasso = as_adjacency_matrix(graph, sparse=FALSE)


lambda1_range = c(0.03)


L = matrix(-1,d,d)
U = matrix(Inf,d,d)
diag(L) = diag(U) = rep(0,d)
fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                           Lpen = L, Upen = U)
A <- 1 * (abs(fit_golazo[[1]]$S) > tol)
graph = igraph::graph_from_adjacency_matrix(A, mode = "undirected", 
                                            diag = FALSE)
gg <- plotFlights(
  IATAs,
  graph = graph,
  xyRatio = 1,
  clipMap = 1.3,
  returnGGPlot = TRUE,
  useAirportNFlights = TRUE,
)
gg

graph_mixed = as_adjacency_matrix(graph, sparse=FALSE)

lambda1_range = c(0.07)

L = matrix(0,d,d)
U = matrix(Inf,d,d)
diag(L) = diag(U) = rep(0,d)
fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                           Lpen = L, Upen = U)
A <- 1 * (abs(fit_golazo[[1]]$S) > tol)
graph = igraph::graph_from_adjacency_matrix(A, mode = "undirected", 
                                            diag = FALSE)
gg <- plotFlights(
  IATAs,
  graph = graph,
  xyRatio = 1,
  clipMap = 1.3,
  returnGGPlot = TRUE,
  useAirportNFlights = TRUE,
)
gg

graph_positivity = as_adjacency_matrix(graph, sparse=FALSE)

hamming_distance = function(adj1, adj2){
  res = sum(adj1 != adj2)/2
  return(res)
}

graphs = list(as_adjacency_matrix(flight_graph,sparse=FALSE),graph_lasso,graph_mixed,graph_positivity)

hamming = matrix(-1,4,4)
for(i in 1:4){
  for(j in (i:4)){
    hamming[i,j] = hamming [j,i] = hamming_distance(graphs[[i]],graphs[[j]])
  }
}

hamming




### Simulation

set.seed(0)

d = 30

mixed_test = 0
lasso_test = 0
positivity_test = 0
posgraphical_test = 0

lasso_edges = 0
mixed_edges = 0
positivity_edges = 0
posgraphical_edges = 0

lasso_rank = 0
mixed_rank = 0
positivity_rank = 0
posgraphical_rank = 0

#lambda1_range = seq(0.0001, .02, length.out = 20)
lambda1_range = seq(0.0001,0.4,length.out = 20)
ll = length(lambda1_range)
lambda2_range = c(4) # 4
max_iter = 10


for (n in c(100)){
  print(n)
  for (h in c(1)) {
  
    print(h)
    latent_model <- twocycles_lap(d)
    G = Theta2Gamma(latent_model)[1:d,1:d]
    
    for (num_iter in 1:max_iter){
      print(num_iter)
      train_data = rmpareto(n = n, d = d, model = "HR", par = G)
      test_data = rmpareto(n = n, d = d, model = "HR", par = G)
      # empirical variogram on training and test sets
      Gamma_train <- emp_vario(data = train_data)
      Gamma_test <- emp_vario(data = test_data)
      
      
      L = matrix(-1,d,d)
      U = matrix(1,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      lasso_test <- lasso_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      lasso_edges <- lasso_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-10)/2 - d/2
      })/max_iter
      lasso_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(-1,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      mixed_test <- mixed_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      mixed_edges <- mixed_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-10)/2 - d/2
      })/max_iter
      mixed_rank <- mixed_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(0,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      positivity_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      positivity_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-10)/2 - d/2
      })/max_iter
      positivity_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L=matrix(0,d,d)
      U=matrix(Inf,d,d)
      L[1:(d/2),(d/2+1):d] = L[(d/2+1):d,1:(d/2)] = -Inf
      U[1:(d/2),(d/2+1):d] = U[(d/2+1):d,1:(d/2)] = Inf
      diag(U)=diag(L)=rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      posgraphical_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      posgraphical_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-10)/2 - d/2
      })/max_iter
      posgraphical_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
    }
    
    
  }
}


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_test,  mixed_test, positivity_test, posgraphical_test)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_test,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_test,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_test,color="C"),linewidth=0.8)+
  geom_line(aes(y=posgraphical_test,color="D"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.3))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint"),TeX("${EMTP}_2$ with graphical model constraints")),values = c("black","red","blue","green"))


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_edges,  mixed_edges, positivity_edges, posgraphical_edges)
colnames(df) = c("lambda","lasso","mixed","positivity")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_edges,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_edges,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_edges,color="C"),linewidth=0.8)+
  geom_line(aes(y=posgraphical_edges,color="D"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.7))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint"),TeX("${EMTP}_2$ with graphical model constraints")),values = c("black","red","blue","green"))

df = data.frame(lambda1_range*lambda2_range[[1]], lasso_rank,  mixed_rank, positivity_rank, posgraphical_rank)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_rank,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_rank,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_rank,color="C"),linewidth=0.8)+
  geom_line(aes(y=posgraphical_rank,color="D"),linewidth=0.8)+
  geom_hline(yintercept = 1, color = "black", linewidth = 0.3, linetype = "solid") +
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.6))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint"),TeX("${EMTP}_2$ with graphical model constraints")),values = c("black","red","blue","green"))




set.seed(4134999)


d = 30

mixed_test = 0
lasso_test = 0
positivity_test = 0

lasso_edges = 0
mixed_edges = 0
positivity_edges = 0

lasso_rank = 0
mixed_rank = 0
positivity_rank = 0

#lambda1_range = seq(0.0001, .02, length.out = 20)
lambda1_range = seq(0.0001,0.04,length.out = 20)
ll = length(lambda1_range)
lambda2_range = c(4) # 4
max_iter = 10


for (n in c(100)){
  print(n)
  for (h in c(3)) {
    
    print(h)
    latent_model <- generate_latent_model_cycle(p = d, h)
    G = latent_model$Gamma[1:d,1:d]
    
    for (num_iter in 1:max_iter){
      print(num_iter)
      train_data = rmpareto(n = n, d = d, model = "HR", par = G)
      test_data = rmpareto(n = n, d = d, model = "HR", par = G)
      # empirical variogram on training and test sets
      Gamma_train <- emp_vario(data = train_data)
      Gamma_test <- emp_vario(data = test_data)
      
      
      L = matrix(-1,d,d)
      U = matrix(1,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      lasso_test <- lasso_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      lasso_edges <- lasso_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      lasso_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(-1,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      mixed_test <- mixed_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      mixed_edges <- mixed_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      mixed_rank <- mixed_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(0,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      positivity_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      positivity_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      positivity_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
    }
    
    
  }
}



df = data.frame(lambda1_range*lambda2_range[[1]], lasso_test,  mixed_test, positivity_test)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_test,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_test,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_test,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.3))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_edges,  mixed_edges, positivity_edges)
colnames(df) = c("lambda","lasso","mixed","positivity")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_edges,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_edges,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_edges,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.7))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))

df = data.frame(lambda1_range*lambda2_range[[1]], lasso_rank,  mixed_rank, positivity_rank)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_rank,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_rank,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_rank,color="C"),linewidth=0.8)+
  geom_hline(yintercept = 3, color = "black", linewidth = 0.3, linetype = "solid") +
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.6))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))



set.seed(4134999)


d = 30


mixed_test = 0
lasso_test = 0
positivity_test = 0

lasso_edges = 0
mixed_edges = 0
positivity_edges = 0

lasso_rank = 0
mixed_rank = 0
positivity_rank = 0

#lambda1_range = seq(0.0001, .02, length.out = 20)
lambda1_range = seq(0.0001,0.04,length.out = 20)
ll = length(lambda1_range)
lambda2_range = c(4) # 4
max_iter = 10


for (n in c(100)){
  print(n)
  for (h in c(5)) { 
    
    print(h)
    latent_model <- generate_latent_model_cycle(p = d, h)
    G = latent_model$Gamma[1:d,1:d]
    
    for (num_iter in 1:max_iter){
      print(num_iter)
      train_data = rmpareto(n = n, d = d, model = "HR", par = G)
      test_data = rmpareto(n = n, d = d, model = "HR", par = G)
      # empirical variogram on training and test sets
      Gamma_train <- emp_vario(data = train_data)
      Gamma_test <- emp_vario(data = test_data)
      
      
      L = matrix(-1,d,d)
      U = matrix(1,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      lasso_test <- lasso_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      lasso_edges <- lasso_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      lasso_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(-1,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      mixed_test <- mixed_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      mixed_edges <- mixed_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      mixed_rank <- mixed_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(0,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      positivity_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      positivity_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      positivity_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
    }
    
    
  }
}



df = data.frame(lambda1_range*lambda2_range[[1]], lasso_test,  mixed_test, positivity_test)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_test,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_test,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_test,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.3))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_edges,  mixed_edges, positivity_edges)
colnames(df) = c("lambda","lasso","mixed","positivity")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_edges,color="Graphical lasso"),linewidth=0.8)+
  geom_line(aes(y=mixed_edges,color="Graphical lasso with positivity"),linewidth=0.8)+
  geom_line(aes(y=positivity_edges,color="Positivity"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.7))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))

df = data.frame(lambda1_range*lambda2_range[[1]], lasso_rank,  mixed_rank, positivity_rank)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_rank,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_rank,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_rank,color="C"),linewidth=0.8)+
  geom_hline(yintercept = 5, color = "black", linewidth = 0.3, linetype = "solid") +
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.6))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))



set.seed(4134999)


d = 30

mixed_test = 0
lasso_test = 0
positivity_test = 0

lasso_edges = 0
mixed_edges = 0
positivity_edges = 0

lasso_rank = 0
mixed_rank = 0
positivity_rank = 0

#lambda1_range = seq(0.0001, .02, length.out = 20)
lambda1_range = seq(0.0001,0.04,length.out = 20)
ll = length(lambda1_range)
lambda2_range = c(4) # 4
max_iter = 10


for (n in c(100)){
  print(n)
  for (h in c(10)) {
    
    print(h)
    latent_model <- generate_latent_model_cycle(p = d, h)
    G = latent_model$Gamma[1:d,1:d]
    
    for (num_iter in 1:max_iter){
      print(num_iter)
      train_data = rmpareto(n = n, d = d, model = "HR", par = G)
      test_data = rmpareto(n = n, d = d, model = "HR", par = G)
      # empirical variogram on training and test sets
      Gamma_train <- emp_vario(data = train_data)
      Gamma_test <- emp_vario(data = test_data)
      
      
      L = matrix(-1,d,d)
      U = matrix(1,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      lasso_test <- lasso_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      lasso_edges <- lasso_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      lasso_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(-1,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      mixed_test <- mixed_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      mixed_edges <- mixed_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      mixed_rank <- mixed_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
      
      L = matrix(0,d,d)
      U = matrix(Inf,d,d)
      diag(L) = diag(U) = rep(0,d)
      fit_golazo = laplsadmm_fit(-Gamma_train/2, lambda1_list = lambda1_range, lambda2_list = lambda2_range,
                                 Lpen = L, Upen = U)
      positivity_test <- positivity_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        loglik_HR(data = test_data, Gamma = Theta2Gamma(fit_golazo[[i]]$X), cens = FALSE)[1]
      })/max_iter
      positivity_edges <- positivity_edges + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(abs(fit_golazo[[i]]$S)>10e-5)/2 - d/2
      })/max_iter
      positivity_rank <- positivity_rank + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
        sum(eigen(fit_golazo[[i]]$L)$values>10e-5)
      })/max_iter
    }
    
    
  }
}



df = data.frame(lambda1_range*lambda2_range[[1]], lasso_test,  mixed_test, positivity_test)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_test,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_test,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_test,color="C"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Log-likelihood",color="Legend")+
  theme(legend.position = c(.7,.3))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))


df = data.frame(lambda1_range*lambda2_range[[1]], lasso_edges,  mixed_edges, positivity_edges)
colnames(df) = c("lambda","lasso","mixed","positivity")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_edges,color="Graphical lasso"),linewidth=0.8)+
  geom_line(aes(y=mixed_edges,color="Graphical lasso with positivity"),linewidth=0.8)+
  geom_line(aes(y=positivity_edges,color="Positivity"),linewidth=0.8)+
  labs(x=TeX("Value of $\\lambda$"),y="Average number of edges",color="Legend")+
  theme(legend.position = c(.7,.7))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))

df = data.frame(lambda1_range*lambda2_range[[1]], lasso_rank,  mixed_rank, positivity_rank)
colnames(df) = c("lambda","lasso","mixed")
ggplot(data=df, aes(x=lambda)) +
  geom_line(aes(y=lasso_rank,color="A"),linewidth=0.8)+
  geom_line(aes(y=mixed_rank,color="B"),linewidth=0.8)+
  geom_line(aes(y=positivity_rank,color="C"),linewidth=0.8)+
  geom_hline(yintercept = 10, color = "black", linewidth = 0.3, linetype = "solid") +
  labs(x=TeX("Value of $\\lambda$"),y="Average estimated rank",color="Legend")+
  theme(legend.position = c(.7,.6))+
  scale_color_manual(labels = c(TeX("Standard $\\ell_1$-penalty"),TeX("$\\ell_1$-penalty with ${EMTP}_2$"),TeX("${EMTP}_2$ constraint")),values = c("black","red","blue"))













