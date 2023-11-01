#loading required packages
library(stats)
library(parallel)
library(pbapply)
library(spGARCH)
library(Matrix)
library(Rsolnp)
library(rlist)
library(future.apply)


#inverse distance matrix as alternative - but not used in the current version of the paper!
# weight_matrix<-function(covariates){
#   covariates=as.matrix(covariates)
#   n=nrow(covariates)
#   weight_matrix=matrix(0,ncol = n,nrow=n)
#   for(i in 1:n){
#     for(j in 1:n){
#       weight_matrix[i,j] = 1/sqrt(sum((covariates[i,] - covariates[j,])^2))
#       
#       if(i==j){
#         weight_matrix[i,j]=0
#       }
#     }
#     
#   }
#   #standardizing
#   
#   #weight_matrix=(1/max(weight_matrix[weight_matrix!=Inf]))*weight_matrix
#   # weight_matrix[weight_matrix==Inf]=1
#   weight_matrix=weight_matrix*array(rep(1/colSums(weight_matrix),n),dim = c(n,n))
#   
#   return(weight_matrix)
# }


weight_matrix<- function(covariates, knn){
  covariates=as.matrix(covariates)
  
  
  n=nrow(covariates)
  weight_matrix=matrix(0,ncol = n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      weight_matrix[i,j] = sqrt(sum((covariates[i,] - covariates[j,])^2))
      
      if(i==j){
        weight_matrix[i,j]=Inf
      }
    }
    
  }
  
  weight_matrix <- apply(weight_matrix, 1, function(x) ifelse(x <= sort(x)[knn], 1, 0))
  #standardizing
  
  #weight_matrix=(1/max(weight_matrix[weight_matrix!=Inf]))*weight_matrix
  # weight_matrix[weight_matrix==Inf]=1
  weight_matrix=weight_matrix*array(rep(1/colSums(weight_matrix),n),dim = c(n,n))
  
  return(weight_matrix)
}

choose_functions <- function(model){
  
  if (model == "spGARCH")
  {
    
    f_inv <- function(x){
      return(x)
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- eps^2
      n     <- length(eps)
      return(diag(eps_2) %*% solve(diag(n) - Matrix(W_1) %*% diag(eps_2) - Matrix(W_2)) %*% alpha)
    }
    tau_y <- function(y){
      return(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, rhoW_1, lambdaW_2){
      n <- length(eps)
      aux <- solve(diag(n) - rhoW_1 %*% diag(as.vector(eps^2)) - lambdaW_2)
      result <- array(, dim = c(n,n))
      for(j in 1:n){
        aux_rhoW_1 <- rhoW_1
        aux_rhoW_1[,-j] <- array(0, dim = c(n, n-1))
        eps_j <- eps[j]
        result[, j] <- 2 * eps_j * aux %*% aux_rhoW_1 %*% aux %*% alpha
      }
      return(result)
    }

  } 
  else if (model == "spEGARCH")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(log(abs(eps)^b))
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      # is only valid if g(eps) = (ln(|eps_1|^b), ..., ln(|eps_n|^b))
      n     <- length(y)
      # b     <- 2 
      # X     <- solve(diag(n) + Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% log(abs(y)^b))
      X     <- solve(diag(n) + 0.5 * b * Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + b * Matrix(rhoW_1) %*% log(abs(y)))
      eps   <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps){
        return(b / eps) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j) * h_i)
    }
    
  }
  else if (model == "spEGARCH2")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      return(theta * eps)
    }
    tau_y <- function(y){
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n     <- length(y)
      # X     <- solve(diag(n) + Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + Matrix(rhoW_1) %*% log(abs(y)^b))
      # mapping from y to eps (or rather X = log(eps))
      function_y_eps <- function(x, y, alpha, rhoW_1, lambdaW_2, theta, zeta){
        n <- length(x)
        eps <- x
        return(    (diag(n) - lambdaW_2) %*% log(abs(y)^2) - alpha     -     ( rhoW_1 %*% (theta * eps) + (diag(n) - lambdaW_2) %*% log(abs(eps)^2) )  )
      }
      # function_y_eps(x = eps, y = Y, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2)
      out <- nleqslv(x = y, fn = function_y_eps, alpha = alpha, rhoW_1 = rhoW_1, lambdaW_2 = lambdaW_2, theta = theta, zeta = zeta, y = y, control = list(ftol = 1e-10))
      # plot(eps, out$x)
      # plot((out$x), ylim = c(-4,4))
      # points((eps), col = "red")
      # X     <- solve(diag(n) + 0.5 * b * Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% (alpha + b * Matrix(rhoW_1) %*% log(abs(y)))
      # eps   <- y / sqrt(as.vector(f_inv(X)))
      eps <- out$x
      X <- log(abs(y)^2) - log(abs(eps)^2)
      return(list(eps = eps, h =  f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      tau_eps_prime <- function(eps, theta){
        return(theta) # Achtung evtl. zu ändern, hard coded!
      }
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j, theta) * h_i)
    }
    
  }
  else if (model == "spHARCH")
  {
    
    f_inv <- function(x){
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b){  # function tau in paper
      eps_2 <- log(eps^2)
      n     <- length(eps)
      return( solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% alpha + (diag(n) + solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% Matrix(W_1)) %*% eps_2  )
    } 
    tau_y <- function(y){
      log(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv){ # mapping from y to eps
      n   <- length(y)
      X   <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + rhoW_1 %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2){
      n     <- length(eps)
      aux   <- solve(diag(n) - Matrix(rhoW_1) - Matrix(lambdaW_2)) %*% Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i   <- array(rep(h, n), dim = c(n, n))
      return(2 * aux / eps_j * h_i)
    }
  }
  
  return(list(f_inv = f_inv, tau_eps = tau_eps, tau_y = tau_y, g = g, d_h_d_eps = d_h_d_eps))
}



model="spHARCH"

functions <- choose_functions(model)
f_inv     <- functions$f_inv
g         <- functions$g
tau_eps   <- functions$tau_eps
tau_y     <- functions$tau_y
d_h_d_eps <- functions$d_h_d_eps

if(is.element(model, c("spGARCH", "spHARCH"))){
  pars     <- c(0.5, 0.5, 1)
  LB       <- c(0, 0, 0.00001)
  UB       <- c(Inf, Inf, Inf)
} else if(is.element(model, c("spEGARCH"))) {
  pars     <- c(0.5, 0.5, 1)
  LB       <- c(0, 0, 0.00001)
  UB       <- c(Inf, Inf, Inf)
} else if(is.element(model, c("spEGARCH2"))) {
  pars     <- c(0.5, 0.5, 1, 0.5)
  LB       <- c(0, 0, 0.00001, 0.00001)
  UB       <- c(Inf, Inf, Inf, Inf)
}

Least_Squares <- function(pars_K, param){
  
  rho    <- pars_K[1]
  lambda <- pars_K[2]
  alpha  <- pars_K[3]
  theta  <- 1
  zeta   <- 1
  b      <- 1
  
  y         <- param$y
  f_inv     <- param$f_inv
  tau_y     <- param$tau_y
  W_1       <- param$W_1
  W_2       <- param$W_2
  g         <- param$g
  
  n            <- length(y)
  rhoW_1       <- rho * W_1
  lambdaW_2    <- lambda * W_2
  
  result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  h            <- as.vector(result_g[[2]])
  
  z <- log(y^2) - (digamma(1) - log(2))
  h_tilde <- log(h)
  
  squares <- (z - h_tilde)^2    
  
  return(sum(squares))
  
}

LL <- function(pars, param){
  
  model  <- param$model
  
  if(is.element(model, c("spGARCH", "spHARCH"))){
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- 1
    zeta   <- 1
    b      <- 1
  } else if(is.element(model, c("spEGARCH"))) {
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- NULL
    zeta   <- NULL
    b      <- pars[4]
  } else if(is.element(model, c("spEGARCH2"))) {
    rho    <- pars[1]
    lambda <- pars[2]
    alpha  <- pars[3]
    theta  <- pars[4]
    zeta   <- 0 # pars[5]
    b      <- NULL
  }
  
  
  y         <- param$y
  f_inv     <- param$f_inv
  tau_y     <- param$tau_y
  W_1       <- param$W_1
  W_2       <- param$W_2
  g         <- param$g
  d_h_d_eps <- param$d_h_d_eps
  
  n            <- length(y)
  rhoW_1       <- rho * W_1
  lambdaW_2    <- lambda * W_2
  alpha        <- alpha * rep(1, n)
  
  result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  eps          <- as.vector(result_g[[1]])
  h            <- as.vector(result_g[[2]])
  
  J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
  
  # print(pars)
  
  log_det_J <- determinant(J, logarithm = TRUE)$modulus  
  
  return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
  
}

#set number of replications
m = 2000



#Wrapper function for simulation
sim_wrapper <- function(i,n, omega, rho, lambda, alpha){
  #compute true locations
  
#adapt space of with amount of firms/objects in space (see paper for more details)   
  if(n==50){
    S = matrix(runif(n*2,0,1),nrow = n, ncol = 2)
  }else if(n==100){
    S = matrix(runif(n*2,0,sqrt(2)),nrow = n, ncol = 2)
  }else{
    S = matrix(runif(n*2,0,2),nrow = n, ncol = 2)
  }
  knn= 10
  WM_true_1 = weight_matrix(S, knn)
  
  WM_true_2 = weight_matrix(S, 5)
  
  #simulate series
  
  eps   <- rnorm(n)
  
  tau   <- tau_eps(eps, rho * WM_true_1, lambda * WM_true_2, rep(alpha, n), theta, zeta, b)
  Inv   <- solve(diag(n) - lambda * WM_true_2)
  X     <- Inv %*% (alpha + rho * WM_true_1 %*% tau)
  
  sim_series     <- diag(as.vector(f_inv(X)))^(0.5) %*% eps
  
  
  #transform locations with 
  X_star = (1-omega)*S + omega * matrix(rnorm(n*2),nrow = n, ncol = 2)
  
  
  #getting R^2
  # df_lm = data.frame(S,X_star)
  # r1 = lm(X1~X1.1, df_lm)
  # r2 = lm(X2~X2.1, df_lm)
  # r1=summary(r1)$r.squared 
  # r2=summary(r2)$r.squared 
  WM_estimate_1 = weight_matrix(X_star, knn)
  
  WM_estimate_2 = weight_matrix(X_star, 5)
  #estimating parameter
  param    <- list(y = sim_series, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM_estimate_1, W_2 = WM_estimate_2, model = model)
  #Einmal ohne Startwerte und mit Startwerten versuchen
  #Sonst solnp() mit startwerten rho=0.4, lambda=0.4, alpha=0.5
  #out      <- Rsolnp::gosolnp(fun = LL, n.restarts=4, n.sim =30, LB = c(0, 0, 0.00001), UB = c(1, 1, 5), param = param, control = list(trace = FALSE))
  out      <- Rsolnp::solnp(c(rho+0.003,lambda+0.003,alpha+0.005),fun = LL, LB = c(0, 0, 0.00001), UB = c(2, 2, 5), param = param, control = list(trace = FALSE))
  
  # out$pars / sqrt(diag(solve(out$hessian)))
  # round(out$pars, 3)
  
  res = c(out$pars)
  return(res)
}
#set number of cores
nc=6

cl <- makeCluster(nc)
setDefaultCluster(cl)
#setting seed for reproducibility


sim_save <- function(i,omega,rho,lambda, alpha){
  set.seed(i)
  sim_results_setup1 = future.apply::future_lapply(FUN=sim_wrapper, X= 1:m, n=50, omega=omega, rho=rho, lambda=lambda, alpha=alpha, future.seed = (i+1))
  sim_results_setup2 = future.apply::future_lapply(FUN=sim_wrapper, X= 1:m, n=100, omega=omega, rho=rho, lambda=lambda, alpha=alpha, future.seed = (i+1)^2)
  sim_results_setup3 = future.apply::future_lapply(FUN=sim_wrapper, X= 1:m, n=200, omega=omega, rho=rho, lambda=lambda, alpha=alpha, future.seed = (i+1)^3)

  sim_results_setup1=list.rbind(sim_results_setup1)
  sim_results_setup2=list.rbind(sim_results_setup2)
  sim_results_setup3=list.rbind(sim_results_setup3)
  
  assign(paste0("sim_results_setup_",i,"_","50","_",omega), sim_results_setup1)
  assign(paste0("sim_results_setup_",i,"_","100","_",omega), sim_results_setup2)
  assign(paste0("sim_results_setup_",i,"_","200","_",omega), sim_results_setup3)
  
  
  save(list=paste0("sim_results_setup_",i,"_","50","_",omega),file=paste0("sim_results_setup_",i,"_","50","_",omega,".RData"))
  save(list=paste0("sim_results_setup_",i,"_","100","_",omega),file=paste0("sim_results_setup_",i,"_","100","_",omega,".RData"))
  save(list=paste0("sim_results_setup_",i,"_","200","_",omega),file=paste0("sim_results_setup_",i,"_","200","_",omega,".RData"))
  
}

i=1
sim_save(i, omega = 0, rho = 0.1 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.1, rho = 0.1 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.2, rho = 0.1 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.5, rho = 0.1 , lambda = 0.2, alpha = 0.5)

i=2
sim_save(i, omega = 0, rho = 0.2 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.1, rho = 0.2 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.2, rho = 0.2 , lambda = 0.2, alpha = 0.5)
sim_save(i, omega = 0.5, rho = 0.2 , lambda = 0.2, alpha = 0.5)

i=3
rho=0.3
lambda =0.2
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)


i=4
rho=0.4
lambda =0.2
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)



i=5
rho=0.5
lambda =0.2
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)

i=6
rho=0.6
lambda =0.2
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)



i=7
rho=0.7
lambda =0.2
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)



i=8
rho=0.1
lambda =0.4
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)



i=9
rho=0.2
lambda =0.4
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)



i=10
rho=0.3
lambda =0.4
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)


i=11
rho=0.4
lambda =0.4
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)


i=12
rho=0.5
lambda =0.4
sim_save(i, omega = 0, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.1, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.2, rho = rho , lambda = lambda, alpha = 0.5)
sim_save(i, omega = 0.5, rho = rho , lambda = lambda, alpha = 0.5)





ende

set.seed(4572)
sim_results_setup1 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup2 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup3 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1=list.rbind(sim_results_setup1)
sim_results_setup2=list.rbind(sim_results_setup2)
sim_results_setup3=list.rbind(sim_results_setup3)

set.seed(4573)
sim_results_setup1_50_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.1, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_100_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.1, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_150_omega_01 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.1, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1_50_omega_01=list.rbind(sim_results_setup1_50_omega_01)
sim_results_setup1_100_omega_01=list.rbind(sim_results_setup1_100_omega_01)
sim_results_setup1_150_omega_01=list.rbind(sim_results_setup1_150_omega_01)

set.seed(4574)
sim_results_setup1_50_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.2, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_100_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.2, rho=0.7, lambda=0.2, alpha=1)
sim_results_setup1_150_omega_02 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.2, rho=0.7, lambda=0.2, alpha=1)

sim_results_setup1_50_omega_02=list.rbind(sim_results_setup1_50_omega_02)
sim_results_setup1_100_omega_02=list.rbind(sim_results_setup1_100_omega_02)
sim_results_setup1_150_omega_02=list.rbind(sim_results_setup1_150_omega_02)



save(sim_results_setup1_50_omega_02,file="sim_results_setup1_50_omega_02.RData")
save(sim_results_setup1_100_omega_02,file="sim_results_setup1_100_omega_02.RData")
save(sim_results_setup1_150_omega_02,file="sim_results_setup1_150_omega_02.RData")

save(sim_results_setup1_50_omega_01,file="sim_results_setup1_50_omega_01.RData")
save(sim_results_setup1_100_omega_01,file="sim_results_setup1_100_omega_01.RData")
save(sim_results_setup1_150_omega_01,file="sim_results_setup1_150_omega_01.RData")



set.seed(4572)
sim_results_setup2_50 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_100= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_150 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup2_50=list.rbind(sim_results_setup2_50)
sim_results_setup2_100=list.rbind(sim_results_setup2_100)
sim_results_setup2_150=list.rbind(sim_results_setup2_150)

sim_results_setup2_50_omega01 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_100_omega01= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_150_omega01 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup2_50_omega01=list.rbind(sim_results_setup2_50_omega01)
sim_results_setup2_100_omega01=list.rbind(sim_results_setup2_100_omega01)
sim_results_setup2_150_omega01=list.rbind(sim_results_setup2_150_omega01)

sim_results_setup2_50_omega02 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_100_omega02= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup2_150_omega02 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup2_50_omega02=list.rbind(sim_results_setup2_50_omega02)
sim_results_setup2_100_omega02=list.rbind(sim_results_setup2_100_omega02)
sim_results_setup2_150_omega02=list.rbind(sim_results_setup2_150_omega02)

save(sim_results_setup2_50,file="sim_results_setup2_50.RData")
save(sim_results_setup2_100,file="sim_results_setup2_100.RData")
save(sim_results_setup2_150,file="sim_results_setup2_150.RData")


save(sim_results_setup2_50_omega02,file="sim_results_setup2_50_omega_02.RData")
save(sim_results_setup2_100_omega02,file="sim_results_setup2_100_omega_02.RData")
save(sim_results_setup2_150_omega02,file="sim_results_setup2_150_omega_02.RData")

save(sim_results_setup2_50_omega01,file="sim_results_setup2_50_omega_01.RData")
save(sim_results_setup2_100_omega01,file="sim_results_setup2_100_omega_01.RData")
save(sim_results_setup2_150_omega01,file="sim_results_setup2_150_omega_01.RData")


set.seed(4572)
sim_results_setup3_50 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_100= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_150 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup3_50=list.rbind(sim_results_setup2_50)
sim_results_setup3_100=list.rbind(sim_results_setup2_100)
sim_results_setup3_150=list.rbind(sim_results_setup2_150)

sim_results_setup3_50_omega01 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_100_omega01= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_150_omega01 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.1, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup3_50_omega01=list.rbind(sim_results_setup2_50_omega01)
sim_results_setup3_100_omega01=list.rbind(sim_results_setup2_100_omega01)
sim_results_setup3_150_omega01=list.rbind(sim_results_setup2_150_omega01)

sim_results_setup3_50_omega02 = pblapply(FUN=sim_wrapper, X= 1:m, n=50, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_100_omega02= pblapply(FUN=sim_wrapper, X= 1:m, n=100, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)
sim_results_setup3_150_omega02 = pblapply(FUN=sim_wrapper, X= 1:m, n=200, omega=0.2, rho=0.3, lambda=0.4, alpha=0.5)

sim_results_setup3_50_omega02=list.rbind(sim_results_setup3_50_omega02)
sim_results_setup3_100_omega02=list.rbind(sim_results_setup3_100_omega02)
sim_results_setup3_150_omega02=list.rbind(sim_results_setup3_150_omega02)

save(sim_results_setup3_50,file="sim_results_setup3_50.RData")
save(sim_results_setup3_100,file="sim_results_setup3_100.RData")
save(sim_results_setup3_150,file="sim_results_setup3_150.RData")


save(sim_results_setup3_50_omega02,file="sim_results_setup3_50_omega_02.RData")
save(sim_results_setup3_100_omega02,file="sim_results_setup3_100_omega_02.RData")
save(sim_results_setup3_150_omega02,file="sim_results_setup3_150_omega_02.RData")

save(sim_results_setup3_50_omega01,file="sim_results_setup3_50_omega_01.RData")
save(sim_results_setup3_100_omega01,file="sim_results_setup3_100_omega_01.RData")
save(sim_results_setup3_150_omega01,file="sim_results_setup3_150_omega_01.RData")


stopCluster(cl)


#evaluation for table in ascending order by n = 50, 100, 200 and mean, sd, mse (respectively)
sim_evaluation<- function(sim_results_50, sim_results_100, sim_results_200, true_val){
mean_50 = sapply(1:3,function(x){round(mean(sim_results_50[,x]),3)})
mean_100 = sapply(1:3,function(x){round(mean(sim_results_100[,x]),3)})
mean_200 = sapply(1:3,function(x){round(mean(sim_results_200[,x]),3)})

sd_50 = sapply(1:3,function(x){round(sd(sim_results_50[,x]),3)})
sd_100 = sapply(1:3,function(x){round(sd(sim_results_100[,x]),3)})
sd_200 = sapply(1:3,function(x){round(sd(sim_results_200[,x]),3)})

mse_50 = c(round(sqrt(mean((sim_results_50[,1]-true_val[1])^2)),3),round(sqrt(mean((sim_results_50[,2]-true_val[2])^2)),3),round(sqrt(mean(sim_results_50[,3]-true_val[3])^2),3))
mse_100 = c(round(sqrt(mean((sim_results_100[,1]-true_val[1])^2)),3),round(sqrt(mean((sim_results_100[,2]-true_val[2])^2)),3),round(sqrt(mean(sim_results_100[,3]-true_val[3])^2),3))
mse_200 = c(round(sqrt(mean((sim_results_200[,1]-true_val[1])^2)),3),round(sqrt(mean((sim_results_200[,2]-true_val[2])^2)),3),round(sqrt(mean(sim_results_200[,3]-true_val[3])^2),3))

print(paste(mean_50[1],sd_50[1], mse_50[1], mean_100[1],sd_100[1],mse_100[1],mean_200[1], sd_200[1], mse_200[1], sep = " & "))
paste("\n")
print(paste(mean_50[2],sd_50[2], mse_50[2], mean_100[2],sd_100[2],mse_100[2],mean_200[2], sd_200[2], mse_200[2], sep = " & "))
paste("\n")
paste(mean_50[3],sd_50[3], mse_50[3], mean_100[3],sd_100[3],mse_100[3],mean_200[3], sd_200[3], mse_200[3], sep = " & ")
}

#setup 1
true_val = c(0.1, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_1_50_0.RData")
load("~/spatialGARCH/sim_results_setup_1_100_0.RData")
load("~/spatialGARCH/sim_results_setup_1_200_0.RData")
sim_evaluation(sim_results_setup_1_50_0,sim_results_setup_1_100_0,sim_results_setup_1_200_0, true_val)

data_tmp_1 = data.frame(sim_results_setup_1_50_0, rep(50,1000) )
data_tmp_2 = data.frame(sim_results_setup_1_100_0, rep(100,1000) )
data_tmp_3 = data.frame(sim_results_setup_1_200_0, rep(200,1000) )
colnames(data_tmp_1) =  c("rho", "lambda", "alpha", "n")
colnames(data_tmp_2) =  c("rho", "lambda", "alpha", "n")
colnames(data_tmp_3) =  c("rho", "lambda", "alpha", "n")
data = rbind(data_tmp_1,data_tmp_2,data_tmp_3)
data$n = as.factor(data$n)
library(ggplot2)
library(gridExtra)
# Basic box plot
p1 <- ggplot(data, aes(x=n, y=rho)) + 
  geom_boxplot() + geom_hline(yintercept=true_val[1], linetype="dashed", color = "red")
p2 <- ggplot(data, aes(x=n, y=lambda)) + 
  geom_boxplot() + geom_hline(yintercept=true_val[2], linetype="dashed", color = "red")
p3 <- ggplot(data, aes(x=n, y=alpha)) + 
  geom_boxplot() + geom_hline(yintercept=true_val[3], linetype="dashed", color = "red")


grid.arrange(p1,p2,p3)


load("~/spatialGARCH/sim_results_setup_1_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_1_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_1_200_0.1.RData")
sim_evaluation(sim_results_setup_1_50_0.1,sim_results_setup_1_100_0.1,sim_results_setup_1_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_1_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_1_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_1_200_0.2.RData")
sim_evaluation(sim_results_setup_1_50_0.2, sim_results_setup_1_100_0.2, sim_results_setup_1_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_1_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_1_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_1_200_0.5.RData")
sim_evaluation(sim_results_setup_1_50_0.5, sim_results_setup_1_100_0.5, sim_results_setup_1_200_0.5, true_val)

#setup 2
true_val = c(0.2, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_2_50_0.RData")
load("~/spatialGARCH/sim_results_setup_2_100_0.RData")
load("~/spatialGARCH/sim_results_setup_2_200_0.RData")
sim_evaluation(sim_results_setup_2_50_0,sim_results_setup_2_100_0,sim_results_setup_2_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_2_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_2_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_2_200_0.1.RData")
sim_evaluation(sim_results_setup_2_50_0.1,sim_results_setup_2_100_0.1,sim_results_setup_2_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_2_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_2_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_2_200_0.2.RData")
sim_evaluation(sim_results_setup_2_50_0.2, sim_results_setup_2_100_0.2, sim_results_setup_2_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_2_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_2_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_2_200_0.5.RData")
sim_evaluation(sim_results_setup_2_50_0.5, sim_results_setup_2_100_0.5, sim_results_setup_2_200_0.5, true_val)

#setup 3
true_val = c(0.3, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_3_50_0.RData")
load("~/spatialGARCH/sim_results_setup_3_100_0.RData")
load("~/spatialGARCH/sim_results_setup_3_200_0.RData")
sim_evaluation(sim_results_setup_3_50_0,sim_results_setup_3_100_0,sim_results_setup_3_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_3_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_3_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_3_200_0.1.RData")
sim_evaluation(sim_results_setup_3_50_0.1,sim_results_setup_3_100_0.1,sim_results_setup_3_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_3_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_3_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_3_200_0.2.RData")
sim_evaluation(sim_results_setup_3_50_0.2, sim_results_setup_3_100_0.2, sim_results_setup_3_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_3_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_3_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_3_200_0.5.RData")
sim_evaluation(sim_results_setup_3_50_0.5, sim_results_setup_3_100_0.5, sim_results_setup_3_200_0.5, true_val)


#setup 4
true_val = c(0.4, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_4_50_0.RData")
load("~/spatialGARCH/sim_results_setup_4_100_0.RData")
load("~/spatialGARCH/sim_results_setup_4_200_0.RData")
sim_evaluation(sim_results_setup_4_50_0,sim_results_setup_4_100_0,sim_results_setup_4_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_4_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_4_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_4_200_0.1.RData")
sim_evaluation(sim_results_setup_4_50_0.1,sim_results_setup_4_100_0.1,sim_results_setup_4_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_4_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_4_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_4_200_0.2.RData")
sim_evaluation(sim_results_setup_4_50_0.2, sim_results_setup_4_100_0.2, sim_results_setup_4_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_4_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_4_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_4_200_0.5.RData")
sim_evaluation(sim_results_setup_4_50_0.5, sim_results_setup_4_100_0.5, sim_results_setup_4_200_0.5, true_val)


#setup 5
true_val = c(0.5, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_5_50_0.RData")
load("~/spatialGARCH/sim_results_setup_5_100_0.RData")
load("~/spatialGARCH/sim_results_setup_5_200_0.RData")
sim_evaluation(sim_results_setup_5_50_0,sim_results_setup_5_100_0,sim_results_setup_5_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_5_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_5_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_5_200_0.1.RData")
sim_evaluation(sim_results_setup_5_50_0.1,sim_results_setup_5_100_0.1,sim_results_setup_5_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_5_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_5_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_5_200_0.2.RData")
sim_evaluation(sim_results_setup_5_50_0.2, sim_results_setup_5_100_0.2, sim_results_setup_5_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_5_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_5_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_5_200_0.5.RData")
sim_evaluation(sim_results_setup_5_50_0.5, sim_results_setup_5_100_0.5, sim_results_setup_5_200_0.5, true_val)



#setup 6
true_val = c(0.6, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_6_50_0.RData")
load("~/spatialGARCH/sim_results_setup_6_100_0.RData")
load("~/spatialGARCH/sim_results_setup_6_200_0.RData")
sim_evaluation(sim_results_setup_6_50_0,sim_results_setup_6_100_0,sim_results_setup_6_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_6_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_6_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_6_200_0.1.RData")
sim_evaluation(sim_results_setup_6_50_0.1,sim_results_setup_6_100_0.1,sim_results_setup_6_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_6_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_6_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_6_200_0.2.RData")
sim_evaluation(sim_results_setup_6_50_0.2, sim_results_setup_6_100_0.2, sim_results_setup_6_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_6_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_6_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_6_200_0.5.RData")
sim_evaluation(sim_results_setup_6_50_0.5, sim_results_setup_6_100_0.5, sim_results_setup_6_200_0.5, true_val)

#setup 7
true_val = c(0.7, 0.2, 0.5)
load("~/spatialGARCH/sim_results_setup_7_50_0.RData")
load("~/spatialGARCH/sim_results_setup_7_100_0.RData")
load("~/spatialGARCH/sim_results_setup_7_200_0.RData")
sim_evaluation(sim_results_setup_7_50_0,sim_results_setup_7_100_0,sim_results_setup_7_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_7_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_7_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_7_200_0.1.RData")
sim_evaluation(sim_results_setup_7_50_0.1,sim_results_setup_7_100_0.1,sim_results_setup_7_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_7_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_7_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_7_200_0.2.RData")
sim_evaluation(sim_results_setup_7_50_0.2, sim_results_setup_7_100_0.2, sim_results_setup_7_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_7_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_7_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_7_200_0.5.RData")
sim_evaluation(sim_results_setup_7_50_0.5, sim_results_setup_7_100_0.5, sim_results_setup_7_200_0.5, true_val)

#setup 8
true_val = c(0.1, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_8_50_0.RData")
load("~/spatialGARCH/sim_results_setup_8_100_0.RData")
load("~/spatialGARCH/sim_results_setup_8_200_0.RData")
sim_evaluation(sim_results_setup_8_50_0,sim_results_setup_8_100_0,sim_results_setup_8_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_8_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_8_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_8_200_0.1.RData")
sim_evaluation(sim_results_setup_8_50_0.1,sim_results_setup_8_100_0.1,sim_results_setup_8_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_8_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_8_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_8_200_0.2.RData")
sim_evaluation(sim_results_setup_8_50_0.2, sim_results_setup_8_100_0.2, sim_results_setup_8_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_8_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_8_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_8_200_0.5.RData")
sim_evaluation(sim_results_setup_8_50_0.5, sim_results_setup_8_100_0.5, sim_results_setup_8_200_0.5, true_val)

#setup 9
true_val = c(0.2, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_9_50_0.RData")
load("~/spatialGARCH/sim_results_setup_9_100_0.RData")
load("~/spatialGARCH/sim_results_setup_9_200_0.RData")
sim_evaluation(sim_results_setup_9_50_0,sim_results_setup_9_100_0,sim_results_setup_9_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_9_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_9_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_9_200_0.1.RData")
sim_evaluation(sim_results_setup_9_50_0.1,sim_results_setup_9_100_0.1,sim_results_setup_9_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_9_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_9_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_9_200_0.2.RData")
sim_evaluation(sim_results_setup_9_50_0.2, sim_results_setup_9_100_0.2, sim_results_setup_9_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_9_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_9_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_9_200_0.5.RData")
sim_evaluation(sim_results_setup_9_50_0.5, sim_results_setup_9_100_0.5, sim_results_setup_9_200_0.5, true_val)

#setup 10
true_val = c(0.3, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_10_50_0.RData")
load("~/spatialGARCH/sim_results_setup_10_100_0.RData")
load("~/spatialGARCH/sim_results_setup_10_200_0.RData")
sim_evaluation(sim_results_setup_10_50_0,sim_results_setup_10_100_0,sim_results_setup_10_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_10_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_10_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_10_200_0.1.RData")
sim_evaluation(sim_results_setup_10_50_0.1,sim_results_setup_10_100_0.1,sim_results_setup_10_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_10_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_10_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_10_200_0.2.RData")
sim_evaluation(sim_results_setup_10_50_0.2, sim_results_setup_10_100_0.2, sim_results_setup_10_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_10_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_10_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_10_200_0.5.RData")
sim_evaluation(sim_results_setup_10_50_0.5, sim_results_setup_10_100_0.5, sim_results_setup_10_200_0.5, true_val)

#setup 11
true_val = c(0.4, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_11_50_0.RData")
load("~/spatialGARCH/sim_results_setup_11_100_0.RData")
load("~/spatialGARCH/sim_results_setup_11_200_0.RData")
sim_evaluation(sim_results_setup_11_50_0,sim_results_setup_11_100_0,sim_results_setup_11_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_11_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_11_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_11_200_0.1.RData")
sim_evaluation(sim_results_setup_11_50_0.1,sim_results_setup_11_100_0.1,sim_results_setup_11_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_11_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_11_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_11_200_0.2.RData")
sim_evaluation(sim_results_setup_11_50_0.2, sim_results_setup_11_100_0.2, sim_results_setup_11_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_11_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_11_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_11_200_0.5.RData")
sim_evaluation(sim_results_setup_11_50_0.5, sim_results_setup_11_100_0.5, sim_results_setup_11_200_0.5, true_val)


#setup 12
true_val = c(0.5, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_12_50_0.RData")
load("~/spatialGARCH/sim_results_setup_12_100_0.RData")
load("~/spatialGARCH/sim_results_setup_12_200_0.RData")
sim_evaluation(sim_results_setup_12_50_0,sim_results_setup_12_100_0,sim_results_setup_12_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_12_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_12_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_12_200_0.1.RData")
sim_evaluation(sim_results_setup_12_50_0.1,sim_results_setup_12_100_0.1,sim_results_setup_12_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_12_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_12_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_12_200_0.2.RData")
sim_evaluation(sim_results_setup_12_50_0.2, sim_results_setup_12_100_0.2, sim_results_setup_12_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_12_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_12_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_12_200_0.5.RData")
sim_evaluation(sim_results_setup_12_50_0.5, sim_results_setup_12_100_0.5, sim_results_setup_12_200_0.5, true_val)

#setup 13
true_val = c(0.6, 0.4, 0.5)
load("~/spatialGARCH/sim_results_setup_13_50_0.RData")
load("~/spatialGARCH/sim_results_setup_13_100_0.RData")
load("~/spatialGARCH/sim_results_setup_13_200_0.RData")
sim_evaluation(sim_results_setup_13_50_0,sim_results_setup_13_100_0,sim_results_setup_13_200_0, true_val)

load("~/spatialGARCH/sim_results_setup_13_50_0.1.RData")
load("~/spatialGARCH/sim_results_setup_13_100_0.1.RData")
load("~/spatialGARCH/sim_results_setup_13_200_0.1.RData")
sim_evaluation(sim_results_setup_13_50_0.1,sim_results_setup_13_100_0.1,sim_results_setup_13_200_0.1, true_val)

load("~/spatialGARCH/sim_results_setup_13_50_0.2.RData")
load("~/spatialGARCH/sim_results_setup_13_100_0.2.RData")
load("~/spatialGARCH/sim_results_setup_13_200_0.2.RData")
sim_evaluation(sim_results_setup_13_50_0.2, sim_results_setup_13_100_0.2, sim_results_setup_13_200_0.2, true_val)

load("~/spatialGARCH/sim_results_setup_13_50_0.5.RData")
load("~/spatialGARCH/sim_results_setup_13_100_0.5.RData")
load("~/spatialGARCH/sim_results_setup_13_200_0.5.RData")
sim_evaluation(sim_results_setup_13_50_0.5, sim_results_setup_13_100_0.5, sim_results_setup_13_200_0.5, true_val)



