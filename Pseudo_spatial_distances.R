#spatial relation by differences 
#Input vector of covariates (here col vector)
library("Rsolnp")
library("Matrix")
library("Rcpp")
library("xts")
#set directory may to be changed
path.expand("~")
load("~/spatialGARCH/workspace.RData")

#inverse distance matrix as alternative - but not used in the current version of the paper!
weight_matrix<-function(covariates,knn){
  covariates=as.matrix(covariates)
  n=nrow(covariates)
  weight_matrix=matrix(0,ncol = n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      weight_matrix[i,j] = 1/sqrt(sum((covariates[i,] - covariates[j,])^2))

      if(i==j){
        weight_matrix[i,j]=0
      }
    }

  }
  weight_matrix[weight_matrix==Inf]=max(weight_matrix[weight_matrix!=Inf])
  #standardizing

  #weight_matrix[weight_matrix==Inf]=max(weight_matrix[weight_matrix!=Inf])
  
  
  weight_matrix=1/max(eigen(weight_matrix)$values)*weight_matrix
  #weight_matrix=(1/max(weight_matrix[weight_matrix!=Inf]))*weight_matrix
  #weight_matrix[weight_matrix==Inf]=1
  #print(max(weight_matrix))
  return(weight_matrix)
}

#KNN weight matrix - currently used
weight_matrix <- function(covariates,knn){
  covariates=as.matrix(covariates)
  knn = knn
  
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


#data cleanage of covariates (normalization)
A=data.matrix(financials)
#receive odd indices
index=1:nrow(A)
index=index[index%%2 !=0]
A_2019=A[index,]

A_2019=A_2019[,2:68]
A_2019=A_2019[,c(1:60,62:67)]

A_2019["Samsung_Electronics",]=A_2019["Samsung_Electronics",]/1000
A_2019[18,] = A_2019[18,]/1000
A_2019[30,] = A_2019[30,]/1000
A_2019[47,] = A_2019[47,]/1000
A_2019[27,] = A_2019[27,]/1000
A_2019[23,] = A_2019[23,]/1000
A_2019[12,] = A_2019[12,]/1000
A_2019[8,] = A_2019[8,]/1000
A_2019[7,] = A_2019[7,]/1000
max(A_2019,na.rm = T)

for(j in 1:66){
print(c(round(mean(A_2019[,j] /10^9,na.rm = T),3),round(median(A_2019[,j] /10^9,na.rm = T),3),sd(A_2019[,j] /10^9,na.rm = T), round(range(A_2019[,j] /10^9,na.rm = T),3)))
}
plot(A_2019[,3])
A_stand <- apply(A_2019, 2,function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})

A_stand[is.na(A_stand)]=0
A_stand=as.matrix(A_stand)




#change to complete data
#Log returns generation
stock_prices=function(x){as.matrix(c(Allianz[272-x,1],Abbott_Laboratories[272-x,1],Apple[272-x,1], ATT[272-x,1], Banco_Santander[272-x,1],Bank_of_America[272-x,1],BHP_Billiton[272-x,1],BNP_Paribas[(length(BNP_Paribas[,1]))-x,1],BP[272-x,1],Chevron[272-x,1],Cisco[272-x,1],Citigroup[272-x,1],Coca_Cola[272-x,1],ConocoPhillips[272-x,1],E.ON[272-x,1],Eni[272-x,1],Exxon_Mobil[272-x,1],Gazprom[272-x,1],General_Electric[272-x,1],GlaxoSmithKline[272-x,1],Google[272-x,1],Hewlett_Packard[272-x,1],HSBC[272-x,1],IBM[272-x,1],Intel[272-x,1],Johnson_Johnson[272-x,1],JPMorgan_Chase_Co.[272-x,1],Merck_Co_Inc.[272-x,1],Microsoft[272-x,1],Mitsubishi_UFJ_Financial_Group[272-x,1],Nestlé[272-x,1],Novartis[272-x,1],Oracle[272-x,1],PepsiCo[272-x,1],Petrobras[272-x,1],Pfizer[272-x,1],Philip_Morris_International[272-x,1],Procter_Gamble[272-x,1],Roche[272-x,1],Royal_Dutch_Shell[272-x,1],Samsung_Electronics[length(Samsung_Electronics[,1])-x,1],Sanofi_Aventis[272-x,1],Schlumberger[272-x,1],Siemens[272-x,1],Telefónica[272-x,1],Total[272-x,1],Toyota[272-x,1],Verizon_Communications[272-x,1],Vodafone_Group[272-x,1],Wal_Mart[272-x,1]))}
y=stock_prices(1)
y=log(stock_prices(0))-log(stock_prices(15))
y=as.numeric(y-mean(y))
returns_plot=data.frame(y, row.names(A_2019))


#Plotting of returns
require(ggplot2)
require(reshape2)
library("igraph")


ggplot(data = returns_plot, aes(y=row.names.A_2019., x=y)) + 
  geom_point()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Returns") +ylab("") +
  theme(
        panel.border = element_blank())
melt(returns_plot)



#new optimization
LLV=numeric(66)
#set knn to 5. May play around a bit..
knn = 5
#optimization using all 66 covariates for WM 
for(i in 1:66){
  #getting WM
  WM=weight_matrix(A_stand[,c(i)],knn)
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
    pars     <- c(0.3, 0.3, 1)
    LB       <- c(0, 0, 0.00001)
    UB       <- c(1, 1, Inf)
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
  
 # LL <- function(pars, param){
 #    
 #    model  <- param$model
 #    
 #    if(is.element(model, c("spGARCH", "spHARCH"))){
 #      rho    <- pars[1]
 #      lambda <- pars[2]
 #      alpha  <- pars[3]
 #      theta  <- 1
 #      zeta   <- 1
 #      b      <- 1
 #    } else if(is.element(model, c("spEGARCH"))) {
 #      rho    <- pars[1]
 #      lambda <- pars[2]
 #      alpha  <- pars[3]
 #      theta  <- NULL
 #      zeta   <- NULL
 #      b      <- pars[4]
 #    } else if(is.element(model, c("spEGARCH2"))) {
 #      rho    <- pars[1]
 #      lambda <- pars[2]
 #      alpha  <- pars[3]
 #      theta  <- pars[4]
 #      zeta   <- 0 # pars[5]
 #      b      <- NULL
 #    }
 #    
 #    
 #    y         <- param$y
 #    f_inv     <- param$f_inv
 #    tau_y     <- param$tau_y
 #    W_1       <- param$W_1
 #    W_2       <- param$W_2
 #    g         <- param$g
 #    d_h_d_eps <- param$d_h_d_eps
 #    
 #    n            <- length(y)
 #    rhoW_1       <- rho * W_1
 #    lambdaW_2    <- lambda * W_2
 #    alpha        <- alpha * rep(1, n)
 #    
 #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
 #    eps          <- as.vector(result_g[[1]])
 #    h            <- as.vector(result_g[[2]])
 #    
 #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
 #    
 #    # print(pars)
 #    
 #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
 #    
 #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
 #    
 #  }
  
#   model="spGARCH"
#   functions
#   functions <- choose_functions(model)
#   f_inv     <- functions$f_inv
#   g         <- functions$g
#   tau_eps   <- functions$tau_eps
#   tau_y     <- functions$tau_y
#   d_h_d_eps <- functions$d_h_d_eps
#   
#   if(is.element(model, c("spGARCH", "spHARCH"))){
#     pars     <- c(0.5, 0.5, 1)
#     LB       <- c(0, 0, 0.00001)
#     UB       <- c(Inf, Inf, Inf)
#   } else if(is.element(model, c("spEGARCH"))) {
#     pars     <- c(0.5, 0.5, 1)
#     LB       <- c(0, 0, 0.00001)
#     UB       <- c(Inf, Inf, Inf)
#   } else if(is.element(model, c("spEGARCH2"))) {
#     pars     <- c(0.5, 0.5, 1, 0.5)
#     LB       <- c(0, 0, 0.00001, 0.00001)
#     UB       <- c(Inf, Inf, Inf, Inf)
#   }
#   
   #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
  #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
#   if(out$values[out$outer.iter]<max_lik){
#     max_index=i
#     max_lik=out$values[out$outer.iter]
#     max_obj=out
#   }
  param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
    
    if((rho+lambda)>=1){
      return(1e24)
    }
    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    
    # print(pars)
    
    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    
    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    
  }
  if( is.na(LL(pars,param))){
    pars= c(0.45,0.45,2)
  }
  #Optim step
  out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
  LL(pars,param)
  #LLV[i] = 2*out$values[out$outer.iter] + 4* log(50)
  LLV[i] = -out$values[out$outer.iter] 
  
}
round(LLV,3)
max(LLV)
colnames(A_stand)[which.max(LLV)]
colnames(A_stand)[which.max(LLV)]
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

#getting index of covariate with highest LLV
max_index = which.max(LLV)

#max_index = 60:66

############# Second step combining optimal single covariate pairwise with all remaining
LLV_2 = numeric(66)
for(i in 1:66){
  if (i==max_index){
    LLV_2[i]=-2500
  }else{
  WM=weight_matrix(A_stand[,c(max_index,i)],knn)
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
    pars     <- c(0.3, 0.3, 1)
    LB       <- c(0, 0, 0.00001)
    UB       <- c(1, 1, Inf)
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
  
  # LL <- function(pars, param){
  #    
  #    model  <- param$model
  #    
  #    if(is.element(model, c("spGARCH", "spHARCH"))){
  #      rho    <- pars[1]
  #      lambda <- pars[2]
  #      alpha  <- pars[3]
  #      theta  <- 1
  #      zeta   <- 1
  #      b      <- 1
  #    } else if(is.element(model, c("spEGARCH"))) {
  #      rho    <- pars[1]
  #      lambda <- pars[2]
  #      alpha  <- pars[3]
  #      theta  <- NULL
  #      zeta   <- NULL
  #      b      <- pars[4]
  #    } else if(is.element(model, c("spEGARCH2"))) {
  #      rho    <- pars[1]
  #      lambda <- pars[2]
  #      alpha  <- pars[3]
  #      theta  <- pars[4]
  #      zeta   <- 0 # pars[5]
  #      b      <- NULL
  #    }
  #    
  #    
  #    y         <- param$y
  #    f_inv     <- param$f_inv
  #    tau_y     <- param$tau_y
  #    W_1       <- param$W_1
  #    W_2       <- param$W_2
  #    g         <- param$g
  #    d_h_d_eps <- param$d_h_d_eps
  #    
  #    n            <- length(y)
  #    rhoW_1       <- rho * W_1
  #    lambdaW_2    <- lambda * W_2
  #    alpha        <- alpha * rep(1, n)
  #    
  #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
  #    eps          <- as.vector(result_g[[1]])
  #    h            <- as.vector(result_g[[2]])
  #    
  #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
  #    
  #    # print(pars)
  #    
  #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
  #    
  #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
  #    
  #  }
  
  #   model="spGARCH"
  #   functions
  #   functions <- choose_functions(model)
  #   f_inv     <- functions$f_inv
  #   g         <- functions$g
  #   tau_eps   <- functions$tau_eps
  #   tau_y     <- functions$tau_y
  #   d_h_d_eps <- functions$d_h_d_eps
  #   
  #   if(is.element(model, c("spGARCH", "spHARCH"))){
  #     pars     <- c(0.5, 0.5, 1)
  #     LB       <- c(0, 0, 0.00001)
  #     UB       <- c(Inf, Inf, Inf)
  #   } else if(is.element(model, c("spEGARCH"))) {
  #     pars     <- c(0.5, 0.5, 1)
  #     LB       <- c(0, 0, 0.00001)
  #     UB       <- c(Inf, Inf, Inf)
  #   } else if(is.element(model, c("spEGARCH2"))) {
  #     pars     <- c(0.5, 0.5, 1, 0.5)
  #     LB       <- c(0, 0, 0.00001, 0.00001)
  #     UB       <- c(Inf, Inf, Inf, Inf)
  #   }
  #   
  #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
  #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
  #   if(out$values[out$outer.iter]<max_lik){
  #     max_index=i
  #     max_lik=out$values[out$outer.iter]
  #     max_obj=out
  #   }
  param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
    
    if((rho+lambda)>=1){
      return(1e24)
    }
    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    
    # print(pars)
    
    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    
    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    
  }
  if( is.na(LL(pars,param))){
    pars= c(0.45,0.45,1)
  }  
  out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
  
  LLV_2[i] = -out$values[out$outer.iter]
  }
}
max(LLV_2)
round(LLV_2,3)
which.max(LLV_2)
colnames(A_stand)[which.max(LLV_2)]
#getting index of covariate with highest LLV
max_index=c(max_index,which.max(LLV_2))

############# Third step combining optimal previous covariate combination pairwise with all remaining
LLV_3=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_3[i]=-2500
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)],knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
   
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_3[i] = -out$values[out$outer.iter]
  }
}

round(LLV_3,3)
max(LLV_3)
which.max(LLV_3)

max_index=c(max_index,which.max(LLV_3))
colnames(A_stand)[max_index]

############# Third step combining optimal previous covariate combination pairwise with all remaining

LLV_4=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_4[i]=-2000
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)],knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
    
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_4[i] = -out$values[out$outer.iter]
  }
}
max(LLV_4)
round(LLV_4,3)
max_index=c(max_index,which.max(LLV_4))
max_obj$pars 
colnames(A_stand)[max_index]

############# Third step combining optimal previous covariate combination pairwise with all remaining

LLV_5=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_5[i]=-2000
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)],knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
    
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_5[i] = -out$values[out$outer.iter]
  }
}
max(LLV_5)
round(LLV_5,3)
max_index=c(max_index,which.max(LLV_5))

############# Third step combining optimal previous covariate combination pairwise with all remaining

LLV_6=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_6[i]=-2000
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)], knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
    
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_6[i] = -out$values[out$outer.iter]
  }
}
round(LLV_6,3)
max(LLV_6)
max_index=c(max_index,which.max(LLV_6))

############# Third step combining optimal previous covariate combination pairwise with all remaining

LLV_7=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_7[i]=-1000
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)], knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
    
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_7[i] = -out$values[out$outer.iter]
  }
}
round(LLV_7,3)
which.max(LLV_7)
colnames(A_stand)[38]
max_index=c(max_index,which.max(LLV_7))

############# Third step combining optimal previous covariate combination pairwise with all remaining

LLV_8=numeric(66)
for(i in 1:66){
  if (any(i==max_index)){
    LLV_8[i]=-1000
  }else{
    WM=weight_matrix(A_stand[,c(max_index,i)], knn)
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
      pars     <- c(0.2, 0.5, 1)
      LB       <- c(0, 0, 0)
      UB       <- c(1, 1, Inf)
    } else if(is.element(model, c("spEGARCH"))) {
      pars     <- c(0.5, 0.5, 1)
      LB       <- c(0, 0, 0.00001)
      UB       <- c(Inf, Inf, Inf)
    } else if(is.element(model, c("spEGARCH2"))) {
      pars     <- c(0.5, 0.5, 1, 0.5)
      LB       <- c(0, 0, 0.00001, 0.00001)
      UB       <- c(Inf, Inf, Inf, Inf)
    }
    
    
    
    # LL <- function(pars, param){
    #    
    #    model  <- param$model
    #    
    #    if(is.element(model, c("spGARCH", "spHARCH"))){
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- 1
    #      zeta   <- 1
    #      b      <- 1
    #    } else if(is.element(model, c("spEGARCH"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- NULL
    #      zeta   <- NULL
    #      b      <- pars[4]
    #    } else if(is.element(model, c("spEGARCH2"))) {
    #      rho    <- pars[1]
    #      lambda <- pars[2]
    #      alpha  <- pars[3]
    #      theta  <- pars[4]
    #      zeta   <- 0 # pars[5]
    #      b      <- NULL
    #    }
    #    
    #    
    #    y         <- param$y
    #    f_inv     <- param$f_inv
    #    tau_y     <- param$tau_y
    #    W_1       <- param$W_1
    #    W_2       <- param$W_2
    #    g         <- param$g
    #    d_h_d_eps <- param$d_h_d_eps
    #    
    #    n            <- length(y)
    #    rhoW_1       <- rho * W_1
    #    lambdaW_2    <- lambda * W_2
    #    alpha        <- alpha * rep(1, n)
    #    
    #    result_g     <- g(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv)
    #    eps          <- as.vector(result_g[[1]])
    #    h            <- as.vector(result_g[[2]])
    #    
    #    J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
    #    
    #    # print(pars)
    #    
    #    log_det_J <- determinant(J, logarithm = TRUE)$modulus  
    #    
    #    return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
    #    
    #  }
    
    #   model="spGARCH"
    #   functions
    #   functions <- choose_functions(model)
    #   f_inv     <- functions$f_inv
    #   g         <- functions$g
    #   tau_eps   <- functions$tau_eps
    #   tau_y     <- functions$tau_y
    #   d_h_d_eps <- functions$d_h_d_eps
    #   
    #   if(is.element(model, c("spGARCH", "spHARCH"))){
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH"))) {
    #     pars     <- c(0.5, 0.5, 1)
    #     LB       <- c(0, 0, 0.00001)
    #     UB       <- c(Inf, Inf, Inf)
    #   } else if(is.element(model, c("spEGARCH2"))) {
    #     pars     <- c(0.5, 0.5, 1, 0.5)
    #     LB       <- c(0, 0, 0.00001, 0.00001)
    #     UB       <- c(Inf, Inf, Inf, Inf)
    #   }
    #   
    #param    <- list(y = y, f_inv = f_inv, tau_y = tau_y, d_h_d_eps = d_h_d_eps, g = g, W_1 = WM, W_2 = WM, model = model, pars_K_minus_1 = pars)
    #out      <- solnp(pars, fun = Least_Squares, LB = LB, UB = UB, param = param, control = list(trace = FALSE))
    #   if(out$values[out$outer.iter]<max_lik){
    #     max_index=i
    #     max_lik=out$values[out$outer.iter]
    #     max_obj=out
    #   }
    param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
      
      if((rho+lambda)>=1){
        return(1e24)
      }
      J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
      
      # print(pars)
      
      log_det_J <- determinant(J, logarithm = TRUE)$modulus  
      
      return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
      
    }
    if( is.na(LL(pars,param)) || is.infinite(LL(pars,param))){
      pars= c(0.25,0.45,0.5)
    }  
    out      <- solnp(pars, fun = LL, LB = c(0, 0, 0.00001), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))
    
    LLV_8[i] = -out$values[out$outer.iter]
  }
}



max(LLV_7)
colnames(A_stand)[max_index]
round(sqrt(diag(solve(max_obj$hessian))),digits=3)

max_obj$pars/sqrt(diag(solve(max_obj$hessian)))



residuals_spGARCH <- function(pars_K, param){
  
  model  <- param$model
  pars_K_minus_1 <- param$pars_K_minus_1
  
  if(is.element(model, c("spGARCH", "spHARCH"))){
    rho    <- pars_K[1]
    lambda <- pars_K[2]
    alpha  <- pars_K[3]
    theta  <- 1
    zeta   <- 1
    b      <- 1
    # 
    rho_k    <- pars_K_minus_1[1]
    lambda_k <- pars_K_minus_1[2]
    alpha_k  <- pars_K_minus_1[3]
    theta_k  <- 1
    zeta_k   <- 1
    b_k      <- 1
  } else if(is.element(model, c("spEGARCH"))) {
    rho    <- pars_K[1]
    lambda <- pars_K[2]
    alpha  <- pars_K[3]
    theta  <- NULL
    zeta   <- NULL
    b      <- 2 # pars[4]
    # 
    rho_k    <- pars_K_minus_1[1]
    lambda_k <- pars_K_minus_1[2]
    alpha_k  <- pars_K_minus_1[3]
    theta_k  <- NULL
    zeta_k   <- NULL
    b_k      <- 2 # pars[4]
  } else if(is.element(model, c("spEGARCH2"))) {
    rho    <- pars_K[1]
    lambda <- pars_K[2]
    alpha  <- pars_K[3]
    theta  <- 0.5 # pars[4]
    zeta   <- NULL # pars[5]
    b      <- NULL
    #
    rho_k    <- pars_K_minus_1[1]
    lambda_k <- pars_K_minus_1[2]
    alpha_k  <- pars_K_minus_1[3]
    theta_k  <- 0.5 # pars[4]
    zeta_k   <- NULL # pars[5]
    b_k      <- NULL
  }
  
  
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
  eps          <- as.vector(result_g[[1]])
  h            <- as.vector(result_g[[2]])
  
  z <- log(y^2) - (digamma(1) - log(2)) # log(y_i^2) - E(log(eps_i^2))
  h_tilde <- log(h) # + (alpha - alpha_k) * d_h_d_alpha + (rho - rho_k) * d_h_d_rho + (lambda - lambda_k) * d_h_d_lambda
  
  #squares <- (z - h_tilde)^2
  
  # return(list(eps = eps, y = y, h = h, z = z))
  return(z - h_tilde)
  
}


#setting i to index combination that yields highest LLV
i=c(max_index)


#Now getting all parameter estimates for each covariate combination of length 1:8 with highest LLV compared to covariate combinations with equal length
WM=weight_matrix(A_stand[,c(i)], knn)
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
  pars     <- c(0.3, 0.5, 0.5)
  LB       <- c(0, 0, 0)
  UB       <- c(1, 1, Inf)
} else if(is.element(model, c("spEGARCH"))) {
  pars     <- c(0.5, 0.5, 1)
  LB       <- c(0, 0, 0.00001)
  UB       <- c(Inf, Inf, Inf)
} else if(is.element(model, c("spEGARCH2"))) {
  pars     <- c(0.5, 0.5, 1, 0.5)
  LB       <- c(0, 0, 0.00001, 0.00001)
  UB       <- c(Inf, Inf, Inf, Inf)
}


#First
WM=weight_matrix(A_stand[,i[1]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
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
  
  if((rho+lambda)>=1){
    return(1e24)
  }
  J <- 0.5 * eps / sqrt(h) * d_h_d_eps(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) + diag(sqrt(h))
  
  # print(pars)
  
  log_det_J <- determinant(J, logarithm = TRUE)$modulus  
  
  return(as.numeric((-1) * (sum(log(dnorm(eps))) -  log_det_J)))
  
}
if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)

#Second estimates and standard errors
WM=weight_matrix(A_stand[,i[1:2]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)

#Third estimates and standard errors
WM=weight_matrix(A_stand[,i[1:3]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)

if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)
#Fourth estimates and standard errors
WM=weight_matrix(A_stand[,i[1:4]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)

#Fifth estimates and standard errors
WM=weight_matrix(A_stand[,i[1:5]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)

if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))


round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)


#Sixth estimates and standard errors
WM=weight_matrix(A_stand[,i[1:6]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)

#Seventh estimates and standard errors
WM=weight_matrix(A_stand[,i[1:7]], knn)
param    <- list(y = y, f_inv = functions$f_inv, tau_y = functions$tau_y, d_h_d_eps = functions$d_h_d_eps, g = functions$g, W_1 = WM, W_2 = WM, model = model)
if( is.na(LL(pars,param))){
  pars= c(0.45,0.45,2)
}
out      <- solnp(pars, fun = LL, LB = c(0, 0, 0), UB = c(1, 1, 100), param = param, control = list(trace = FALSE))

round(1/sqrt((diag(out$hessian))),3)

round(out$pars, 3)
round(sqrt(diag(solve(-out$hessian)%*%grad(LL,out$pars,param = param)%*%t(grad(LL,out$pars,param = param))%*%solve(-out$hessian))),3)






#Weight matrix plotting

WM=weight_matrix(A_stand[,c(15,38,43,50,53,60,65)],5)
boxplot(WM)
rownames(WM)=row.names(A_2019)
colnames(WM)=row.names(A_2019)

rownames(WM)[37]="Philip Morris"
colnames(WM)[37]="Philip Morris"

rownames(WM)[30]="Mitsubishi"
colnames(WM)[30]="Mitsubishi"
rownames(WM)[27]="JPMorgan Chase"
colnames(WM)[27]="JPMorgan Chase"

rownames(WM)[41]="Samsung"
colnames(WM)[41]="Samsung"

rownames(WM)[48]="Verizon Com"
colnames(WM)[48]="Verizon Com"

rownames(WM)[2]="Abbott Lab"
colnames(WM)[2]="Abbott Lab"

heatmap(WM, symm =T, cexRow = 0.55, cexCol = 0.55)

isSymmetric.matrix(WM)
matrix_obj = Matrix(WM)

names(matrix_obj)="s"

image(matrix_obj,axes=FALSE
      
      ,xlab=rownames(WM), ylab=rownames(WM))

plot(c(1,2),xlab=rownames(WM)[1:2],las=2)

image(Matrix(WM), sub =rownames(WM) )
axis(1, at = seq(0, 1, length = nrow(WM)), labels = rownames(WM))
axis(2, at = seq(0, 1, length = ncol(WM)), labels = colnames(WM))


