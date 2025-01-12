#PMF of ZOMP innovation
zompPMF <- function(k, phi, lam, p){
  if(k == 0){
    return( phi*(1-p) + (1-phi)*exp(-lam) )
  }else if(k == 1){
    return( phi*p + (1-phi)*lam*exp(-lam) )
  }else{
    return( (1-phi)*dpois(k,lam, log=FALSE) )
  }
}

#Simulate ZOMP random variable
rzomp <- function(phi, lam, p){
  ber <- rbinom(1, 1, phi)
  if(ber == 1){
    return(rbinom(1, 1, p))
  }else{
    return(rpois(1, lam))
  }
}


#One-step PMF for process
pij <- function(i, j, alpha, phi, lam, p){
  sum <- 0
  for(k in 0:min(i,j)){
    if(j-k==0){
      sum <- sum + dbinom(k, i, alpha)*zompPMF(0, phi, lam, p)
    }else if(j-k==1){
      sum <- sum + dbinom(k, i, alpha)*zompPMF(1, phi, lam, p)
    }else if(j-k>=2){
      sum <- sum + dbinom(k, i, alpha)*zompPMF((j-k), phi, lam, p)
    }
  }
  return(sum)
}

pij(5, 6,  0.3, 0.6, 2, 0.7)



#PMF of thinned innovation
inn_thin <- function(k, alpha, phi, lam, p){
  s <- 0
  for(i in 0:1e+4){
    s <- s + dbinom(k, i, alpha)*zompPMF(i, phi, lam, p)
  }
  return(s)
}




#Simulating INAR(1) process with ZOMP innovations
sim <- function(n, alpha, phi, lam, p){
  X <- c()
  X[1] <- rzomp(phi, lam, p)
  if(n>1){
    for(i in 2:n){
      X[i] <- rbinom(1, X[i-1], alpha) + rzomp(phi, lam, p)
    }
  }
  return(X)
}






#Simulate process
sample <- sim(500, 0.2, 0.8, 5, 0.6)

#Estimation by Conditional Maximum Likelihood approach
umat<- matrix(c(1,0,0,0,-1,0,0,0,0,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,-1), nrow=7, ncol=4, byrow=T)
cvec <- matrix(c(0,-1,0,-1,0,0,-1),nrow=7, ncol=1, byrow=T)
init <- c(0.5,0.5,10,0.5)
LLH_wrapper <- function(params){
  alpha <- params[1]
  phi <- params[2]
  lam <- params[3]
  p <- params[4]
  #Log-Likelihood function
  LLH <- function(alpha, phi, lam, p){
    llh <- log(zompPMF(sample[1], phi, lam, p))
    if(length(sample)>1){
      for(i in 2:length(sample)){
        llh <-  llh + log(pij(sample[i-1], sample[i], alpha, phi, lam, p))
      }
    }
    return(llh)
  }
  return(LLH(alpha, phi, lam, p))
}
CML <- constrOptim(init, LLH_wrapper, grad=NULL, ui=umat, ci=cvec, control=list(fnscale=-1))
CML$par


alpha <- c()
phi <- c()
lam <- c()
p <- c()
for(iter in 1:1){
  sample <- sim(1000, 0.8, 0.7, 10, 0.9)
  #Estimation by Conditional Least Squares approach
  umat<- matrix(c(1,0,0,0,-1,0,0,0,0,1,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,-1), nrow=7, ncol=4, byrow=T)
  cvec <- matrix(c(0,-1,0,-1,0,0,-1),nrow=7, ncol=1, byrow=T)
  init <- c(0.5,0.5,5,0.5)
  CLS_func <- function(alpha, phi, lam, p){
    Q <- c()
    Q[1] <- 0
    for(i in 2:length(sample)){
      Q[i] <- (sample[i] - alpha*sample[i-1] - phi*p - (1-phi)*lam)^2
    }
    return(sum(Q))
  }
  CLS_wrapper <- function(params){
    alpha <- params[1]
    phi <- params[2]
    lam <- params[3]
    p <- params[4]
    return(CLS_func(alpha, phi, lam, p))
  }
  CLS <- constrOptim(init, CLS_wrapper, grad=NULL, ui=umat, ci=cvec)
  alpha <- CLS$par[1]
  phi <- CLS$par[2]
  lam <- CLS$par[3]
  p <- CLS$par[4]
}
mean(alpha)
mean(phi)
mean(lam)
mean(p)







sample <- sim(500, 0.2, 0.2, 10, 0.5)
#Moment based estimation
num <- 0
for(i in 2:length(sample)){
  num <- num + (sample[i]-mean(sample[-1]))*(sample[i-1]-mean(sample[-1]))
}
denom <- 0
for(i in 2:length(sample)){
  denom <- denom + (sample[i]-mean(sample[-1]))^2
}


alpha_est <- num/denom
xbar <- mean(sample)
u2 <- sum(sample^2)/length(sample)
u3 <- sum(sample^3)/length(sample)

loss <- function(X) {
  phi_est = X[1]
  p_est = X[2]
  lambda_est = (xbar - (phi_est*p_est)/(1-phi_est))
  
  A1 <- alpha_est*xbar - (alpha_est^2)*xbar + (alpha_est^2)*u2 + phi_est*p_est
  B1 <- 2*alpha_est*(phi_est*p_est + (xbar^2)*(1-phi_est) - xbar*phi_est*p_est)
  C1 <- (xbar*(1-phi_est) - phi_est*p_est)^2
  D1 <- 4*phi_est*p_est*(xbar*(1-phi_est) - phi_est*p_est)
  E1 <- (xbar*(1-phi_est) - phi_est*p_est)*(1-phi_est*lambda_est)
  
  
  A2 <- (alpha_est^3)*(u3 - 3*u2 + 2*xbar) + 3*(alpha_est^2)*(u2-xbar) + alpha_est*xbar + (lambda_est^3)*(1-phi_est)
  B2 <- 3*(lambda_est^2)*(1-phi_est)^2 + 12*phi_est*p_est*lambda_est*(1-phi_est) + phi_est*p_est 
  C2 <- 3*lambda_est*(1-phi_est)*(1-phi_est*lambda_est) - 2*lambda_est*(1-phi_est)
  D2 <- 3*alpha_est*xbar*((lambda_est*(1-phi_est))^2 + 4*phi_est*lambda_est*p_est*(1-phi_est) + phi_est*p_est + lambda_est*(1-phi_est)*(1-phi_est*lambda_est))
  E2 <- 3*(phi_est*p_est + lambda_est*(1-phi_est))*(alpha_est*(alpha_est-1)*xbar - u2*alpha_est^2)
    
    
  EQ1 <- A1 + B1 + C1 + D1 + E1 - u2
  EQ2 <- A2 + B2 + C2 + D2 + E2 - u3
  return(sum(abs(EQ1),abs(EQ2)))
}

alpha_est
phi <- nlm(loss, c(0.5, 0.5))$estimate[1];phi
p <- nlm(loss, c(0.5,0.5))$estimate[2];p
lam <- xbar - (phi*p)/(1-phi);lam


#Monthly number of U.S. cases of poliomyelitis for 1970 to 1983; Zeger(1988)
sample <- c(0,1,0,0,1,3,9,2,3,5,3,5,2,2,0,1,0,1,3,3,
            2,1,1,5,0,3,1,0,1,4,0,0,1,6,14,1,1,0,0,1,
            1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,2,
            0,1,0,1,0,0,1,2,0,0,1,2,0,3,1,1,0,2,0,4,
            0,2,1,1,1,1,0,1,1,0,2,1,3,1,2,4,0,0,0,1,
            0,1,0,2,2,4,2,3,3,0,0,2,7,8,2,4,1,1,2,4,0,
            1,1,1,3,0,0,0,0,1,0,1,1,0,0,0,0,0,1,2,0,2,
            0,0,0,1,0,1,0,1,0,2,0,0,1,2,0,1,0,0,0,1,2,
            1,0,1,3,6)


p1 <- hist(sample, breaks=c(0:14))
p2 <- hist(sim(168, 0.13, 0.7, 2.9, 0.42))
plot(p1, col=rgb(0,0,1,1/5))  
plot(p2, col=rgb(1,0,0,1/5), add=T) 
estimates

