#Function to evaluate higher order derivatives
DD<-function(expr,name,order=1){
  if(order<1) stop("Order must be >=1")
  if(order==1)D(expr,name)
  else DD(D(expr,name),name,order-1)
}

# 5-step PGF
five_pgf <- expression( ((1- alpha^5 + (alpha^5)*s)^i) *
                     (phi - (s-1)*phi*p*alpha^(4) + (1-phi)*exp((s-1)*lambda*alpha^(4))) *
                     (phi - (s-1)*phi*p*alpha^(3) + (1-phi)*exp((s-1)*lambda*alpha^(3))) *
                     (phi - (s-1)*phi*p*alpha^(2) + (1-phi)*exp((s-1)*lambda*alpha^(2))) *
                     (phi - (s-1)*phi*p*alpha^(1) + (1-phi)*exp((s-1)*lambda*alpha^(1))) *
                     (phi*(1-p) + s*phi*p + (1-phi)*exp(lambda*(s-1))) )

# 5-step PMF
five_PMF <- function(i, j, alpha, phi, lambda, p){
  pgf <- five_pgf
    if(j==0){
    s=0
    return(eval(pgf))
  }else{
    diff <- DD(pgf,'s',j)
    s=0
    i=i
    alpha=alpha
    phi=phi
    lambda=lambda
    p=p
    return(eval(diff)/factorial(j))
  }
}
five_PMF(2,3,0.1,0.3,4,0.7)

#k-step PGF, defined only for k>1
k_step_PGF <- function(m){
  p1 <- paste0('((1- alpha^',m,' + (alpha^',m,')*s)^i)*') #Part A as annotated in PDF file
  p3 <- '(phi*(1-p) + s*phi*p + (1-phi)*exp(lambda*(s-1)))' #Part C
  str <- ''
  for(i in (m-1):1){ #For Part B, see page 37 on pdf, see highlighted area.
    str <- paste0(str, paste0('(phi - (s-1)*phi*p*alpha^(',i,') + (1-phi)*exp((s-1)*lambda*alpha^(',i,')))*'))
  }
  char <- paste0(p1,str,p3)
  return(parse(text=char)) #This is an expression object which can be differentiated(The function DD works similarly to the one in MATLAB)
}

# DD(k_step_PGF(6), 's', 1)

# k-step PMF
k_step_PMF <- function(m, i, j, alpha, phi, lambda, p){
  if(m==1){
    print('Error: PMF only defined for m>1')
  }else{
    pgf <- k_step_PGF(m)
    if(j==0){
      s=0
      return(eval(pgf))
    }else{
      diff <- DD(pgf,'s',j)
      s=0
      i=i
      alpha=alpha
      phi=phi
      lambda=lambda
      p=p
      # return(eval(diff))
      return(eval(diff)/factorial(j))
    }
  }
}


k_step_PMF(3, 4 ,6 , 0.5, 0.8, 10, 0.1) #P(Xt+3 = 6 | Xt = 4) for alpha=0.5, phi=0.8, lambda=10, p=0.1



pij <- c()
for(j in 0:5){
  pij[j+1] <- k_step_PMF(3, 4 ,j , 0.9, 0.8, 10, 0.8)
}
barplot(pij, names.arg=c(0:5), ylab='Probability', xlab='X')
sum(pij)
#Observe sum tends to 1, but computations are slow since for each X=i, pgf needs to be differentiated i times
#So sum is restricted to X=c(0:10)


forecast <- c()
for(j in 0:5){
  forecast[j+1] <- k_step_PMF(6,6,j,0.13,0.7,2.9,0.42)
}
plot(x=c(0:5),y=forecast,type='b')
plot(x=(which(forecast==max(forecast))-1), y=forecast[which(forecast==max(forecast))], col='red',pch=19)
points(x=(which(forecast==max(forecast))-1), y=forecast[which(forecast==max(forecast))], col='blue',pch=19)
points(x=(which(forecast==max(forecast))-1), y=forecast[which(forecast==max(forecast))], col='green',pch=19)
points(x=(which(forecast==max(forecast))-1), y=forecast[which(forecast==max(forecast))], col='purple',pch=19)
points(x=(which(forecast==max(forecast))-1), y=forecast[which(forecast==max(forecast))], col='black',pch=19)

