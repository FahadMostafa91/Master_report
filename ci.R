bayesreg = function(y, X, b0, B0, n0, S0, xstar=NA, alpha=0.05)
{
  # Perform Bayes regression.  NO ERROR CHECKING.
  n = length(X[,1])
  p = length(X[1,]) - 1
  
  # Find posterior parameters.
  B1 = solve(solve(B0)+t(X)%*%X)
  b1 = B1 %*% (solve(B0)%*%b0+t(X)%*%y)
  n1 = n0 + n
  S1 = as.numeric((n0*S0 + (t(b0)%*%solve(B0)%*%b0 + t(y)%*%y - t(b1)%*%solve(B1)%*%b1))/n1)
  
  # Find Bayesian intervals.
  lbci = NA
  ubci = NA
  lbpi = NA
  ubpi = NA
  if (length(xstar)>1)
  {
    btcrit = qt(1-alpha/2,df=n1)
    lbci = as.numeric(xstar%*%b1 - btcrit*sqrt(S1)*sqrt(xstar%*%B1%*%t(xstar)))
    ubci = as.numeric(xstar%*%b1 + btcrit*sqrt(S1)*sqrt(xstar%*%B1%*%t(xstar)))
    lbpi = as.numeric(xstar%*%b1 - btcrit*sqrt(S1)*sqrt(1+xstar%*%B1%*%t(xstar)))
    ubpi = as.numeric(xstar%*%b1 + btcrit*sqrt(S1)*sqrt(1+xstar%*%B1%*%t(xstar)))
  }
  
  ansout = list (b1 = b1, B1=B1, n1=n1, S1=S1, bayesCI=c(lbci,ubci), bayesPI=c(lbpi,ubpi))
  return(ansout)
}

freqreg = function(y, X, xstar=NA, alpha=0.05)
{
  # Perform classical regression.  NO ERROR CHECKING.
  n = length(X[,1])
  p = length(X[1,]) - 1
  
  # Find estimates of parameters.
  XpXinv = solve(t(X)%*%X) # Since we use this a lot, give it a name.
  betahat = XpXinv %*% t(X) %*% y
  MSE = 1/(n-(p+1)) * t(y - X%*%betahat) %*% (y - X%*%betahat)
  MSE = as.numeric(MSE) # Make it a scalar.
  
  # Find Classical intervals.
  lfci = NA
  ufci = NA
  lfpi = NA
  ufpi = NA
  if (length(xstar)>1)
  {
    ftcrit = qt(1-alpha/2,df=n-(p+1))
    lfci = as.numeric(xstar%*%betahat - ftcrit*sqrt(MSE)*sqrt(xstar%*%XpXinv%*%t(xstar)))
    ufci = as.numeric(xstar%*%betahat + ftcrit*sqrt(MSE)*sqrt(xstar%*%XpXinv%*%t(xstar)))
    lfpi = as.numeric(xstar%*%betahat - ftcrit*sqrt(MSE)*sqrt(1+xstar%*%XpXinv%*%t(xstar)))
    ufpi = as.numeric(xstar%*%betahat + ftcrit*sqrt(MSE)*sqrt(1+xstar%*%XpXinv%*%t(xstar)))
  }
  
  ansout = list (betahat = betahat, MSE=MSE, freqCI=c(lfci,ufci), freqPI=c(lfpi,ufpi))
  return(ansout)
}

# Read in the data.
data = read_csv("C:/Users/gmostafa/Downloads/MS analysis/data.csv")
n = length(data[,1])
X = as.matrix(data[,1:3])
y = as.matrix(data[,4])

# Take log transformations.
X[,2] = log(X[,2])
X[,3] = log(X[,3])
y = log(y)

b0 = matrix(c(0.01605,2,1),nrow=3,ncol=1)
#V0 = 256
#B0 = matrix(c(V0,0,0, 0,V0,0, 0,0,V0),nrow=3,ncol=3)
n0 = 10
S0 = 0.0194
xstar = matrix(c(1,log(45),log(23)),nrow=1,ncol=3)

V0 = c(1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256)
CI_lower = zeros(11,1)
CI_upper = zeros(11,1)
HPD_lower = zeros(11,1)
HPD_upper = zeros(11,1)
### loop
for(i in 1:11){
  V = V0[i]
  B0 <- matrix(c(V, 0.0, 0.0, 0.0, V, 0.0, 0.0, 0.0, V), 
               nrow = 3, ncol = 3)   # B0
  breg[i] = bayesreg(y,X,b0,B0,n0,S0,xstar)
  exp(breg$bayesCI)
  exp(breg$bayesPI)
}
freg = freqreg(y,X,xstar)
exp(freg$freqCI)
exp(freg$freqPI)

