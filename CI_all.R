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

plotci = function(lower, upper, names)
{
  numci = length(lower)
  y.plot = numci:1
  plot(10, 20, xlim = c(min(lower),max(upper)), ylim=c(1,numci), ylab="V0", xlab="Weight (lbs)", yaxt="n")
  for (i in 1:numci)
  {
    lines(c(lower[i],upper[i]), c(y.plot[i],y.plot[i]))
    lines(c(lower[i],lower[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
    lines(c(upper[i],upper[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
  }
  axis(side=2, at=y.plot, labels=names)
}


# Read in the data.
data =  read_csv("data.csv")
n = length(data[,1])
X = as.matrix(data[,1:3])
y = as.matrix(data[,4])

# Take log transformations.
X[,2] = log(X[,2])
X[,3] = log(X[,3])
y = log(y)

b0 = matrix(c(0.01605,2,1),nrow=3,ncol=1)
n0=50
S0 = 0.1355
xstar = matrix(c(1,log(45),log(23)),nrow=1,ncol=3)

V0 = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)
lowerci = vector(mode="numeric", length=length(V0))
upperci = vector(mode="numeric", length=length(V0))
lowerpi = vector(mode="numeric", length=length(V0))
upperpi = vector(mode="numeric", length=length(V0))
for (i in 1:length(V0))
{
  B0 = matrix(c(V0[i],0,0, 0,V0[i],0, 0,0,V0[i]),nrow=3,ncol=3)
  breg = bayesreg(y,X,b0,B0,n0,S0,xstar)
  lowerci[i] = exp(breg$bayesCI[1])
  upperci[i] = exp(breg$bayesCI[2])
  lowerpi[i] = exp(breg$bayesPI[1])
  upperpi[i] = exp(breg$bayesPI[2])
}

freg = freqreg(y,X,xstar)
lowerci = c(lowerci, exp(freg$freqCI[1]))
upperci = c(upperci, exp(freg$freqCI[2]))
lowerpi = c(lowerpi, exp(freg$freqPI[1]))
upperpi = c(upperpi, exp(freg$freqPI[2]))

ynames = c(as.character(V0), "Freq")

plotci(lower=lowerci, upper=upperci, names=ynames)
plotci(lower=lowerpi, upper=upperpi, names=ynames)

