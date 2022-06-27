
install.packages("pracma") # package for inverse of a mtx
library('pracma')
# data input
#data <- read_csv("data.csv") # data set with design matrix
#head(data)
#data2 <- log(data) # data with log transform
#write.csv(data2, "/Users/Fahad/Downloads/data2.csv", row.names = FALSE)

data2 <- read_csv("data2.csv")
X <- matrix(as.matrix(data2[,-4]),ncol=3, dimnames=NULL)
head(X)
#X[,2] <- log(X[,2])
#X[,3] <- log(X[,3])

plot(data2[,-1])
y = (data2$Weight)
n = 101 # sample size
n0 = 10; # 10, 25, 50 ; subject to change 
S0 = 0.1355 #0.0194 # inch^3, sigma^2, used from frequentist MSE
shape = n0/2
scale = n0*S0/2
b0 = c(16.03*10^(-3), 2, 1) #b0
n1 = n0 + n
x_star = cbind(1, log(45), log(23))
critical_value = 1.984


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
  # posterior parameters are
  #A = inv(B0)+t(X)%*%X
  B1 = inv(inv(B0)+t(X)%*%X)
  b1 = B1%*%(inv(B0)%*%b0+t(X)%*%y)
  n1S1 = n0*S0 + (t(b0)%*%inv(B0)%*%b0 + t(y)%*%y - t(b1)%*%inv(B1)%*%b1)
  S1 = n1S1/n1
  # (1-alpha)*100% prediction interval #----------------------------
  CI_lower[i] = (x_star%*%b1 - critical_value*sqrt(S1%*%(1+x_star%*%B1%*%(t(x_star))))) # lower
  CI_upper[i] = (x_star%*%b1 + critical_value*sqrt(S1%*%(1+x_star%*%B1%*%(t(x_star))))) # upper
  # (1-alpha)*100% HPD Interval #------------------------------------
  HPD_lower[i] = (x_star%*%b1 - critical_value*sqrt(S1%*%(x_star%*%B1%*%(t(x_star)))))  # lower
  HPD_upper[i] = (x_star%*%b1 + critical_value*sqrt(S1%*%(x_star%*%B1%*%(t(x_star)))))  # upper
  ##########-------------------------------------------------------------------
}

# prediction interval
CI = data.frame(CI_lower, CI_upper)
CI

# draw confidence interval in R

CI = data.frame(V0,CI_lower,CI_upper)
cis = CI

plotci = function(lower, upper, names)
{
  numci = length(lower)
  y.plot = numci:1
  plot(10, 20, xlim = c(min(lower),max(upper)), ylim=c(1,numci), ylab="", xlab="Weight (lbs)", yaxt="n")
  for (i in 1:numci)
  {
    lines(c(lower[i],upper[i]), c(y.plot[i],y.plot[i]))
    lines(c(lower[i],lower[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
    lines(c(upper[i],upper[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
  }
  axis(side=2, at=y.plot, labels=names)
}

#plotci(lower=cis$CI_lower, upper=cis$CI_upper, names=cis$V0)

# conversion to original data 

L = exp(cis$CI_lower)
U = exp(cis$CI_upper)
plotci(lower=L, upper=U, names=cis$V0)
########
HPD = data.frame(HPD_lower[], HPD_upper[])
HPD

CI = data.frame(V0,HPD_lower,HPD_upper)
cis = CI

plotci = function(lower, upper, names)
{
  numci = length(lower)
  y.plot = numci:1
  plot(10, 20, xlim = c(min(lower),max(upper)), ylim=c(1,numci), ylab="", xlab="Weight (lbs)", yaxt="n")
  for (i in 1:numci)
  {
    lines(c(lower[i],upper[i]), c(y.plot[i],y.plot[i]))
    lines(c(lower[i],lower[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
    lines(c(upper[i],upper[i]), c(y.plot[i]-0.1, y.plot[i]+0.1))
  }
  axis(side=2, at=y.plot, labels=names)
}

#plotci(lower=cis$HPD_lower, upper=cis$HPD_upper, names=cis$V0)

# conversion to original data 
CL = exp(cis$HPD_lower)
CU = exp(cis$HPD_upper)
plotci(lower=CL, upper=CU, names=cis$V0)

