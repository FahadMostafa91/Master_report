
install.packages("pracma") # package for inverse of a mtx
library('pracma')

data2 <- read_csv("C:/Users/gmostafa/Downloads/MS analysis/data.csv")
X <- matrix(as.matrix(data2[,-4]),ncol=3, dimnames=NULL)
head(X)

#X[,2] <- log(X[,2])
#X[,3] <- log(X[,3])

plot(data2[,-1])
y = (data2$Weight)
n = 101 # sample size
xstar = cbind(1, log(45), log(23))
critical_value = 1.984
Sigma_hat = 0.1355
n = length(X[,1])
p = length(X[1,]) - 1
XpXinv = solve(t(X)%*%X) # Since we use this a lot, give it a name.
betahat = XpXinv %*% t(X) %*% y
MSE = 1/(n-(p+1)) * t(y - X%*%betahat) %*% (y - X%*%betahat)
MSE = as.numeric(MSE) # Make it a scalar.
ftcrit =1.984
#confidence interval 
lower =  as.numeric(xstar%*%betahat - ftcrit*sqrt(MSE)*sqrt(xstar%*%XpXinv%*%t(xstar)))
upper =  as.numeric(xstar%*%betahat + ftcrit*sqrt(MSE)*sqrt(xstar%*%XpXinv%*%t(xstar)))
df1= data.frame(t((lower)),t((upper)))
df1
### prediction
p = x_star%*%(inv(t(X)%*%X))%*%t(x_star)
lower = x_star%*%t(beta_hat) - critical_value*Sigma_hat*sqrt(1+p)
upper = x_star%*%t(beta_hat) + critical_value*Sigma_hat*sqrt(1+p)

lower
upper


