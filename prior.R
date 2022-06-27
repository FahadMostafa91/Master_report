
install.packages("pracma") # package for inverse of a mtx
library('pracma')
# data input
data <- read_csv("data.csv") # data set with design matrix
### find prediction and HPD interval
n0 = 10; # 10, 25, 50 ; subject to change 
S0 = 0.0194 # inch^3, sigma^2, used from frequentist MSE
V0 = c(1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256)
V = 2
shape = n0/2
scale = n0*S0/2

b0 = c(16.03*10^(-3), 2, 1) #b0
B0 <- matrix(c(V, 0.0, 0.0, 0.0, V, 0.0, 0.0, 0.0, V), 
            nrow = 3, ncol = 3)   # B0
head(data)
X <- matrix(as.matrix(data[,-4]),ncol=3, dimnames=NULL)
X
head(X)
plot(data[,-1])
y = (data$Weight)
# posterior parameters are

#A = inv(B0)+t(X)%*%X
n = 101
B1 = inv(inv(B0)+t(X)%*%X)
b1 = B1%*%(inv(B0)%*%b0+t(X)%*%y)
n1 = n0 + n
n1S1 = n0*S0 + (t(b0)%*%inv(B0)%*%b0 + t(y)%*%y - t(b1)%*%inv(B1)%*%b1)
S1 = n1S1/n1
# Finally we will get the posterior parameter 
b1
B1
n1/2
n1S1/2
# The posterior predictive distribution 


# (1-alpha)*100% prediction interval Ynew|y at x* in Bayesian setting 
b1
x_star = cbind(1, 45, 23)
critical_value = 1.984
#MSE = S1
x_star%*%b1 - critical_value*sqrt(S1%*%(1+x_star%*%B1%*%(t(x_star)))) # lower
x_star%*%b1 + critical_value*sqrt(S1%*%(1+x_star%*%B1%*%(t(x_star)))) # upper
# (1-alpha)*100% HPD
x_star%*%b1 - critical_value*sqrt(S1%*%(x_star%*%B1%*%(t(x_star))))  # lower
x_star%*%b1 + critical_value*sqrt(S1%*%(x_star%*%B1%*%(t(x_star))))  # upper
#---------------------------------------------------------------------


# draw confidence interval in R

require(readxl)
cis = read_excel(path="CI.xlsx")

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

plotci(lower=cis$lower, upper=cis$upper, names=cis$V0)

#############################################################
#-------------    data visualization -------------------#####
# whole data set by scatter plot 
plot(data[,-1])
# summary statistics using box plots 
Length <- data$Length
Chestgirth <- data$Chestgirth
df1 = data.frame(Length,Chestgirth)
colors = c( rep("blue",1), rep("red",1))
boxplot(df1, col = colors, ylab="inches")
# response variables 
Weight <- data$Weight
colors = c( rep("green",1))
boxplot(Weight, col = colors, ylab="lbs", xlab="Weight")
################################################################
lnr <- log(data$Chestgirth)
lnL <- log(data$Length)
lnW <- log(data$Weight)
par(mfrow = c(1, 3))
plot(lnW,lnr,cex.lab = 1.3,cex.axis = 1.2,lwd=1.2)
plot(lnW,lnL, col="blue",cex.lab = 1.3,cex.axis = 1.2,lwd=1.2)
plot(lnW,2*lnr+lnL, col="red",cex.lab = 1.3,cex.axis = 1.2,lwd=1.2)



#lm model

model=lm(log(Weight)~log(Length)+log(Chestgirth), data = data )
summary(model)
anova(model)
#install.packages("performance")
#install.packages("see")
#library(performance)
check_model(model)

library(tidyverse)
library(broom)
library(glmnet)

X = data%>%select(Chestgirth,Length)%>% data.matrix()
y = data$Weight
lambdas <- 10^seq(3, 0, by = .1)
cv_fit <- cv.glmnet(X, y, alpha = 0, lambda = lambdas)

plot(cv_fit)

opt_lambda <- cv_fit$lambda.min
opt_lambda
fit <- glmnet(X, y, alpha = 0, lambda = lambdas)

y_predicted <- predict(fit, s = opt_lambda, newx = X)

# Sum of Squares Total and Error
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

# R squared
rsq <- 1 - sse / sst
rsq

#find coefficients of best model
best_model <- glmnet(X, y, alpha = 0, lambda = opt_lambda)
coef(best_model)
#produce Ridge trace plot
model <- glmnet(X, y, alpha = 0)
plot(model, xvar = "lambda")

y_predicted <- predict(model, s = opt_lambda, newx = X)

#find SST and SSE
sst <- sum((y - mean(y))^2)
sse <- sum((y_predicted - y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq

mdl=lm(Weight ~Chestgirth+Length, data=data)
summary(mdl)

df = data.frame(Chestgirth,Length,Weight)
summary(df)


