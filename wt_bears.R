library(readxl)
library(caret)
library(e1071)
install.packages("MASS")
library(MASS)
library(lattice)
library(ggplot2)
install.packages("car")
library(car)
install.packages("lmtest")
library(lmtest)
install.packages("sandwich")
library(sandwich)
#https://www.statology.org/weighted-least-squares-in-r/

Bears1 <- read_csv("Bears1.csv")
#View(Bears1)
data <- Bears1
head(data)
# Plot the original data
plot(data)
# missing values visualization
#install.packages("naniar")
library(naniar)
vis_miss(data)

# data visualization
boxplot(data[,-2],col=c("red", "limegreen"))
plot(data$Headlen,data$Weight, xlab="Head Length", ylab="Weight")

#3#########################-------original model ------######################
df = na.omit(data2)
x1= ((df$Neckgirth))
x2=(df$Neckgirth)^(2)
y = (log((df$Weight)))
######## model ###############
model3 <- lm(y~ x1)
summary(model3)
#par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(model3)

### goodness of fit test ##############################
plot(x1,y)
abline(model3$coefficients)
## Compute the test statistic and the p-value
# v = Residual standard error x degrees of freedom
v = 0.2796*55
x = 1-pchisq(v,55)

###### normal MLR

x1= (df$Length)
x2=((df$Chestgirth)^2)  # 
y = (((df$Weight)))   # response 
model <- lm(y~ x1+x2)
summary(model)
anova(model)

#par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(model)
##########################################################################
## Frequentist linear Regression ################################
#Let us think the valume of a cylinder or a cone
# volume of cone = 1/3 * pi * r^2 *h
# volume of cylinder = pi * r^2 *h
# r = chest girth and h = length of bear
df = na.omit(data) # remove missing data
x1= log(df$Length)
x2=2*log((df$Chestgirth))  # 
y = (log((df$Weight)))   # response 
model <- lm(y~ x1+x2)
summary(model)
anova(model)
# confidence 
df =data.frame(x1,x2,y)
confint(model)
# bayesian cross validation
# install.packages("blmeco")
loo.cv(model, bias.corr=TRUE)
loo.cv(model, bias.corr=TRUE)  # increase nsim!!

########################
# High Leverage (hat) points
influenceIndexPlot(model, id.n=3)
#Influence bubble plot
influencePlot(model, id.n=3)
#####################################################
### cross validation
df1= data.frame(y,x1,x2)
## 75% of the sample size
smp_size <- floor(0.75 * nrow(df))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(df1)), size = smp_size)

train <- mtcars[train_ind, ]
test <- mtcars[-train_ind, ]
train3.control=trainControl(method="cv",number=5)
model3<- train(y~., data = df1, method = "lm", trControl = train3.control)
print(model3)
model3$results
#########################
par(mfrow = c(1, 2))
plot(model)
## goodness of fit test and p value 
#test_stat = Residual standard error * degrees of freedom
test_stat = 0.1393* 53
p_value = 1 - pchisq (test_stat, 53)
p_value
## watch plot to see it
plot(x1*x2,y)
abline(model$coefficients)
#par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(model)














