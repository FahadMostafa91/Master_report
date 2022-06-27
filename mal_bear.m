% log scale 
X = [data2(:,2),data2(:,3)];
y = data2(:,4);
mdl = fitlm(X,y);
disp(mdl)
anova(mdl,'summary')
plot(mdl)
weightPred = predict(mdl);
evaluateFit(y,weightPred,"Linear Model")
%% on original data data(:,2),
X = [data(:,2),data(:,3)];
y = data(:,4);
mdl = fitlm(X,y);
disp(mdl)
anova(mdl,'summary')
plot(mdl)
weightPred = predict(mdl);
evaluateFit(y,weightPred,"Linear Model")
% correlation plot 
dat=[X,y];
corrplot(dat)
%
plotAdded(mdl)
x1=data(:,2);
x2=data(:,3);
y=data(:,4);
surffit = fit([x1,x2],y,'poly23','normalize','on');
plot(surffit,[x1,x2],y,'Style','Residuals')

%% ML implementation 
dat=[X,y];
cv = cvpartition(size(dat,1),'HoldOut',0.3);
idx = cv.test;
% Separate to training and test data
dataTrain = dat(~idx,:);
dataTest  = dat(idx,:);

%Create a vector lambda of integers from 0 to 100. Then perform ridge regression on the training data for the values in lambda. 
%Return the coefficients in the scale of the original data, and name the coefficients b.
scatter3(dataTrain(:,1),dataTrain(:,2),dataTrain(:,3),"red")
scatter3(x1,x2,y)
lambda = 0:0.1:0.8;
b = ridge(dataTrain(:,3),dataTrain(:,1:2),lambda,0);
%Plot b against lambda to see how the coefficients change as λ increases.
%This plot is called the ridge trace.
plot(lambda,b)
%Predict the response for the test data XTest. Name the result yPred.
yPred = b(1,:) + dataTest(:,1:2)*b(2:end,:);
%Calculate mean squared error and name the result mdlMSE. Plot mdlMSE against lambda.
err = dataTest(:,3) - yPred;
mdlMSE = mean(err.^2);
plot(lambda,mdlMSE)
xlabel("\lambda")
ylabel("MSE")
%Find the smallest MSE and the index where it occurs. Name the results minMSE and idx, respectively.
[minMSE,idx] = min(mdlMSE);
%You can find the coefficients which minimize MSE by indexing into b.
bMin = b(:,idx);
%Try plotting the actual response with the predicted response.
plot(dataTest,"o",, LineWidth=2)
hold on
plot(yPred(:,idx),".", LineWidth=2)
hold off

plot(mdl)
weightPred = predict(mdl);
evaluateFit(y,dataTest(:,3),yPred,"Linear Model")



                       SumSq      DF       MeanSq        F         pValue  
                     _________    ___    __________    ______    __________

    Total               42.248    100       0.42648                        
    Model               40.703      2        20.352    1291.2    3.8754e-71
    Residual            1.5447     98      0.019400                        
    . Lack of fit        1.543     96      0.016073    19.187      0.050758
    . Pure error     0.0016754      2    0.00083771 
