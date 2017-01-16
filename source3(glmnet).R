######################################################
### Packages
######################################################

library(glmnet)





######################################################
### Program
######################################################


##### Sample data generation 
data(iris)

preObj = preProcess(iris[,-5],method=c("center","scale"))
iris_input = predict(preObj,iris)
iris_input[,5] = as.factor(ifelse(iris_input[,5]=="setosa",1,0))
colnames(iris_input)[5] = "setosa"


##### Arguments definition
name		= "setosa"						# the name of the response variable
data		= iris_input					# preprocessed version of input
k			= 3								# fold number for CV
resp.ind	= which(colnames(data) == name)
x 			= as.matrix(data[,-resp.ind])	# matrix of explanatory variables
y 			= as.matrix(data[,resp.ind])	# vector of response variable
LASSO		= T 							# Choose b/w LASSO / Ridge
thres		= 10^(-5)						# threshold for discarding variables w/ too-small coefficients

##### Model generation
model	= cv.glmnet(x,y,alpha=ifelse(LASSO,1,0),family="binomial",nfolds=k,nlambda=50)
plot(model,xvar="lambda",label=T)

##### Variable selection based on regularized regression
coeffs	= as.vector(coef(model,s="lambda.min"))[-1]
#params = coef(model,s="lambda.1se")

var.ind	= ifelse(abs(coeffs) > thres,1,0)
