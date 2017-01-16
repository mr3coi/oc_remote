######################################################
### Packages
######################################################

library(glmnet)
library(caret)



######################################################
### Functions
######################################################

##### Returns variables selected by regularized regression (LASSO/Ridge) through CV
# input		: input data (df)
# resp.ind	: column index of response variable
# k			: fold number for CV
# LASSO		: Choose b/w LASSO(T) and Ridge(F)
# thres		: threshold for discarding variables w/ too-small coefficients
# crit		: criterion for finding the best choice of 'lambda'
regular.CV = function(input, resp.ind, k=3, LASSO=T, thres=10^(-5), crit="deviance") {
	##### Preprocessing
		### Remove rows w/ NA
	row_NA	= apply(input,1,function(v) {sum(is.na(v)) == 0})
	input	= input[row_NA,]
	
		### Separate training and testing datasets
	inTrain		= createDataPartition(y,p=0.7,list=F)
	x 			= as.matrix(input[,-resp.ind])	# matrix of explanatory variables
	y 			= as.matrix(input[,resp.ind])	# vector of response variable
	x.train		= x[inTrain,]
	y.train		= y[inTrain,]
	x.test		= x[-inTrain,]
	y.test		= y[-inTrain,]
	
	##### Model generation
	model	= cv.glmnet(x.train, y.train, alpha=ifelse(LASSO,1,0),
					  family="binomial", nfolds=k, nlambda=50,
					  type.measure=crit)
	plot(model,xvar="lambda",label=T)
	
	##### Variable selection based on regularized regression
	var.ind	= rep(0,ncol(x))
	var.ind[unlist(predict(model, newx=x.test, type="nonzero", s=c("lambda.min")))] = 1
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s="lambda.min"))[-1]
	
	return(list(var.ind,coeffs))
}


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
input		= iris_input					# preprocessed version of input
resp.ind	= which(colnames(input) == name)

result		= regular.CV(input, resp.ind, crit="auc")