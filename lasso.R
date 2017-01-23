
##### ============================================================================================
##### ============================================================================================
##### ============== Variable Selection : Regularized Regression (LASSO) =============
##### ============================================================================================
##### ============================================================================================


##### Returns variables selected by regularized regression (LASSO/Ridge) through CV.
##### It is assumed that all rows w/ NA's are removed via preprocessing beforehand.
# input		: input data (df)
# resp.ind	: column index of response variable
# k			: fold number for CV
# LASSO		: Choose b/w LASSO(T) and Ridge(F)
# thres		: threshold for discarding variables w/ too-small coefficients
# crit		: criterion for finding the best choice of 'lambda'
regular.CV = function(input, resp.ind, k=3, LASSO=T, thres=10^(-5), crit="deviance") {
	##### Partitioning
	input.x	= data.matrix(input[,-resp.ind])
	input.y	= data.matrix(input[,resp.ind])
	
	inTrain	= createDataPartition(input.y,p=0.7,list=F)
	x.train	= input.x[inTrain,]
	y.train	= input.y[inTrain,]
	x.test	= input.x[-inTrain,]
	y.test	= input.y[-inTrain,]
	
	##### Model generation
	model	= cv.glmnet(x.train, y.train, alpha=ifelse(LASSO,1,0),
					  family="binomial", nfolds=k, nlambda=50, type.measure=crit)
	
	##### Plot the generated model
	#plot(model,xvar="lambda",label=T)
	
	##### Variable selection based on regularized regression
	l 		= model$lambda.min
	var.ind	= rep(0,ncol(input.x))
	var.ind[unlist(predict(model, newx=x.test, type="nonzero", s=l))] = 1
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
	
	return(list(lambda=l, vars=var.ind, coeffs=coeffs))
}