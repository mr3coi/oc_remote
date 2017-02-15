
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
	#plot(model)
	
	##### Variable selection based on regularized regression
	l 		= model$lambda.1se 		### use '1se' instead of 'min' to regularize as much as possible
	var.ind	= rep(1,ncol(input.x))
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
	var.ind[abs(coeffs) < thres] = 0
	
	if (all(abs(coeffs) < thres)) {			### regularized too much =>  try 'min'
		cat("Trying w/ 'min' instead of '1se'...")
		plot(model)
		l 		= model$lambda.min
		var.ind	= rep(1,ncol(input.x))
		coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
		var.ind[abs(coeffs) < thres] = 0
	}
	
	if (all(abs(coeffs) < thres)) {			### still regularized too much =>  try EN w/ alpha=0.8
		### Plot the previous LASSO and new EN models to see what lambda values are taken
		cat("Trying EN w/ alpha = 0.8")
		par(mfrow=c(2,1)); plot(model)
		model	= cv.glmnet(x.train, y.train, alpha=0.8,
					  family="binomial", nfolds=k, nlambda=50, type.measure=crit)
		plot(model)
		
		l 		= model$lambda.min
		var.ind	= rep(1,ncol(input.x))
		coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
		var.ind[abs(coeffs) < thres] = 0
	}
	
	return(list(lambda=l, vars=var.ind, coeffs=coeffs))
}