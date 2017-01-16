######################################################
### Packages & Settings
######################################################

library(glmnet)
library(caret)

##### Load previously saved 'ova.db.csv'(preprocessed) from 'datafile', and 'vars' from 'varfile'
wd <- "/Users/SJC/Documents/practice/internship"
setwd(wd)
load("ova.db.csv.RData")





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
regular.CV = function(input.x, input.y, resp.ind, k=3, LASSO=T, thres=10^(-5), crit="deviance") {
	##### Preprocessing
		### Remove rows w/ NA
	row_NA	= apply(input,1,function(v) {sum(is.na(v)) == 0})
	input	= input[row_NA,]
	
		### Separate training and testing datasets
	inTrain		= createDataPartition(y,p=0.7,list=F)
	x.train		= input.x[inTrain,]
	y.train		= input.y[inTrain,]
	x.test		= input.x[-inTrain,]
	y.test		= input.y[-inTrain,]
	
	##### Model generation
	model	= cv.glmnet(x.train, y.train, alpha=ifelse(LASSO,1,0),
					  family="binomial", nfolds=k, nlambda=50,
					  type.measure=crit)
	#plot(model,xvar="lambda",label=T)
	
	##### Variable selection based on regularized regression
	l 		= model$lambda.min
	var.ind	= rep(0,ncol(x))
	var.ind[unlist(predict(model, newx=x.test, type="nonzero", s=l))] = 1
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
	
	return(list(lambda=l, vars=var.ind, coeffs=coeffs))
}






######################################################
### Program
######################################################


# ##### Sample data generation 
# data(iris)
# 
# preObj = preProcess(iris[,-5],method=c("center","scale"))
# iris_input = predict(preObj,iris)
# iris_input[,5] = as.factor(ifelse(iris_input[,5]=="setosa",1,0))
# colnames(iris_input)[5] = "setosa"


##### Arguments definition
outcomes = c("Recurrence",
			 "Platinum_resistance_6mo", "Platinum-resistance_group",
			 "End_Date_1st_regimen", "Op_date",
			 "Dx_date", "Death", "Last_FU_Date", "Expired_Date",
			 "Residual_tumor_site_1st_debulking", "Residual_tumor_size_1st_debulking",
			 "PLN_status", "PALN_status")
name		= outcomes[1]					# the name of the response variable
input		= ova.db.csv					# preprocessed version of input
resp.ind	= which(colnames(input) == name)

##### Preprocess : remove NA's
row_NA	= apply(input,1,function(v) {sum(is.na(v)) == 0})
input	= input[row_NA,]

##### Partition
	# Avoid all columns corresponding to post-surgery
avoid = search_colname(ova.db.csv,"Start_Date_1st_regimen"):length(colnames(ova.db.csv))
	# Add other pre/during-surgery columns that are obviously meaningless or irrelevant
avoid = c(avoid, 1:3, 5, 43)
x 	= as.matrix(input[,-c(avoid,resp.ind)])	# matrix of explanatory variables
y 	= as.matrix(input[,resp.ind])	# vector of response variable

##### Function call
result_auc	= regular.CV(x,y, resp.ind, crit="auc"); result_auc
result_dev	= regular.CV(x,y, resp.ind, crit="dev"); result_dev
result_cls	= regular.CV(x,y, resp.ind, crit="class"); result_cls





######################################################
### Utility Functions
######################################################

search_colname = function(df, name) { return(grep(name, colnames(df))) }