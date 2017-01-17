
##### ============================================================================================
##### ============================================================================================
##### =============== Packages & Initial Settings =================
##### ============================================================================================
##### ============================================================================================

library(e1071)
library(MASS)
library(ROCR)
library(nnet)
library(glmnet)
library(caret)
library(AUC)

##### Locate & read source data
wd 		 <- "/Users/SJC/Documents/practice/internship"
datafile <- "ova_db.csv" 
varfile  <- "ova_variable.csv" 

setwd(wd) 
ova.db.csv = as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 


##### ============================================================================================
##### ============================================================================================
##### =============== Utility Functions ===============
##### ============================================================================================
##### ============================================================================================

##### Date-to-int conversion function 				##### TODO convert to # of days since 1900-01-01
date_conv = function(date_string) {
	date = as.vector(sapply(strsplit(date_string,split="\\."),strtoi))
	date = sum(date*c(10000,100,1))
	return(date)
}

##### Returns the column index in 'df' of the variable w/ specified 'name'
search_colname = function(df, name) { return(grep(name, colnames(df))) }


##### ============================================================================================
##### ============================================================================================
##### ============ datafile : Preprocess ==============
##### ============================================================================================
##### ============================================================================================

##### Change 'unknown' values to NA
NAval <- sapply(strsplit(colnames(ova.db.csv),"\\n"), 
                function(v) { q  <- trimws(v[-1])           # Generate character vector of items in colnames
                                                            # w/ whitespaces removed
                qq <- grep("unknown", q, ignore.case=T)     # Find all indices in list w/ unknown items
                ifelse(length(qq),strsplit(q[qq[1]], "\\.")[[1]][1],NA) } ) 
for (i in 1:ncol(ova.db.csv)) {
  if (!is.na(NAval[i])) {
    ova.db.csv[which(NAval[i] == ova.db.csv[,i]),i] <- NA
  }
}

##### Specific deletions (flaws)
	### Delete unnamed columns
idx = which(colnames(ova.db.csv)=="")
	### Remove last two rows (empty rows)
ova.db.csv = ova.db.csv[ 1:236, -idx ]
	### Delete the repeated patient number column
ova.db.csv = ova.db.csv[,-which(colnames(ova.db.csv)=="Pt_No.1")]

##### Alter column names to main titles (removing choices) and remove spaces within the main titles
colnames(ova.db.csv) = gsub(" ", "_", sapply(strsplit(colnames(ova.db.csv),"\\n"),function(v)trimws(v[1]))) 

##### Find indices of columns w/ (almost) uniform values and remove such columns				### TODO almost
idx <- which( apply(ova.db.csv, 2, function(v) length(table(v)) == 1) )
ova.db.csv <- ova.db.csv[,-idx] 

##### Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
options(warn=2);	# to exploit errors

##### Find and remove/modify columns w/ non-numeric data
	### idx : non-numeric indices
idx <- which(apply(ova.db.csv, 2, function(v) class(try(as.numeric(v), T))=='try-error'))
options(warn=1);

	##### Convert date info into numeric data (YYYYMMDD)
pattern = "\\d+\\.\\d+\\.\\d+"
idx2 = c()

for (col in idx) {
	if (sum(is.na(ova.db.csv[,col])) == 0 && length(grep(pattern,ova.db.csv[,col])) != 0) {		### TODO Optimize
		ova.db.csv[,col] = sapply(ova.db.csv[,col], function(str) {ifelse(str != "", date_conv(str), NA)} )
		idx2 = c(idx2,col)
	}
}
idx = idx[!idx %in% idx2]

	##### Additional manipulation for 'Stage' column
ref = c(rep(1,3),rep(2,3),rep(3,3),4)
names(ref) = c("1a","1b","1c","2a","2b","2c","3a","3b","3c","4")
for (i in 1:nrow(ova.db.csv)) {
	ova.db.csv[i,"Stage"] = ref[ova.db.csv[i,"Stage"]]
}
idx = idx[names(idx)!="Stage"]

	##### Additional manipulation for 'Recur_Site' column
ova.db.csv[,84] = sapply(ova.db.csv[,84],function(v) {paste(sort(unlist(strsplit(v,split=","))),collapse="")})
idx = idx[names(idx)!="Recur_Site"]

	##### Extract numeric values at the beginning
pattern = "(\\d).*"
idx2 = c()
for (col in idx) {
	if (length(grep(pattern,ova.db.csv[,col])) != 0) {				### TODO Optimize
		ova.db.csv[,col] = sub(pattern,"\\1",ova.db.csv[,col])
		idx2 = c(idx2,col)
	}
}
idx = idx[!idx %in% idx2]

	### Delete the remaining non-numeric columns
ova.db.csv <- ova.db.csv[,-idx]

##### Convert character values in ova.db.csv into numeric values
ova.db.csv <- as.data.frame(apply(ova.db.csv, 2, as.numeric))

##### Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
# options(warn=2)

##### Check preprocessing results
apply(ova.db.csv,2,function(v) class(v))

##### Convert 'Refractory progression during CTx' to 'Yes'
ova.db.csv$Recurrence[ova.db.csv$Recurrence==2] <- 1


##### ============================================================================================
##### ============== varfile : Preprocess =============
##### ============================================================================================
##### ============================================================================================

# Read varfile
ova.var.csv = read.csv(varfile)

# Extract variable names from specific rows in varfile into a matrix and remove choices
vars <- apply(ova.var.csv[seq(2,18,4),], 1, function(v)sapply(strsplit(v,"\\n"),
                                            function(vv)trimws(vv[1]))
               )[-1,] 
rownames(vars) <- NULL
colnames(vars) = c("Recurrence","Platinum_resistance","PFS","OS","Op. Debulking")


##### ============================================================================================
##### ============================================================================================
##### ============== Variables & Arguments Setting =============
##### ============================================================================================
##### ============================================================================================

##### Interested response variables
outcomes = c("Recurrence",
			  "Platinum_resistance_6mo", "Platinum-resistance_group",
			  "End_Date_1st_regimen", "Op_date",
			  "Dx_date", "Death", "Last_FU_Date", "Expired_Date",
			  "Residual_tumor_site_1st_debulking", "Residual_tumor_size_1st_debulking",
			  "PLN_status", "PALN_status")
##### Reference
outcome_index = c("CR",
				  "CS", "CT",
				  "CO", "BF",
				  "E", "DC", "DD", "DE",
				  "CI", "CJ",
				  "AX", "AY")
outcome_idx = list(1,2:3,4:5,6:9,10:11,12:13)

##### Choose which set of response variable and explanatory variables to analyze
set_num 	= 1
var_num 	= 1
resp.var 	= outcomes[outcome_idx[[set_num]]]
resp.var 	= resp.var[resp.var %in% colnames(ova.db.csv)][var_num]
resp.ind	= which(colnames(ova.db.csv) == resp.var)

##### Extract relevant sub-data.frame from given data	=>	'ova.db.csv.sub'
	### Avoid all columns corresponding to post-surgery (specify indices of variables to avoid here)
avoid 		= search_colname(ova.db.csv,"Start_Date_1st_regimen"):length(colnames(ova.db.csv))
	### Add other pre/during-surgery columns that are obviously meaningless or irrelevant
avoid 		= c(avoid, 1:3, 5, 43)
avoid 		= avoid[avoid != resp.ind]		# keep response variable if included in 'avoid'
ova.db.csv.sub = ova.db.csv[,-avoid]		# exclude all variables w/ column indices in 'avoid'
resp.ind	= which(colnames(ova.db.csv.sub) == resp.var)	# update 'resp.ind'

##### Preprocess : remove NA's	=>	'input' (final data form)
row_NA		= apply(ova.db.csv.sub,1,function(v) {sum(is.na(v)) == 0})
input 		= ova.db.csv.sub[row_NA,]

##### Convert the column for response variable into a factor
input[,resp.ind] = as.factor(input[,resp.ind])



##### ============================================================================================
##### ============================================================================================
##### ============== Variable Selection : Stepwise AIC =============
##### ============================================================================================
##### ============================================================================================


stepwiseAIC = function(input,resp.var) {
	##### Preprocessing	
	exp.vars = colnames(input)[colnames(input)!=resp.var]
	
	##### Generate model using logistic regression
	f 	= as.formula(paste(resp.var, "~", paste(exp.vars, collapse="+")))
	mod = glm(f, data = input, family = "binomial")
	nul = glm(paste(resp.var, "~1"), data = input, family = "binomial")
	
	##### Conduct stepwise AIC both forwards and backwards
	step.0	<- stepAIC(nul, scope=list(lower=nul,upper=mod), direction="both")
	step.f	<- stepAIC(mod, direction="both")
	#step.0x <- step(nul, scope=list(lower=nul,upper=mod), direction="both")
	#step.fx <- step(mod, direction="both")
	
	# ##### Show ANOVA tables from the above regressions
	# step.0$anova
	# step.f$anova
	
	step0.coeff = ifelse(exp.vars %in% names(step.0$coefficients),1,0)
	stepF.coeff = ifelse(exp.vars %in% names(step.f$coefficients),1,0)
	
	##### Picked
	# var.full <- setdiff(exp_vars[[set_num]], substr(as.character(step.f$anova$Step), 3, 1000))
	# mean(do.cv(var.full, 3)$auc)
	# var.full <- substr(as.character(step.0$anova$Step), 3, 1000)[-1]
	# mean(do.cv(var.full, 3)$auc)
	# 
	# var.full <- setdiff(exp_vars[[set_num]], substr(as.character(step.fx$anova$Step), 3, 1000))
	# mean(do.cv(var.full, 3)$auc)
	# var.full <- substr(as.character(step.0x$anova$Step), 3, 1000)[-1]
	# mean(do.cv(var.full, 3)$auc)
	
	return(list(step0=step0.coeff,stepF=stepF.coeff))
}


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
	input.x	= as.matrix(input[,-resp.ind])
	input.y	= as.matrix(input[,resp.ind])
	
	inTrain	= createDataPartition(input.y,p=0.7,list=F)
	x.train	= input.x[inTrain,]
	y.train	= input.y[inTrain,]
	x.test	= input.x[-inTrain,]
	y.test	= input.y[-inTrain,]
	
	##### Model generation
	model	= cv.glmnet(x.train, y.train, alpha=ifelse(LASSO,1,0),
					  family="binomial", nfolds=k, nlambda=50,
					  type.measure=crit)
	
	##### Plot the generated model
	#plot(model,xvar="lambda",label=T)
	
	##### Variable selection based on regularized regression
	l 		= model$lambda.min
	var.ind	= rep(0,ncol(input.x))
	var.ind[unlist(predict(model, newx=x.test, type="nonzero", s=l))] = 1
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
	
	return(list(lambda=l, vars=var.ind, coeffs=coeffs))
}


##### ============================================================================================
##### ============================================================================================
##### ============== Performance evaluation : RF, kNN, SVM, Logistic Reg. ===============
##### ============================================================================================
##### ============================================================================================


######################################################
### fitting functions
######################################################

##### Wrapper for RF algorithm
# train.dat : dataset for train.dat the machine
# test.dat	   : dataset for test.dat the machine
# resp.ind : index of the response variable in 'train.dat' and 'test.dat' datasets
# marker   : vector indicating which columns to use (1) and not to use (0) as
#				explanatory variables in 'train.dat' and 'test.dat' datasets
doRF = function(train.dat, test.dat, resp.ind, marker) {
	require(caret)
	resp 	= colnames(train.dat)[resp.ind]
	vars	= colnames(train.dat)
	f 		= as.formula( paste(resp, "~", paste(vars[which(marker==1)],collapse="+")) )
	model 	= train(f, data=train.dat, method="rf")
#	model 	= train(train.dat[,-resp.ind],as.factor(train.dat[,resp.ind]), method="rf")		# factor here
	pred	= predict(model,test.dat)
	AUC 	= AUC::auc( roc(pred, as.factor(test.dat[,resp.ind])) )
	
	return(AUC)
}

##### Wrapper for SVM algorithm
# kernel : type of kernel function to use
# etc.	 : refer to 'doRF'
doSVM = function(train.dat, test.dat, resp.ind, marker, kernel="radial") {
	require(e1071)
	resp 	= colnames(train.dat)[resp.ind]
	vars	= colnames(train.dat)
	f 		= as.formula(paste(resp, "~", paste(vars[which(marker==1)],collapse="+")))
	model 	= svm(f, data=train.dat, kernel=kernel, prob=T)
	pred	= predict(model,test.dat, decision.values=T)
	AUC 	= AUC::auc( roc(pred, as.factor(test.dat[,resp.ind])) )
	
	return(AUC)
}

##### wrapper for k Nearest Neighbors algorithm
# NN 	: number of neighbor items
# etc.	: refer to 'doRF'
doKNN = function(train.dat, test.dat, resp.ind, marker, NN) {
	# Using 'caret' package
	require(caret)
	resp 	= colnames(train.dat)[resp.ind]
	vars	= colnames(train.dat)
	f 		= as.formula(paste(resp, "~", paste(vars[which(marker==1)],collapse="+")))
	#model 	= train(f, data=train.dat,method="knn3",k=NN)
	#model 	= train(train.dat[,-5],train.dat[,5],method="knn",k=NN,prob=T)
	model 	= knn3(f,data=train.dat,k=NN,prob=T)
	pred	= predict(model,test.dat)
	pred 	= as.numeric(apply(pred,1,function(v) { colnames(pred)[which.max(v)] }))
	
	# Using 'knn' function directly
# 	require(class)
#  	pred 	= knn(train.dat[,-5],test.dat[,-5],train.dat[,5],k=NN,prob=T)

 	AUC 	= AUC::auc( roc(pred, as.factor(test.dat[,resp])) )
 	
 	return(AUC)
}

##### Wrapper for Generalized Logistic Regression (glm) function
# etc.	 : refer to 'doRF'
doLOGIT = function(train.dat, test.dat, resp.ind, marker) {
	#require(caret)
	require(ROCR)
	resp 	= colnames(train.dat)[resp.ind]
	vars	= colnames(train.dat)
	f 		= as.formula(paste(resp, "~", paste(vars[which(marker==1)],collapse="+")))
# 	model 	= train(f, data=train.dat,method="glm")
	model 	= glm(f,family=binomial(link=logit),data=train.dat)
	pred	= predict(model, test.dat, type="response")
# 	pred 	= prediction(pred, te.grp, label.ordering = NULL)			### TODO fix to remove te.grp
	AUC 	= AUC::auc( roc(pred, as.factor(test.dat[,resp.ind])) )
	
	return(AUC)
}


######################################################
### Wrapper function
######################################################

##### Generate a matrix w/ each row containing AUC values for RF, SVM, kNN, LOGIT
##### 	b/w the response variable and variables w/ indices specified in the
#####	corresponding column of 'marker.mat'.
##### CV is done in this function to generate 'train.dat' and 'test.dat'.
##### Given datasets & 'marker.mat' are preprocessed so that all rows w/ NA in data are removed.
# input			: input data.frame
# resp.ind		: Column index of the response variable in both 'train.dat' and 'test.dat'
# marker.mat	: A matrix w/ each row specifying which variables(cols) are included (1) or not (0)
# k				: # of folds for CV
# knn.NN		: # of nearest neighbors for kNN algorithm
# svm.kernel	: type of kernel function to be used for SVM
performance = function(input, resp.ind, marker.mat, k, knn.NN=3, svm.kernel='radial') {
	### Assert that the response variable column is a factor
	if(!is.factor(input[,resp.ind])) stop("Make sure the column for response variable is a factor.")
	
	### Definition
	result	= list()								# container for computed AUC results
	algs	= c("RF", "SVM", "kNN", 'LOGISTIC')		# ML algorithms for computing AUC
 	
	##### k-fold CV
	folds	= createFolds(input[,resp.ind],k=k,returnTrain=F)			# create folds w/ specified k

	for (i in 1:k) {
		train.dat 	= input[-folds[[i]],]
		test.dat 	= input[folds[[i]],]

	 	### Create 'res.all' to contain results
	 	cname = paste(algs, 'AUC', sep='_')
	 	res.all = matrix(nrow=nrow(marker.mat), ncol=length(algs))
	 	colnames(res.all) = cname
	 	rownames(res.all) = apply(marker.mat, 1, function(mm) {paste(which(mm==1), collapse = ',')})
	 	
	 	### Call each function to compute AUC's
	 	for( j in 1:nrow(marker.mat) ) {
	 		# Console progress indicator
	 		# cat('/ marker set', j, ' / ', nrow(marker.mat), ' total /')
	 		
	 		# Generate AUC values using the above submatrices and store in 'res.all'
	 		tmp.rf		<- doRF(train.dat, test.dat, resp.ind, marker.mat[j,])
	 		tmp.svm 	<- doSVM(train.dat, test.dat, resp.ind, marker.mat[j,], kernel=svm.kernel)
	 		tmp.knn 	<- doKNN(train.dat, test.dat, resp.ind, marker.mat[j,], NN=knn.NN)
	 		tmp.logit 	<- doLOGIT(train.dat, test.dat, resp.ind, marker.mat[j,])
	 		res.all[j,] <- c(tmp.rf, tmp.svm, tmp.knn, tmp.logit)
	 	}
	 	
	 	### Store AUC values to 'result'
 		result[[i]]	= res.all
	}

	### Compute mean AUC	
	result[["MEAN"]] = Reduce("+",result) / k
	
	return(result)
}


##### ============================================================================================
##### ============================================================================================
##### ============== Main Program ==================
##### ============================================================================================
##### ============================================================================================

##### Function call for stepwise AIC
result_AIC	= stepwiseAIC(input,resp.var)

##### Function call for regularized regression
result_auc	= regular.CV(input, resp.ind, crit="auc"); 	 result_auc
result_dev	= regular.CV(input, resp.ind, crit="dev"); 	 result_dev
result_cls	= regular.CV(input, resp.ind, crit="class"); result_cls
result_mae	= regular.CV(input, resp.ind, crit="mae"); 	 result_mae

##### 'marker.mat' generation
marker.mat 			 = rbind(result_AIC$step0,result_AIC$stepF,
					   		 result_auc$vars, result_dev$vars,
					   		 result_cls$vars, result_mae$vars)
colnames(marker.mat) = colnames(input[,-resp.ind])
rownames(marker.mat) = c("AIC_0","AIC_F","reg_auc","reg_dev","reg_cls","reg_mae")

##### Run comparison
eval.result = performance(input, resp.ind, marker.mat, k=5); eval.result


