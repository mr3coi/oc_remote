
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
# test.dat	: dataset for test.dat the machine
# resp.ind	: index of the response variable in 'train.dat' and 'test.dat' datasets
# marker  	: vector indicating which columns to use (1) and not to use (0) as
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
##### Given datasets & 'marker.mat' arguments are assummed to be preprocessed so that all rows w/ NA in data are removed.
# input			: input data.frame for testing (evaluation)
# resp.ind		: Column index of the response variable in both 'train.dat' and 'test.dat'
# marker.mat	: A matrix w/ each row specifying which variables(cols) are included (1) or not (0)
# knn.k			: # of nearest neighbors for kNN algorithm
# svm.kernel	: type of kernel function to be used for SVM
performance = function(input, resp.ind, marker.mat, knn.k=5, svm.kernel='radial') {
	### Assert that the response variable column is a factor
	if(!is.factor(input[,resp.ind])) stop("Make sure the column for response variable is a factor.")
	
	### Definition
	algs	= c("RF", "SVM", "kNN", 'LOGISTIC')		# ML algorithms for computing AUC
 	
 	### Create 'res.all' to contain results
 	res.all 		  = matrix(nrow=nrow(marker.mat), ncol=length(algs))
 	colnames(res.all) = paste(algs, 'AUC', sep='_')
 	rownames(res.all) = apply(marker.mat, 1, function(mm) {paste(which(mm==1), collapse = ',')})
 	rownames(res.all) = paste0(rownames(marker.mat),'_',rownames(res.all))
 	
 	### Call each function to compute AUC's
 	for (j in 1:nrow(marker.mat)) {
 		### Console progress indicator
 		cat('/ marker set', j, ' / ', nrow(marker.mat), ' total /\n')
 		
 		### Generate AUC values using the above submatrices and store in 'res.all'
 		tmp.rf		<- doRF(train.dat, test.dat, resp.ind, marker.mat[j,])
 		tmp.svm 	<- doSVM(train.dat, test.dat, resp.ind, marker.mat[j,], kernel=svm.kernel)
 		tmp.knn 	<- doKNN(train.dat, test.dat, resp.ind, marker.mat[j,], NN=knn.k)
 		tmp.logit 	<- doLOGIT(train.dat, test.dat, resp.ind, marker.mat[j,])
 		res.all[j,] <- c(tmp.rf, tmp.svm, tmp.knn, tmp.logit)
 	}
	
	return(res.all)
}

