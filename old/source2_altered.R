

doRF <- function(tr, tr.grp, te, te.grp){
	require(AUC)
	require(randomForest)
	
	print(tr.grp)
	print(te.grp)
	
	if(class(tr)=="numeric" & class(te)=="numeric"){
		model.rf <- randomForest(as.matrix(tr), as.factor(tr.grp))
		res.rf <- predict(model.rf, as.matrix(te), type="Prob")
		res.rf.table <- predict(model.rf, as.matrix(te))
	}else if(class(tr)=="matrix" & class(te)=="matrix"){
		model.rf <- randomForest(t(as.matrix(tr)), as.factor(tr.grp))
		res.rf <- predict(model.rf, t(as.matrix(te)), type="Prob")
		res.rf.table <- predict(model.rf, t(as.matrix(te)))
	}else break("Check data dimension")
	
	table.out <- table(res.rf.table, as.factor(te.grp))
	#   ACC <- sum(diag(table.out))/sum(table.out)
	#   AUC <- auc(roc(res.rf[,2], as.factor(te.grp)))
	#   SEN = table.out[1,1] / (table.out[1,1]+table.out[2,1])
	#   SPE = table.out[2,2] / (table.out[1,2]+table.out[2,2])
	#   return(cbind(ACC, SEN, SPE, AUC))
	AUC <- auc(roc(res.rf[,2], as.factor(te.grp)))
}
# Usage: doRF(exprdat1[1:10,], colnames(exprdat1), exprdat2sub[1:10,], colnames(exprdat2sub))

doSVM <- function(tr, tr.grp, te, te.grp, kernel="radial"){
	require(e1071)
	if(class(tr)=="numeric" & class(te)=="numeric"){
		model.svm <- svm(as.vector(tr), as.factor(tr.grp), kernel=kernel, probability=T)
		res.svm <- predict(model.svm, as.matrix(te), decision.values=T)
	}else if(class(tr)=="matrix" & class(te)=="matrix"){
		model.svm <- svm(matrix(t(tr), dim(tr)[2], dim(tr)[1]), as.factor(tr.grp), kernel=kernel, probability=T)
		res.svm <- predict(model.svm, t(as.matrix(te)), decision.values=T)
	}else break("Check data dimension")
	table.out <- table(res.svm, as.factor(te.grp))
	#   ACC <- sum(diag(table.out))/sum(table.out)
	#   AUC <- auc(roc(attributes(res.svm)$decision.values, as.factor(te.grp)))
	#   SEN = table.out[1,1] / (table.out[1,1]+table.out[2,1])
	#   SPE = table.out[2,2] / (table.out[1,2]+table.out[2,2])
	#   return(cbind(ACC, SEN, SPE, AUC))
	AUC <- auc(roc(attributes(res.svm)$decision.values, as.factor(te.grp)))
}
# Usage: doSVM(exprdat1[peptide.list.ix,], colnames(exprdat1), exprdat2[peptide.list.ix,], colnames(exprdat2))

doKNN <- function(tr, tr.grp, te, te.grp, NN){
	require(class)
	if(class(tr)=="numeric" & class(te)=="numeric"){
		res.knn <- knn(as.vector(tr), as.matrix(te), tr.grp, k=NN, prob=T)
	}else if(class(tr)=="matrix" & class(te)=="matrix"){
		res.knn <- knn(t(as.matrix(tr)), t(as.matrix(te)), tr.grp, k=NN, prob=T)
	}else break("Check data")
	table.out <- table(res.knn, te.grp)
	#   ACC <- sum(diag(table.out))/sum(table.out)
	#   AUC <- auc(roc(attributes(res.knn)$prob, as.factor(te.grp)))
	#   SEN = table.out[1,1] / (table.out[1,1]+table.out[2,1])
	#   SPE = table.out[2,2] / (table.out[1,2]+table.out[2,2])
	#   return(cbind(ACC, SEN, SPE, AUC))
	AUC <- auc(roc(attributes(res.knn)$prob, as.factor(te.grp)))
}

doLOGIT = function(tr, tr.grp, te, te.grp, tr.cov = NULL, te.cov = NULL, mark) {
	require(ROCR)
	if( !is.null(tr.cov) | !is.null(te.cov) ) {
		if(class(tr)=="numeric" & class(te)=="numeric"){
			dat.tr = data.frame(cond = as.numeric(tr.grp), (tr), t(tr.cov))
			dat.te = data.frame(cond = as.numeric(te.grp), (te), t(te.cov))
		}else if(class(tr)=="matrix" & class(te)=="matrix"){
			dat.tr = data.frame(cond = as.numeric(tr.grp), t(tr), t(tr.cov))
			dat.te = data.frame(cond = as.numeric(te.grp), t(te), t(te.cov))
		}else stop("Check data dimension")
		colnames(dat.tr) = colnames(dat.te) = c('cond', mark, 'sex', 'age')
	} else {
		if(class(tr)=="numeric" & class(te)=="numeric"){
			dat.tr = data.frame(cond = as.numeric(tr.grp), (tr))
			dat.te = data.frame(cond = as.numeric(te.grp), (te))
		}else if(class(tr)=="matrix" & class(te)=="matrix"){
			dat.tr = data.frame(cond = as.numeric(tr.grp), t(tr))
			dat.te = data.frame(cond = as.numeric(te.grp), t(te))
		}else stop("Check data dimension")
		colnames(dat.tr) = colnames(dat.te) = c('cond', mark)
	}
	model = glm(cond ~ ., family=binomial(link = logit) , data = dat.tr)
	pred=prediction(predict(model,dat.te,type="response"), te.grp, label.ordering = NULL)
	AUC=performance(pred,measure="auc")@y.values[[1]]
}




###########################
### fitting function
###########################
# train.dat, test.dat : p X n
# 각각의 column names에는 sample id가 아니라 sample label (0 or 1) 이 붙어있는 구조에요.
# res.all			  : comb X p

# machine.learning 쪽 function은 comb * p dimension의 matrix이고요,
# i번째 (i = 1, ..., p) column에는 observation matrix (train.dat 및 test.dat)의
# row index가 들어있어요.

machine.learning = function(train.dat, test.dat, marker.mat, knn.NN=3, svm.kernel='radial') {
	tr.class = colnames(train.dat)
	te.class = colnames(test.dat)
	
	# Create 'res.all' to contain results
	cname = paste(c("RF", "SVM", "kNN", 'LOGISTIC'), 'AUC', sep='_')
	res.all = matrix(NA, nrow=nrow(marker.mat), ncol=length(cname))
	colnames(res.all) = cname
	rownames(res.all) = apply(marker.mat, 1, function(mm) {paste(which(mm==1), collapse = ', ')})
	
	for( i in 1:nrow(marker.mat) ) {
		cat('/ marker set', i, ' / ', nrow(marker.mat), ' total /')

		tr.dat = train.dat[marker.mat[i,],]
		te.dat = test.dat[marker.mat[i,],]

		if( sum(is.na(tr.dat)) == 0 & sum(is.na(te.dat)) == 0 ) {
			tmp.rf <- doRF(tr.dat, tr.class, te.dat, te.class)
			tmp.svm <- doSVM(tr.dat, tr.class, te.dat, te.class, kernel=svm.kernel)
			tmp.knn <- doKNN(tr.dat, tr.class, te.dat, te.class, NN=knn.NN)
			tmp.logit <- doLOGIT(tr.dat, tr.class, te.dat, te.class, mark=marker.mat[,i])
			res.all[i,] <- c(tmp.rf, tmp.svm, tmp.knn, tmp.logit)
		} else next
	}
	cat('\n')
	
	# Return
	res.all
}


###########################
### Main Program (including CV)
###########################
data(iris); library(caret)

preObj = preProcess(iris[,-5],method=c("center","scale"))
input = predict(preObj,iris)
input[,5] = input[,5]=="setosa"
input[,-5] = input[,-5]>=0
input = apply(input,2,function(v) ifelse(v,1,0))

k = 3
folds = createFolds(input[,5],k=k)
marker = diag(x=1,nc=4,nr=4)			# deal w/ marker
result = matrix(nr=ncol(marker),nc=4)

for (i in 1:k) {
	train  = t(input[-folds[[i]],])
	test   = t(input[folds[[i]],])
	result = result + machine.learning(train,test,marker.mat)		# add CV specifics (mean?)
}
result = result / 4