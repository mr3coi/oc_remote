doRF <- function(tr, tr.grp, te, te.grp){
	require(AUC)
	require(randomForest)
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
machine.learning = function(train.dat, test.dat, marker.mat, knn.NN=3, svm.kernel='radial') {
	tr.class = colnames(train.dat)
	te.class = colnames(test.dat)
	cname = paste(c("RF", "SVM", "kNN", 'LOGISTIC'), 'AUC', sep='_')
	res.all = matrix(NA, nrow=dim(marker.mat)[2], ncol=length(cname))
	colnames(res.all) = cname
	rownames(res.all) = apply(marker.mat, 2, function(mm) {paste(mm, collapse = ', ')})
	for( i in 1:dim(marker.mat)[2] ) {
		cat('/ marker ', i, ' / ', dim(marker.mat)[2], ' total /')
		tr.dat = train.dat[marker.mat[,i],]
		te.dat = test.dat[marker.mat[,i],]
		na.1 = is.na(tr.dat)
		na.2 = is.na(te.dat)
		if( sum(na.1) + sum(na.2) == 0 ) {
			tmp.rf <- doRF(tr.dat, tr.class, te.dat, te.class)
			tmp.svm <- doSVM(tr.dat, tr.class, te.dat, te.class, kernel=svm.kernel)
			tmp.knn <- doKNN(tr.dat, tr.class, te.dat, te.class, NN=knn.NN)
			tmp.logit <- doLOGIT(tr.dat, tr.class, te.dat, te.class, mark=marker.mat[,i])
			res.all[i,] <- c(tmp.rf, tmp.svm, tmp.knn, tmp.logit)
		} else next
	}
	cat('\n')
	res.all
}