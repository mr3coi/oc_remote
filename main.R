
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
library(car)

#require(devtools); devtools::install_github("briatte/ggnet")	### Install if necessary
library(ggnet); library(network); library(sna)

wd 		 <- "/Users/SJC/Documents/practice/internship/ovarian_cancer"
datafile <- "ova_db2.csv" 
varfile  <- "ova_variable.csv" 
setwd(wd)
set.seed(51)

source("preproc.R")
source("stepAIC.R")
source("lasso.R")
source("perform_eval.R")
source("util_fn.R")
source("corr_pair.R")

##### ============================================================================================
##### ============================================================================================
##### ============== Variables & Arguments Setting =============
##### ============================================================================================
##### ============================================================================================

##### Read and preprocess source data
ova.db.csv = as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 
ova.db.csv = ova.db.csv[1:235,]

##### Extract relevant portion of the given source data and preprocess it
col_index  = 1:103				### Column indices of explanatory variables to use
resp.var   = "Recurrence"		### Name of response variable
ova.db.csv = preProc(ova.db.csv,col_index,resp.var)

##### Extract relevant sub-data.frame from given data
	### Exclude columns that are meaningless, irrelevant, or redundant
avoid 	   = c("BMI","BMI_category1_2","BMI_category1_3","BMI_category1_4","Origin_2.1","Origin_3.1",
			   "Lymphocyte","Monocyte","Segmented_neutrophil","PLT")
avoid 	   = match(avoid,colnames(ova.db.csv),nomatch=ncol(ova.db.csv)+1)
	### Exclude all variables w/ column indices in 'avoid'
ova.db.csv = ova.db.csv[,-avoid]

##### Find indices of columns w/ (almost : 95%) uniform values and remove such columns
idx = which( apply(ova.db.csv, 2, function(v) { any(table(v)/(length(v)-sum(is.na(v))) > 0.95) }) )
#sapply(idx,function(i) table(ova.db.csv[,i]))			### Check distributions of variables
if (length(idx) > 0) ova.db.csv <- ova.db.csv[,-idx] 

##### Check if variables still have rare values (freq. < 0.01) and remove corresponding rows
# factor_check = unlist(lapply(ova.db.csv,is.factor))			## boolean for factor variables in 'ova.db.csv'
# tmp 	= apply(ova.db.csv,2,function(v) { any(table(v)/length(v) < 0.01) })
# tmp 	= names(tmp)[tmp]
# vars 	= tmp[tmp %in% names(factor_check)[factor_check]]
# 
# 	### Remove the rows containing such rare values
# if (length(vars) > 0) {
# 	idx = c()
# 	for (var in vars) {
# 		tmp 	= ova.db.csv[,var]
# 		tbl 	= table(tmp) / length(tmp); tbl
# 		target 	= names(tbl)[which(tbl < 0.01)]; target
# 		if (tbl[target] > 0) {
# 			#cat(colnames(ova.db.csv)[col], '\n')
# 			#cat(target, '---', which(tmp %in% target), '\n\n')
# 			idx = c(idx, match(target,tmp))
# 		}
# 	}
# 	idx = sort(unique(idx))
# 	ova.db.csv = ova.db.csv[-idx,] 
# }

##### Preprocess : remove NA's	=>	'input' (final data form)
row_NA	= apply(ova.db.csv,1,function(v) {sum(is.na(v)) == 0}); sum(row_NA)
input 	= ova.db.csv[row_NA,]
#input 	= input[,c(1:31,ncol(input))]		### For categories 1-7

##### Check that the response variable still survives, and compute its column index
resp.var %in% colnames(input)
resp.ind   = which(colnames(input) == resp.var)	# index of response variable

##### Partition data into training / test sets
inTrain 	= createDataPartition(input[,resp.ind],p=0.7,list=F)
test.data 	= input[-inTrain,]
train.data 	= input[inTrain,]

which(apply(train.data,2,function(v) any(table(v) == length(v))))
which(apply(test.data,2,function(v) any(table(v) == length(v))))
apply(train.data,2,table,useNA="always")
apply(test.data,2,table,useNA="always")

##### ============================================================================================
##### ============================================================================================
##### ============== Collinearity resolution ==================
##### ============================================================================================
##### ============================================================================================

exp.vars = colnames(train.data)[colnames(train.data)!=resp.var]

##### Correlation b/w surviving explanatory variables
#cor(train.data)

##### VIF calculation
vif_alias  = alias(glm(Recurrence ~ ., family=binomial(link=logit), data=train.data))	### Check for perfect MC

vif_result = vif(glm(Recurrence ~ ., family=binomial(link=logit), data=train.data))
vif_result = cbind(match(names(vif_result),exp.vars),vif_result)
colnames(vif_result) = c("exp.var_index","vif")
vif_result = vif_result[order(vif_result[,2],decreasing=T),]
vif_high = vif_result[vif_result[,2] > 5,1]
#names(vif_high) = NULL

##### Test for and remove collinear variables
	### Generate matrix of p-values from independence tests for each variable pair
corr.pairs  = pair.indep(train.data)

	### Plot # of pairs against thresholds for p-value to find optimum threshold
alpha_seq	= c(1 %o% 10^(-1:-15))
pair_count	= sapply(alpha_seq, function(v) { nrow(corr.pairs[corr.pairs[,3] < v,]) })
plot(-log10(alpha_seq),pair_count,xlab="-log10(alpha)",ylab="# of variable pairs",
	 main="dist. of dependent variable pairs\nper threshold")
abline(h=100)

	### Extract all pairs w/ p-value < threshold(alpha)	and sort by p-value in increasing order
alpha  = 0.05 / nrow(corr.pairs)	### alpha chosen by Bonferroni correction
abline(v = -log10(alpha))
corr.edit 	= corr.pairs[corr.pairs[,3] < alpha,]
corr.edit 	= corr.edit[order(corr.edit[,3]),]

	### Check names of and plot variable pairs to check whether they are truly correlated
# exp.vars = colnames(train.data)[colnames(train.data)!=resp.var]	## potential explanatory variables in 'train.data'
# var1.ind = 12
# var2.ind = 13
# exp.vars[c(var1.ind,var2.ind)]
# plot(train.data[,exp.vars[var1.ind]],train.data[,exp.vars[var2.ind]],xlab=exp.vars[var1.ind],ylab=exp.vars[var2.ind])

	### Plot a network graph for visualization of key variables
	### 	& Remove most highly linked nodes until each node is isolated
nodes	= sort(unique(as.numeric(corr.edit[,-3])))		# Column indices of vars in 'train.data'
net 	= matrix(0,nc=length(nodes),nr=length(nodes),dimnames=list(nodes,nodes))
for (i in 1:nrow(net)) {
	for (j in (i+1):ncol(net)) {
		if( any((corr.edit[,1] %in% nodes[i] & corr.edit[,2] %in% nodes[j]) |
			(corr.edit[,1] %in% nodes[j] & corr.edit[,2] %in% nodes[i])) ) net[i,j] = net[j,i] = 1
	}
}
exclude_node_val = c() 		### Container for indices of variables in 'train.data'
							###		that need to be removed due to dependency
include = (1:length(nodes))[!(1:length(nodes)) %in% match(exclude_node_val,nodes)]
net_obj = network(net[include,include],directed=F)
		### Highlight nodes w/ high VIF value
net_obj %v% "vif"   = ifelse(nodes[include] %in% vif_high, "MC", "non-MC")
net_obj %v% "color" = ifelse(net_obj %v% "vif" == "MC", "steelblue", "grey")
network.vertex.names(net_obj) = nodes[include]
ggnet2(net_obj,size=5,label=T,color="color")

while(T) {
	## include : INDICES of vars in 'nodes' that are currently surviving
	include = (1:length(nodes))[!(1:length(nodes)) %in% match(exclude_node_val,nodes)]
	
	## Draw the network graph b/w vars in 'include'
	net_obj = network(net[include,include],directed=F)
	# net_obj %v% "vif" = ifelse(include %in% vif_high, "MC", "non-MC")
	# network.vertex.names(net_obj) = nodes[include]
	# ggnet2(net_obj,size=5,label=T)
	
	## Check # of remaining links (exit if none)
	row.sum = apply(net[include,include],1,sum)
	if (all(row.sum == 0)) break
	
	## Exclude the 'hub' node (node w/ the most links)
	most_link_node = names(row.sum)[which.max(row.sum)]
	exclude_node_val = c(exclude_node_val,most_link_node)
}

	### Check final (isolated) state of variables
net_obj = network(net[include,include],directed=F)
		### Highlight nodes w/ high VIF value
net_obj %v% "vif"   = ifelse(nodes[include] %in% vif_high, "MC", "non-MC")
net_obj %v% "color" = ifelse(net_obj %v% "vif" == "MC", "steelblue", "grey")
network.vertex.names(net_obj) = nodes[include]
ggnet2(net_obj,size=5,label=T,color="color")
length(exclude_node_val); exclude_node_val

	### Update 'train.data' to exclude variables in 'exclude_node_val'
orig_train.data = train.data; dim(orig_train.data) 								### Dataset w/ variables not removed (for regularization)
train.data 	 	= train.data[,-as.numeric(exclude_node_val)]; dim(train.data)	### Dataset w/ variables removed (for AIC)
orig_resp.ind 	= which(colnames(orig_train.data) == resp.var)
resp.ind 	 	= which(colnames(train.data) == resp.var)



##### ============================================================================================
##### ============================================================================================
##### ============== Main Program ==================
##### ============================================================================================
##### ============================================================================================

##### Function call for stepwise AIC
#result_AIC	= stepwiseAIC(train.data,resp.var)
input = train.data
result_AIC = stepwiseAIC(input,resp.var)

##### Function call for regularized regression
eps = 10^(-3)		### Threshold size of coefficients
k 	= 5
result_auc	= regular.CV(orig_train.data, orig_resp.ind, k, thres=eps, crit="auc"); 	result_auc
result_dev	= regular.CV(orig_train.data, orig_resp.ind, k, thres=eps, crit="dev"); 	result_dev
result_cls	= regular.CV(orig_train.data, orig_resp.ind, k, thres=eps, crit="class");  result_cls
result_mae	= regular.CV(orig_train.data, orig_resp.ind, k, thres=eps, crit="mae"); 	result_mae

##### 'marker.mat' generation
AIC_0 = ifelse(exp.vars %in% colnames(train.data)[-ncol(train.data)][result_AIC$step0 == 1],1,0); AIC_0
AIC_F = ifelse(exp.vars %in% colnames(train.data)[-ncol(train.data)][result_AIC$stepF == 1],1,0); AIC_F

marker.mat 			 = rbind(AIC_0, AIC_F, result_auc$vars,
					   		 result_dev$vars, result_cls$vars, result_mae$vars)
colnames(marker.mat) = exp.vars
rownames(marker.mat) = c("AIC_0","AIC_F","reg_auc","reg_dev","reg_cls","reg_mae")
apply(marker.mat,1,sum)			### # of chosen variables
sort(apply(marker.mat,2,sum),decreasing=T)

##### Run comparison
	### CV.k : # of repetitions for CV
eval.result = performance(test.data, orig_resp.ind, marker.mat, CV.k=c(3,5))
eval.result[[3]]$MEAN_3; eval.result[[5]]$MEAN_5









##########################################################
########## Things to do #############
##########################################################

# - Check validity of independence-test-based variable reduction process (collinearity part)
#	- Add code to allow manually choosing variables based on prior knowledge
#	- Check why some p-values are >1, and whether it is ok to simply ignore them
#	- Discuss how to undo the 'tangle' b/w variables and choose the best ones
# - Data-related issues
#	- Check whether "Parity" and "Other_sites" are actually categorical variables
#	- Fix flawed data points in the given data
#	- Check if too many variables are removed during preprocessing
#	- Find a way to avoid removing rows again in 'Main Program' stage
# - Discuss the validity / significance of the calculated results
# - lasso.R
#	- Check whether separating training/test sets is necessary for variable selection
# - stepAIC.R
#	- Check the necessity of the codes commented out
# - perform_eval.R
#	- (doLOGIT) difference b/w 'predict' and 'prediction' objects?
# 	- Run kNN w/ various k values
