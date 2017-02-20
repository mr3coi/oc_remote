
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

par(mfrow=c(1,1)); graphics.off(); rm(list=ls())
wd 		 <- "/Users/SJC/Documents/practice/internship/ovarian_cancer"
datafile <- "ova_db2.csv" 
varfile  <- "ova_variable.csv" 
setwd(wd)
seed = 51
set.seed(seed)

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
#resp.var   = "Recurrence"		### Name of response variable
resp.var   = "Platinum_user_only_6mo"
ova.db.csv = preProc(ova.db.csv,col_index,resp.var)

##### Extract relevant sub-data.frame from given data
	### Exclude columns that are meaningless, irrelevant, or redundant
avoid 	   = c("BMI","BMI_category1_2","BMI_category1_3","BMI_category1_4","Origin_2.1","Origin_3.1","Origin.1_3",
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
#input 	= input[,c(1:28,ncol(input))]		### For categories 1-7

##### Check that the response variable still survives, and compute its column index
resp.var %in% colnames(input)
resp.ind   = which(colnames(input) == resp.var)	# index of response variable

##### ============================================================================================
##### ============================================================================================
##### ============== Collinearity resolution ==================
##### ============================================================================================
##### ============================================================================================

exp.vars = colnames(input)[colnames(input)!=resp.var]

##### Correlation b/w surviving explanatory variables
#cor(input)

##### VIF calculation
f = paste(resp.var,"~",".",sep=" ")
vif_alias  = alias(glm(f, family=binomial(link=logit), data=input))	### Check for perfect MC

vif_result = vif(glm(f, family=binomial(link=logit), data=input))
vif_result = cbind(match(names(vif_result),exp.vars),vif_result)
colnames(vif_result) = c("exp.var_index","vif")
vif_result = vif_result[order(vif_result[,2],decreasing=T),]
vif_high = vif_result[vif_result[,2] > 5,1]
#names(vif_high) = NULL

##### Test for and remove collinear variables
	### Generate matrix of p-values from independence tests for each variable pair
corr.pairs  = pair.indep(input)

	### Plot # of pairs against thresholds for p-value to find optimum threshold
alpha_seq	= c(1 %o% 10^(-3:-15))
pair_count	= sapply(alpha_seq, function(v) { nrow(corr.pairs[corr.pairs[,3] < v,]) })
plot(-log10(alpha_seq),pair_count,xlab="-log10(alpha)",ylab="# of variable pairs",
	 main="dist. of dependent variable pairs\nper threshold")
#abline(h=100)

	### Extract all pairs w/ p-value < threshold(alpha)	and sort by p-value in increasing order
alpha  = 0.05 / nrow(corr.pairs)	### alpha chosen by Bonferroni correction
abline(v = -log10(alpha))
corr.edit 	= corr.pairs[corr.pairs[,3] < alpha,]
corr.edit 	= corr.edit[order(corr.edit[,3]),]

	### Check names of and plot variable pairs to check whether they are truly correlated
# exp.vars = colnames(input)[colnames(input)!=resp.var]	## potential explanatory variables in 'input'
# var1.ind = 12
# var2.ind = 13
# exp.vars[c(var1.ind,var2.ind)]
# plot(input[,exp.vars[var1.ind]],input[,exp.vars[var2.ind]],xlab=exp.vars[var1.ind],ylab=exp.vars[var2.ind])

	### Plot a network graph for visualization of key variables
	### 	& Remove most highly linked nodes until each node is isolated
nodes	= sort(unique(as.numeric(corr.edit[,-3])))		# Column indices of vars in 'input'
net 	= matrix(0,nc=length(nodes),nr=length(nodes),dimnames=list(nodes,nodes))
for (i in 1:nrow(net)) {
	for (j in (i+1):ncol(net)) {
		if( any((corr.edit[,1] %in% nodes[i] & corr.edit[,2] %in% nodes[j]) |
			(corr.edit[,1] %in% nodes[j] & corr.edit[,2] %in% nodes[i])) ) net[i,j] = net[j,i] = 1
	}
}
exclude_node_val = c() 		### Container for indices of variables in 'input'
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
length(exclude_node_val); colnames(input)[as.numeric(exclude_node_val)]

	### Update 'input' to exclude variables in 'exclude_node_val'
#orig_input 		= input; dim(orig_input) 							### Dataset w/ variables not removed (for regularization)
reduc_input 	= input[,-as.numeric(exclude_node_val)]; dim(reduc_input)	### Dataset w/ variables removed (for AIC)
orig_resp.ind 	= which(colnames(input) == resp.var)
reduc_resp.ind 	= which(colnames(reduc_input) == resp.var)



##### ============================================================================================
##### ============================================================================================
##### ============== Main Program ==================
##### ============================================================================================
##### ============================================================================================

### Function call for stepwise AIC
set.seed(seed)
result_AIC	= stepwiseAIC(reduc_input,resp.var)

### Function call for regularized regression
eps = 10^(-3)		### Threshold size of coefficients
k 	= 3
par(mfrow=c(2,2))
set.seed(seed)
result_auc	= regular.CV(input, orig_resp.ind, k, thres=eps, crit="auc"); 	result_auc
set.seed(seed)
result_dev	= regular.CV(input, orig_resp.ind, k, thres=eps, crit="dev"); 	result_dev
set.seed(seed)
result_cls	= regular.CV(input, orig_resp.ind, k, thres=eps, crit="class"); result_cls
set.seed(seed)
result_mae	= regular.CV(input, orig_resp.ind, k, thres=eps, crit="mae"); 	result_mae

### 'marker.mat' generation
AIC_0 = ifelse(exp.vars %in% colnames(reduc_input)[-ncol(reduc_input)][result_AIC$step0 == 1],1,0); AIC_0
AIC_F = ifelse(exp.vars %in% colnames(reduc_input)[-ncol(reduc_input)][result_AIC$stepF == 1],1,0); AIC_F

marker.mat 			 = rbind(AIC_0, AIC_F, result_auc$vars, result_dev$vars, result_cls$vars, result_mae$vars)
colnames(marker.mat) = exp.vars
rownames(marker.mat) = c("AIC_0","AIC_F","reg_auc","reg_dev","reg_cls","reg_mae")

### Run comparison
	### CV.k : # of repetitions for CV
eval.result = performance(input, orig_resp.ind, marker.mat, CV.k=3)
#eval.result[[5]]$MEAN_5

### Output results
apply(marker.mat,1,sum)							### # of chosen variables per each variable selection method
sort(apply(marker.mat,2,sum),decreasing=T)		### Variables in order of frequent inclusion
markers = apply(marker.mat,1, function(row) colnames(marker.mat)[row==1]); markers
eval.result[[3]]$MEAN_3
