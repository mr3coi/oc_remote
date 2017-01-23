
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

#require(devtools); devtools::install_github("briatte/ggnet")	### Install if necessary
library(ggnet); library(network); library(sna)

wd 		 <- "/Users/SJC/Documents/practice/internship/ovarian_cancer"
datafile <- "ova_db.csv" 
varfile  <- "ova_variable.csv" 
setwd(wd) 

source("preproc.R")
source("stepAIC.R")
source("lasso.R")
source("perform_eval.R")
source("util_fn.R")

##### ============================================================================================
##### ============================================================================================
##### ============== Variables & Arguments Setting =============
##### ============================================================================================
##### ============================================================================================

##### Read and preprocess source data
ova.db.csv = as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 
ova.db.csv = preProc(ova.db.csv)

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
	### Avoid other pre/during-surgery columns that are obviously meaningless or irrelevant
avoid 		= c(avoid, 1:3, 5, 43)
avoid 		= avoid[avoid != resp.ind]		# keep response variable if included in 'avoid'
ova.db.csv.sub = ova.db.csv[,-avoid]		# exclude all variables w/ column indices in 'avoid'
resp.ind	= which(colnames(ova.db.csv.sub) == resp.var)	# update 'resp.ind'

##### Preprocess : remove NA's	=>	'input' (final data form)
row_NA		= apply(ova.db.csv.sub,1,function(v) {sum(is.na(v)) == 0})
input 		= ova.db.csv.sub[row_NA,]

##### temporary storage of 'input' for easier reload
#save(input,file='input.RData')
#load(file='input.RData')




##### ============================================================================================
##### ============================================================================================
##### ============== Collinearity resolution ==================
##### ============================================================================================
##### ============================================================================================

##### Test for and remove collinear variables
corr.pairs = matrix(nc=3,nr=0)
colnames(corr.pairs) = c("var1","var2","p.value")
exp.vars = colnames(input)[colnames(input)!=resp.var]
options(warn=0)
require(car)

for (var1.ind in 1:(length(exp.vars)-1)) {
	var1 = exp.vars[var1.ind]
	
	for (var2.ind in (var1.ind+1):length(exp.vars)) {
		var2 = exp.vars[var2.ind]

		### b/w categorical variables	=> Chi-sq. independence test
		if (prod(factor_check[c(var1,var2)]) == 1) {
			tbl = table(input[,var1],input[,var2])
			corr.pairs = rbind(corr.pairs, c(var1.ind,var2.ind, chisq.test(tbl)$p.value))
		}
		### b/w numerical variables		=> Linear regression independence test
		else if (sum(factor_check[c(var1,var2)]) == 0) {
			#cat(var1,var2,'\n')
			f 		= as.formula(paste(var1,var2,sep="~"))
			model 	= lm(f,data=input)
			p.value	= as.data.frame(anova(model))[1,5]
			corr.pairs = rbind(corr.pairs, c(var1.ind,var2.ind, p.value))
		}
		### b/w a numerical variable and a categorical variable
		###		=> Logistic regression independence test
		else {
			categ = ifelse(factor_check[var1]==1, var1, var2)
			numer = ifelse(factor_check[var1]==1, var2, var1)
			lv	  = length(unique(input[,categ]))
			f 		= as.formula(paste(categ, numer, sep="~"))
			cat("cat : '", categ, "', num : '", numer, "', lv : ", lv, '\n')
			
			if (lv == 1) { p.value = NaN }	# Uniform column => pass
			else if (lv == 2) { 			# Binomial Logistic Regression (assuming normality)
				model 	= glm(f,data=input,family=binomial(link=logit))
				p.value	= as.data.frame(anova(model))[1,5]
			}
			else {							# Multinomial Logistic Regression
				input$relv	= relevel(input[,categ], ref=levels(input[,categ])[1])
				f 			= as.formula(paste("relv", numer, sep="~"))
				model 	  	= multinom(f, data=input)
				## Use Wald Z-test under assumption of normality
				# wald_z		= summary(model)$coefficients / summary(model)$standard.errors
				# p.value 	= (1 - pnorm(abs(wald_z),0,1)) * 2	# Both ends
				p.value 	= unlist(Anova(model))[3]
			}
			corr.pairs = rbind(corr.pairs, c(var1.ind,var2.ind, p.value))
		}
	}
}
	### tidy-up (restoration to original state)
input = input[,-length(colnames(input))]
options(warn=1)

	### Remove defected variables (var.s causing NaN p-values)
if (any(is.nan(corr.pairs[,"p.value"]))) {
	NaN_tbl = table(corr.pairs[which(is.nan(corr.pairs[,"p.value"])),-3])
	ranking = sort(rank(NaN_tbl,ties.method="min"),decreasing=T)
	idx		= as.numeric(names(ranking)[which(ranking == ranking[1])])
	corr.pairs = corr.pairs[!(corr.pairs[,1] %in% idx | corr.pairs[,2] %in% idx),]
}

# if (any(corr.pairs[,"p.value"] > 1)) {
# 	over_tbl = table(corr.pairs[which( corr.pairs[,"p.value"] > 1 ),-3])
# 	ranking = sort(rank(over_tbl,ties.method="min"),decreasing=T)
# 	idx		= as.numeric(names(ranking)[which(ranking == ranking[1])])
# 	corr.pairs = corr.pairs[!(corr.pairs[,1] %in% idx | corr.pairs[,2] %in% idx),]
# }

	### temporary storage of 'corr.pairs'(full) for easier reload
#save(corr.pairs,file='corr.pairs.RData')
#load(file='corr.pairs.RData')

	### Plot # of pairs against thresholds for p-value to find optimum threshold
alpha_seq	= c(1 %o% 10^(-3:-15))
pair_count	= sapply(alpha_seq, function(v) { nrow(corr.pairs[corr.pairs[,3] < v,]) })
plot(-log10(alpha_seq),pair_count)

	### Extract all pairs w/ p-value < threshold(alpha)	and sort by p-value in increasing order
alpha 		= 1e-8
corr.edit 	= corr.pairs[corr.pairs[,3] < alpha,]
corr.edit 	= corr.edit[order(corr.edit[,3]),]

	### Check names of and plot variable pairs to check whether they are truly correlated
# var1.ind = 21
# var2.ind = 22
# exp.vars[c(var1.ind,var2.ind)]
# plot(input[,exp.vars[var1.ind]],input[,exp.vars[var2.ind]])

	### Plot a network graph for visualization of key variables
nodes	= sort(unique(as.numeric(corr.edit[,-3])))
#plot(table(as.numeric(corr.edit[,-3])))
net 	= matrix(0,nc=length(nodes),nr=length(nodes),dimnames=list(nodes,nodes))
for (i in 1:nrow(net)) {
	for (j in (i+1):ncol(net)) {
		if( any((corr.edit[,1] %in% nodes[i] & corr.edit[,2] %in% nodes[j]) |
			(corr.edit[,1] %in% nodes[j] & corr.edit[,2] %in% nodes[i])) ) net[i,j] = net[j,i] = 1
	}
}

	### Remove most highly linked nodes until each node is isolated
exclude_node_val = c()

while(T) {
	include = (1:length(nodes))[!(1:length(nodes)) %in% match(exclude_node_val,nodes)]
	net_obj = network(net[include,include],directed=F)
	network.vertex.names(net_obj) = nodes[include]
	ggnet2(net_obj,size=5,label=T)
	row.sum = apply(net[include,include],1,sum)
	if (all(row.sum == 0)) break
	most_link_node = names(row.sum)[which.max(row.sum)]
	exclude_node_val = c(exclude_node_val,most_link_node)
}

	### Check final (isolated) state of variables
ggnet2(net_obj,size=5,label=T)
length(include); include

	### Update 'input' to only contain the surviving explanatory variables (in 'include')

input 	 = input[,c(include,resp.ind)]; head(input)
resp.ind = which(colnames(input) == resp.var)



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

