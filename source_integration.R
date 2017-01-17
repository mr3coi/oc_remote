
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

##### Locate source data
wd <- "/Users/SJC/Documents/practice/internship"
datafile <- "ova_db.csv" 
varfile <- "ova_variable.csv" 

setwd(wd) 
ova.db.csv = as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 

##### Stop from converting warnings to errors
options(warn=1)

##### ============================================================================================
##### ============================================================================================
##### =============== Utility Functions ===============
##### ============================================================================================
##### ============================================================================================

# Date-to-int conversion function 				##### TODO convert to # of days since 1900-01-01
date_conv = function(date_string) {
	date = as.vector(sapply(strsplit(date_string,split="\\."),strtoi))
	date = sum(date*c(10000,100,1))
	return(date)
}

search_colname = function(df, name) { return(grep(name, colnames(df))) }




##### ============================================================================================
##### ============================================================================================
##### ============ datafile : Preprocess ==============
##### ============================================================================================
##### ============================================================================================

# Change 'unknown' values to NA
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

# Specific deletions (flaws)
	# Delete unnamed columns
idx = which(colnames(ova.db.csv)=="")
	# Remove last two rows (empty rows)
ova.db.csv = ova.db.csv[ 1:236, -idx ]
	# Delete the repeated patient number column
ova.db.csv = ova.db.csv[,-which(colnames(ova.db.csv)=="Pt_No.1")]

# Alter column names to main titles (removing choices) and remove spaces within the main titles
colnames(ova.db.csv) = gsub(" ", "_", sapply(strsplit(colnames(ova.db.csv),"\\n"),function(v)trimws(v[1]))) 

# Find indices of columns w/ (almost) uniform values and remove such columns				### TODO almost
idx <- which( apply(ova.db.csv, 2, function(v) length(table(v)) == 1) )
ova.db.csv <- ova.db.csv[,-idx] 

# Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
options(warn=2);

# Find and remove/modify columns w/ non-numeric data
	#idx : non-numeric indices
idx <- which(apply(ova.db.csv, 2, function(v) class(try(as.numeric(v), T))=='try-error'))
options(warn=1);

	# Convert date info into numeric data (YYYYMMDD)
pattern = "\\d+\\.\\d+\\.\\d+"
idx2 = c()

for (col in idx) {
	if (sum(is.na(ova.db.csv[,col])) == 0 && length(grep(pattern,ova.db.csv[,col])) != 0) {		### TODO Optimize
		ova.db.csv[,col] = sapply(ova.db.csv[,col], function(str) {ifelse(str != "", date_conv(str), NA)} )
		idx2 = c(idx2,col)
	}
}
idx = idx[!idx %in% idx2]

	# Additional manipulation for 'Stage' column
ref = c(rep(1,3),rep(2,3),rep(3,3),4)
names(ref) = c("1a","1b","1c","2a","2b","2c","3a","3b","3c","4")
for (i in 1:nrow(ova.db.csv)) {
	ova.db.csv[i,"Stage"] = ref[ova.db.csv[i,"Stage"]]
}
idx = idx[names(idx)!="Stage"]

	# Additional manipulation for 'Recur_Site' column
ova.db.csv[,84] = sapply(ova.db.csv[,84],function(v) {paste(sort(unlist(strsplit(v,split=","))),collapse="")})
idx = idx[names(idx)!="Recur_Site"]

	# Extract numeric values at the beginning
pattern = "(\\d).*"
idx2 = c()
for (col in idx) {
	if (length(grep(pattern,ova.db.csv[,col])) != 0) {				### TODO Optimize
		ova.db.csv[,col] = sub(pattern,"\\1",ova.db.csv[,col])
		idx2 = c(idx2,col)
	}
}
idx = idx[!idx %in% idx2]

	# Delete the remaining non-numeric columns
ova.db.csv <- ova.db.csv[,-idx]

# Convert character values in ova.db.csv into numeric values
ova.db.csv <- as.data.frame(apply(ova.db.csv, 2, as.numeric))

# Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
options(warn=2)

# Check preprocessing results
apply(ova.db.csv,2,function(v) class(v))

# Convert 'Refractory progression during CTx' to 'Yes'
ova.db.csv$Recurrence[ova.db.csv$Recurrence==2] <- 1











##### ============================================================================================
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
set_num 	= 2
var_num 	= 1
resp_var 	= outcomes[outcome_idx[[set_num]]]
resp_var 	= resp_var[resp_var %in% colnames(ova.db.csv)][var_num]
resp.ind	= which(colnames(ova.db.csv) == resp_var)

##### Extract relevant sub-data.frame from given data	=>	'ova.db.csv.sub'
	### Avoid all columns corresponding to post-surgery
avoid 	= search_colname(ova.db.csv,"Start_Date_1st_regimen"):length(colnames(ova.db.csv))
	### Add other pre/during-surgery columns that are obviously meaningless or irrelevant
avoid 	= c(avoid, 1:3, 5, 43)
avoid 	= avoid[avoid != resp.ind]		# keep response variable if included in 'avoid'
ova.db.csv.sub = ova.db.csv[,-avoid]	# exclude all variables w/ column indices in 'avoid'

##### Preprocess : remove NA's	=>	'input'
row_NA	= apply(ova.db.csv.sub,1,function(v) {sum(is.na(v)) == 0})
input 	= ova.db.csv.sub[row_NA]






##### ============================================================================================
##### ============================================================================================
##### ============== Variable Selection : Stepwise AIC =============
##### ============================================================================================
##### ============================================================================================







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
	input.x	= input[,-resp.ind]
	input.y	= input[,resp.ind]
	
	inTrain	= createDataPartition(y,p=0.7,list=F)
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
	var.ind	= rep(0,ncol(x))
	var.ind[unlist(predict(model, newx=x.test, type="nonzero", s=l))] = 1
	coeffs	= as.vector(predict(model, newx=x.test, type="coefficients", s=l))[-1]
	
	return(list(lambda=l, vars=var.ind, coeffs=coeffs))
}

# ##### Function call
# result_auc	= regular.CV(x,y, resp.ind, crit="auc"); result_auc
# result_dev	= regular.CV(x,y, resp.ind, crit="dev"); result_dev
# result_cls	= regular.CV(x,y, resp.ind, crit="class"); result_cls
# result_mae	= regular.CV(x,y, resp.ind, crit="mae"); result_mae


##### ============================================================================================
##### ============================================================================================
##### ============== Main Program ==================
##### ============================================================================================
##### ============================================================================================



