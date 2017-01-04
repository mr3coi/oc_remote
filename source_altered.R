
library(MASS)
library(ROCR)

wd <- "/Users/SJC/Documents/practice/internship"
datafile <- "ova_db.csv" 
varfile <- "ova_variable.csv" 

setwd(wd) 

##### ============================================================================================
##### ============================================================================================
##### =============== Utility Functions ===============
##### ============================================================================================
##### ============================================================================================

# Date-to-int conversion function (days since 1900-01-01)
date_conv = function(date_string) {
	date = as.vector(sapply(strsplit(date_string,split="\\."),strtoi))
	date = sum(date*c(10000,100,1))
	return(date)
}

##### ============================================================================================
##### ============================================================================================
##### ============ datafile : Preprocess ==============
##### ============================================================================================
##### ============================================================================================

# Read data file
ova.db.csv = as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 

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

# # Change empty strings to NA
# ova.db.csv = apply(ova.db.csv,2,
# 				   function(v) { if(class(try(as.numeric(v), T))=='try-error') 
# 				   					sapply(v, function(vv) { if(vv == "") vv = NA }) })


# Specific deletions (flaws)
	# Delete unnamed columns
idx = which(colnames(ova.db.csv)=="")
	# Remove last two rows (empty rows)
ova.db.csv = ova.db.csv[ 1:236, -idx ]
	# Delete the repeated patient number column
ova.db.csv = ova.db.csv[,-which(colnames(ova.db.csv)=="Pt_No.1")]

# Alter column names to main titles (removing choices) and remove spaces within the main titles
colnames(ova.db.csv) = gsub(" ", "_", sapply(strsplit(colnames(ova.db.csv),"\\n"),function(v)trimws(v[1]))) 

# Find indices of columns w/ (almost) uniform values and remove such columns								### ALMOST
idx <- which( apply(ova.db.csv, 2, function(v) length(table(v)) == 1) )
ova.db.csv <- ova.db.csv[,-idx] 

# Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
options(warn=2);

# Find and remove/modify columns w/ non-numeric data
idx <- which(apply(ova.db.csv, 2, function(v) class(try(as.numeric(v), T))=='try-error'))	# non-numeric indices
options(warn=1);

	# Convert date info into numeric data (YYYYMMDD)
pattern = "\\d+\\.\\d+\\.\\d+"
idx2 = c()

for (col in idx) {
	if (sum(is.na(ova.db.csv[,col])) == 0 && length(grep(pattern,ova.db.csv[,col])) != 0) {
		ova.db.csv[,col] = sapply(ova.db.csv[,col], function(str) {ifelse(str != "", date_conv(str), NA)} )
		idx2 = c(idx2,col)
	}
}
idx = idx[!idx %in% idx2]

#idx[!idx %in% idx2]
#length(colnames(ova.db.csv))

	# Additional manipulation for 'Stage' column
ref = 1:10
names(ref) = c("1a","1b","1c","2a","2b","2c","3a","3b","3c","4")
for (i in 1:nrow(ova.db.csv)) {
	ova.db.csv[i,"Stage"] = ref[ova.db.csv[i,"Stage"]]
}
idx = idx[names(idx)!="Stage"]

	# Extract numeric values at the beginning
pattern = c()			############################ EDIT ################################
	
	# Delete the remaining non-numeric columns
ova.db.csv <- ova.db.csv[,-idx]

#length(colnames(ova.db.csv))

# Find and remove columns w/ NA/total ratio > 0.2
idx <- which(apply(ova.db.csv, 2, function(v) sum(is.na(v)))/nrow(ova.db.csv) > 0.2)
ova.db.csv <- ova.db.csv[,-idx]
options(warn=2)

#length(colnames(ova.db.csv))

# Convert character values in ova.db.csv into numeric values
ova.db.csv <- as.data.frame(apply(ova.db.csv, 2, as.numeric))

# Convert 'Refractory progression during CTx' to 'Yes'
ova.db.csv$Recurrence[ova.db.csv$Recurrence==2] <- 1

# n : Size of data in ova.db.csv.sub
n <- nrow(ova.db.csv.sub)

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
##### ============ Data manipulation ==============
##### ============================================================================================
##### ============================================================================================

# Find variables in ova.db.csv w/ colnames corresponding to those listed in each element of 'vars'
cols <- list()
for (i in 1:ncol(vars)) {
  cols[[i]] <- na.omit(colnames(ova.db.csv)[colnames(ova.db.csv) %in% vars[,i]])
  #ova.db.csv[,cols[[i]]]
}
vars[which(is.na(match(tolower(vars[,1]), tolower(colnames(ova.db.csv))))), 1]    # ???????

# Choose response variable and explanatory variable sets
resp_var = 1

# # Generate a count table of values per each variable listed in element 'resp_var' of 'col'
# tbl <- apply(ova.db.csv[,cols[[resp_var]]], 2, table)
# 
# # newcols : Exclude from tbl the names of binary variables w/ ratio b/w major and minor items < 0.05
# newcols <- names(which(unlist(lapply(tbl, function(v) {
#   length(v) != 2 | min(v[1]/v[2], v[2]/v[1]) >= 0.05
# }))))

# # newcols2 : Choose variables among newcols w/o NA
# newcols2 <- names(which((colSums(is.na(ova.db.csv[,newcols])) / nrow(ova.db.csv)) == 0)) 

##### ============================================================================================
##### ============================================================================================
##### ============ Data Analysis : Correlation check ==============
##### ============================================================================================
##### ============================================================================================

# Calculate inter-variable correlation and plot significant correlations - Unnecessary
# if (0) {
#   # Calculate pairwise correlation b/w variables in 'cols[[resp_var]]'
#   vv <- cor(ova.db.csv[,cols[[resp_var]]], use="pair")
#   
#   # Nullify correlations < 0.5
#   vv[abs(vv) < 0.5] <- 0
#   
#   # Find variables w/ correlation sum > 1 and discard other variables
#   vx <- vv[colSums(abs(vv)) > 1, colSums(abs(vv)) > 1] 
#   
#   # Plot the above data and draw grids
#   image(vx, xaxt='n', yaxt='n', zlim=c(-1,1)); box()
#   axis(1, seq(0, 1, length.out=nrow(vx)), labels=rep("", nrow(vx)))
#   text(seq(0, 1, length.out=nrow(vx)), par("usr")[1]-0.05, srt=-45, adj=0, labels=colnames(vx), xpd=T)
# }

#colSums(is.na(ova.db.csv[,newcols2])) / nrow(ova.db.csv)       # ???????


##### ============================================================================================
##### ============================================================================================
##### ============ Data Analysis : Regression ==============
##### ============================================================================================
##### ============================================================================================

# Interested response variables
outcomes = c("Recurrence",
			  "Platinum_resistance_6mo", "Platinum-resistance_group",
			  "End_Date_1st_regimen", "Op_date",
			  "Dx_date", "Death", "Last_FU_Date", "Expired_Date",
			  "Residual_tumor_site_1st_debulking", "Residual_tumor_size_1st_debulking",
			  "PLN_status", "PALN_status")
outcome_index = c("CR",
				  "CS", "CT",
				  "CO", "BF",
				  "E", "DC", "DD", "DE",
				  "CI", "CJ",
				  "AX", "AY")
outcome_idx = 1

# Submatrix of ova.db.csv of rows w/ "Recurrence" value != NA
ova.db.csv.sub <- ova.db.csv[!is.na(ova.db.csv[outcomes[outcome_idx]]),]

# Run regressions b/w the interested response variable and the explanatory variables in 'vars[,1]'
mod <- glm(paste(outcomes[outcome_idx], "~", paste(cols[[resp_var]],collapse="+")),
		   data = ova.db.csv.sub, family = "binomial")
nul <- glm(paste(outcomes[outcome_idx], "~1"), data = ova.db.csv.sub)
step.0 <- stepAIC(nul, scope=list(lower=nul,upper=mod), direction="both")
#step.0x <- step(nul, scope=list(lower=nul,upper=mod), direction="both")
step.f <- stepAIC(mod, direction="both")
#step.fx <- step(mod, direction="both")

# Show ANOVA tables from the above regressions
step.0$anova
step.f$anova


# Picked 
var.full <- setdiff(cols[[resp_var]], substr(as.character(step.f$anova$Step), 3, 1000))
mean(do.cv(var.full, 3)$auc)
var.full <- substr(as.character(step.0$anova$Step), 3, 1000)[-1]
mean(do.cv(var.full, 3)$auc)


var.full <- setdiff(cols[[resp_var]], substr(as.character(step.fx$anova$Step), 3, 1000))
mean(do.cv(var.full, 3)$auc)
var.full <- substr(as.character(step.0x$anova$Step), 3, 1000)[-1]
mean(do.cv(var.full, 3)$auc)





##### ============================================================================================
##### ============================================================================================
# =============== Evaluate (Cross-Validation) ===============
##### ============================================================================================
##### ============================================================================================

# Function for cross validation
do.cv <- function(var, ncv) {
  ret <- c()
  for (i in 1:100) {
    # Sample round(n/ncv) rows to exclude from ova.db.csv.sub while model generation
    idx <- sample(1:n)[1:round(n/ncv)]
    curmod <- glm(paste(outcomes[outcome_idx], "~",paste(var,collapse="+")),
    			  data=ova.db.csv.sub[-idx,], family="binomial")
    
    # Generate predictions using the above model and previously-excluded rows in ova.db.csv.sub
    curpr <- predict(curmod, ova.db.csv.sub[idx,], type="response")
    curpr <- prediction(curpr, ova.db.csv.sub[idx,outcomes[outcome_idx]])
    
    # Calculate the performance of the above model
    prof <- performance(curpr, measure="tpr", x.measure="fpr")
    auc <- performance(curpr, measure = "auc")@y.values[[1]]
    ret <- c(ret, auc)
  }
  list(
    v = prof,
    auc = ret
  )
}

# Evaluate the performance of individual variables in cols[[resp_var]]
rex <- matrix(nr=length(cols[[resp_var]]),nc=100)
for (i in 1:length(cols[[resp_var]])) {
  ret <- do.cv(cols[[resp_var]][i], 3)$auc
  rex[i,] = ret
}
rownames(rex) <- cols[[resp_var]]

# Discard variables w/ mean AUC <= 0.7
req <- t(rex)[,which(colMeans(t(rex)) > 0.7)]

# Rearrange rows of 'req' in decreasing order (column means)
req <- req[,order(colMeans(req), decreasing=T)]





##### ============================================================================================
##### ============================================================================================
##### =============== Visualize ===============
##### ============================================================================================
##### ============================================================================================

# Define graphic parameters 'xx' and 'xxx'
xxx <- xx <- par("mai")
xxx[1] <- par("mai")[1] + 0.5
xxx[2] <- par("mai")[2] + 1

# Set graphic parameters to 'xxx'
par(mai=xxx)

# Plot the colMeans of 'req' and draw dotted gridlines into the plot
barplot(colMeans(req), ylim=c(0,1), xaxt='n', xpd=F, space=1, ylab="AUC", 
        main="Individual variable prediction performance w/ AUC < 0.6"); box()
grid(NA, 5, lwd=2)

# Plot the colMeans in another way(???)
par(new=T)
barplot(colMeans(req), ylim=c(0,1), xaxt='n', xpd=F, space=1, ylab="AUC", 
        main="Individual variable prediction performance w/ AUC<0.6");box();

# Add variable labels to plot
end_point = 0.5 + ncol(req) + ncol(req) - 1 
text(seq(1.5,end_point,by=2), par("usr")[3]-0.03,  
     srt = 25, adj= 1, xpd = TRUE, 
     labels = paste(colnames(req)), cex=.8)

# Draw CIs for each variable
for (i in seq(1.5,end_point,by=2)) {
  # J : Variable index (1~ncol(req))
  J <- (i+0.5) / 2
  
  # y : 1-sigma CI
  y <- c(mean(req[,J]) + sd(req[,J]), mean(req[,J]) - sd(req[,J]))
  
  # Draw the CI
  lines(x=c(i,i), y=y)
  lines(x=c(i-0.2,i+0.2), y=c(y[1],y[1])) 
  lines(x=c(i-0.2,i+0.2), y=c(y[2],y[2])) 
} 
par(mai=xx)       # Restore 'mai' values
