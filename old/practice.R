library(ROCR)
library(nnet)

wd <- "/Users/SJC/Documents/practice/internship"
setwd(wd)

main = readRDS("main.rds")
testset = "Residual_tumor_site_1st_debulking"
#sub <- main[!is.na(main[,testset]),]

avoid = search_colname(main,"Start_Date_1st_regimen"):length(colnames(main))
avoid = c(avoid, 1:3, 5)

size = 5
dep_var = sample(colnames(main)[-avoid],size)
final = na.omit(main[c(dep_var,testset)])
# check = apply(final,1,function(v) {any(is.na(v))})
# final = final[!check,]
# which(is.na(final))

final[,testset] = relevel(as.factor(final[,testset]), ref=order(unique(final[,testset]))[1])
curmod = multinom(testset ~ paste(dep_var,collapse="+"),data=final)

# ====================
	
search_colname = function(df, name) {
	return(grep(name, colnames(df)))
}