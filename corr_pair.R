
##### ============================================================================================
##### ============================================================================================
##### ============== Calculate p-values for independence tests b/w variables ==================
##### ============================================================================================
##### ============================================================================================

pair.indep = function(input) {
	require(car)
	
	exp.vars = colnames(input)[colnames(input)!=resp.var]	## potential explanatory variables in 'input'
	factor_check = unlist(lapply(input,is.factor))			## boolean for factor variables in 'input'
	
	##### Test for and remove collinear variables
	corr.pairs = matrix(nc=3,nr=0)
	colnames(corr.pairs) = c("var1","var2","p.value")
		
	options(warn=0)
	
	for (var1.ind in 1:(length(exp.vars)-1)) {
		var1 = exp.vars[var1.ind]
		
		for (var2.ind in (var1.ind+1):length(exp.vars)) {
			var2 		 = exp.vars[var2.ind]
	
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

	### Remove variables causing p-value > 1
	# if (any(corr.pairs[,"p.value"] > 1)) {
	# 	over_tbl = table(corr.pairs[which( corr.pairs[,"p.value"] > 1 ),-3])
	# 	ranking = sort(rank(over_tbl,ties.method="min"),decreasing=T)
	# 	idx		= as.numeric(names(ranking)[which(ranking == ranking[1])])
	# 	corr.pairs = corr.pairs[!(corr.pairs[,1] %in% idx | corr.pairs[,2] %in% idx),]
	# }
	
	return(corr.pairs)
}