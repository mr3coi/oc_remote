
##### ============================================================================================
##### ============================================================================================
##### ============== Variable Selection : Stepwise AIC =============
##### ============================================================================================
##### ============================================================================================

stepwiseAIC = function(input,resp.var) {
	##### Preprocessing	
	exp.vars = colnames(input)[colnames(input)!=resp.var]
	
	##### Generate model using logistic regression
	f 	= as.formula(paste(resp.var, "~", paste(exp.vars, collapse=" + ")))
	mod = glm(f, data = input, family = "binomial")
	nul = glm(paste(resp.var, "~1"), data = input, family = "binomial")
	
	##### Conduct stepwise AIC both forwards and backwards
	step.0	<- stepAIC(nul, scope=list(lower=nul,upper=mod), direction="both")
	step.f	<- stepAIC(mod, direction="both")
	#step.0x <- step(nul, scope=list(lower=nul,upper=mod), direction="both")
	#step.fx <- step(mod, direction="both")
	
	# ##### Show ANOVA tables from the above regressions
	# step.0$anova
	# step.f$anova
	
	step0.coeff = ifelse(exp.vars %in% names(step.0$coefficients),1,0)
	stepF.coeff = ifelse(exp.vars %in% names(step.f$coefficients),1,0)
	
	##### Picked
	# var.full <- setdiff(exp_vars[[set_num]], substr(as.character(step.f$anova$Step), 3, 1000))
	# mean(do.cv(var.full, 3)$auc)
	# var.full <- substr(as.character(step.0$anova$Step), 3, 1000)[-1]
	# mean(do.cv(var.full, 3)$auc)
	# 
	# var.full <- setdiff(exp_vars[[set_num]], substr(as.character(step.fx$anova$Step), 3, 1000))
	# mean(do.cv(var.full, 3)$auc)
	# var.full <- substr(as.character(step.0x$anova$Step), 3, 1000)[-1]
	# mean(do.cv(var.full, 3)$auc)
	
	return(list(step0=step0.coeff,stepF=stepF.coeff))
}
