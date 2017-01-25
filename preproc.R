
##### ============================================================================================
##### ============================================================================================
##### ============ datafile : Preprocess ==============
##### ============================================================================================
##### ============================================================================================

preProc = function(dataset) {

	##### Change 'unknown' values to NA
	NAval <- sapply(strsplit(colnames(dataset),"\\n"), 
					function(v) { q  <- trimws(v[-1])           # Generate character vector of items in colnames
					# w/ whitespaces removed
					qq <- grep("unknown", q, ignore.case=T)     # Find all indices in list w/ unknown items
					ifelse(length(qq),strsplit(q[qq[1]], "\\.")[[1]][1],NA) } ) 
	for (i in 1:ncol(dataset)) {
		if (!is.na(NAval[i])) dataset[which(NAval[i] == dataset[,i]),i] <- NA
	}
	
	##### Specific deletions (flaws)
	### Delete unnamed columns & last two rows (empty rows)	(-2)
	idx = which(colnames(dataset)=="")
	dataset = dataset[ 1:236, -idx ]
	### Delete the repeated patient number column (-1)
	dataset = dataset[,-which(colnames(dataset)=="Pt_No.1")]
	
	##### Find and remove columns w/ (NA or "")/total ratio > 0.2 (-22)
	idx <- which(apply(dataset, 2, function(v) sum(is.na(v) | v==""))/nrow(dataset) > 0.2)
	dataset <- dataset[,-idx]
	
	##### Find and remove/modify columns w/ non-numeric data
		### idx : non-numeric indices
	options(warn=2);	# to exploit errors
	idx <- which(apply(dataset, 2, function(v) class(try(as.numeric(v), T))=='try-error'))
	options(warn=1);
	
		### Convert date info into numeric data (YYYYMMDD)
	pattern = "\\d+\\.\\d+\\.\\d+"
	idx2 = c()
	
	for (col in idx) {
		if (length(grep(pattern,dataset[,col])) != 0) {		### TODO Optimize
			dataset[,col] = sapply(dataset[,col], function(str) {ifelse(str != "", date_conv(str), NA)} )
			idx2 = c(idx2,col)
		}
	}
	idx = idx[!idx %in% idx2]
	
		### Additional manipulation for 'Stage' column (transform "1a" ~ "4" => 1,2,3,4)
	ref 		= c(rep(1,3),rep(2,3),rep(3,3),4)
	names(ref)	= c("1a","1b","1c","2a","2b","2c","3a","3b","3c","4")
	stage.ind	= extract_ind(dataset,"Stage")
	dataset[,stage.ind] = ref[dataset[,stage.ind]]
	idx 		= idx[idx != stage.ind]
	
		### Additional manipulation for 'Recur_Site' column
	# recur.ind = extract_ind(dataset,"Recur_Site")
	# dataset[,recur.ind] = sapply(dataset[,recur.ind],
	# 								function(v) {paste(sort(unlist(strsplit(v,split=","))),collapse="")})
	# dataset[,recur.ind]
	# idx = idx[idx != recur.ind]
	
		### Extract numeric values at the beginning
	pattern = "(\\d).*"
	idx2 = c()
	for (col in idx) {
		if (length(grep(pattern,dataset[,col])) != 0) {				### TODO Optimize
			dataset[,col] = sub(pattern,"\\1",dataset[,col])
			idx2 = c(idx2,col)
		}
	}
	idx = idx[!idx %in% idx2]
	
		### Delete the remaining non-numeric columns (-1)
	dataset <- dataset[,-idx]
	
	##### Convert character values in dataset into numeric values
	dataset <- as.data.frame(apply(dataset, 2, as.numeric))
	
	##### Convert 'Refractory progression during CTx' to 'Yes'
	rec.ind = extract_ind(dataset,"Recurrence")
	dataset[,rec.ind][!is.na(dataset[,rec.ind]) & dataset[,rec.ind]==2] <- 1
	
	##### Classify b/w categorical and numerical variables
	factor_check 		= sapply(strsplit(colnames(dataset),split="\n"), 
							function(v) ifelse(length(v) > 1,1,0))
	names(factor_check) = gsub(" ","_",sapply(strsplit(colnames(dataset),"\\n"), function(v) trimws(v[1])))
	
		### Add to 'factor_check' the variables that are categorical but not filtered
	manual_add			= c("Parity","Other_site")						# names of variables to add
	factor_check[manual_add[manual_add %in% colnames(dataset)]] = 1
	
		### Convert the prespecified columns to factors
	for (i in 1:ncol(dataset)) { if (factor_check[i]==1) dataset[,i] = as.factor(dataset[,i]) }
	
	##### Choose main titles as column names (removing choices) and remove spaces within the main titles
	colnames(dataset) = names(factor_check)
	
	##### Check for flawed values in remaining columns
	# lv_data  = unlist(lapply(strsplit(colnames(dataset),split="\n"), length)) - 1
	# lv_count = unlist(lapply(dataset, function(v) length(unique(v[!is.na(v)]))))
	# names(lv_count) = NULL
	# lv_count[factor_check==0] = 0
	# 
	# names(factor_check)[which(lv_data - lv_count < 0)]
	# (lv_data - lv_count)[which(lv_data - lv_count < 0)]
	
	# Menopausal_state : 2
	# Marital_state : 0
	# Educational_state : 0
	# Occupational_state : 2
	# Omentectomy : 2
	# 1st_Regimen : 0
	# Death : 4
	
	##### Remove rows w/ rare values(< 0.01) in factor columns			### TODO fix this part
	while(T) {
		idx = c()
		for (col in which(factor_check==1)) {
			tmp = dataset[,col]
			tbl = table(tmp) / length(tmp)
			target = names(tbl)[which(tbl < 0.01)]
			if (length(target) > 0) {
				#cat(colnames(dataset)[col], '\n')
				#cat(target, '---', which(tmp %in% target), '\n\n')
				idx = c(idx,which(tmp %in% target))
			}
		}
		idx = sort(unique(idx))
		if (length(idx) == 0) break
		dataset = dataset[-idx,] 
	}
	
	##### Find indices of columns w/ (almost : 95%) uniform values and remove such columns
	idx <- which( apply(dataset, 2, function(v) { any(table(v)/length(v) > 0.95) }) )
	dataset <- dataset[,-idx] 

	##### Find and remove columns w/ NA/total ratio > 0.2
	idx <- which(apply(dataset, 2, function(v) sum(is.na(v)))/nrow(dataset) > 0.2)
	if (length(idx) > 0) dataset <- dataset[,-idx]
	
	return(dataset)
}