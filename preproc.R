
##### ============================================================================================
##### ============================================================================================
##### ============ datafile : Preprocess ==============
##### ============================================================================================
##### ============================================================================================

preProc = function(input, col_index, resp.var) {
	##### Check the input contains the response variable
	if (length(extract_ind(input,resp.var)) == 0) stop("Response variable not present in given data.")
	
	#### Extract the subset of data we are interested in
	dataset = input[,c(col_index,extract_ind(input,resp.var))]

	
	
	
	
	#####################################################################################
	##### Insert NA's where necessary
	#####################################################################################
	
	##### Change 'unknown' to NA
	NAval <- sapply(strsplit(colnames(dataset),"\\n"), 
					function(v) { q  <- trimws(v[-1])           # Generate character vector of items in colnames
					# w/ whitespaces removed
					qq <- grep("unknown", q, ignore.case=T)     # Find all indices in list w/ unknown items
					ifelse(length(qq),strsplit(q[qq[1]], "\\.")[[1]][1],NA) } ) 
	for (i in 1:ncol(dataset)) {
		if (!is.na(NAval[i])) dataset[which(NAval[i] == dataset[,i]),i] <- NA
	}
	
	##### Change "" to NA
	dataset[dataset==""] = NA

	
	
	
	
	
	
	
	#####################################################################################
	##### Removal of flawed columns
	#####################################################################################
		
	##### Specific deletions (flaws)
		### Delete unnamed columns & empty rows
	idx = which(colnames(dataset)=="")
	dataset = dataset[, -idx]
		### Delete the repeated patient number column
	useless = c("Pt_No.1","MEMO","Pt_No","Specimen_No","Name")
	dataset = dataset[,-which(colnames(dataset) %in% useless)]
	
	##### Find and remove columns w/ (NA or "")/total ratio > 0.1
	idx <- which(apply(dataset, 2, function(v) sum(is.na(v) | v==""))/nrow(dataset) > 0.1)
	dataset <- dataset[,-idx]
	
	##### Check that the response variable is still present in 'dataset'
	if (length(extract_ind(dataset,resp.var)) == 0) stop("Response variable has too many NA's.")
	
	
	
	
	
	
	
	#####################################################################################
	##### Manipulation for special cases
	#####################################################################################
	
	##### Find and modify(or remove) columns w/ non-numeric data
		### idx : non-numeric indices
	options(warn=2);	# to exploit errors
	idx <- which(apply(dataset, 2, function(v) class(try(as.numeric(v), T))=='try-error'))
	options(warn=1);
	
		### Specific cases
	# missed = c(56)
	# dataset[,missed] = lapply(dataset[,missed],as.numeric)
	# idx = idx[!idx %in% missed]
	
		### No-answer cases => exclude for now 						#### TODO fix this part
	no_answer = c(42,50,52)
	idx 	  = idx[!idx %in% no_answer]
	
		### Convert date info into numeric data (# days since 1901-01-01)
	pattern = "\\d+\\.\\d+\\.\\d+"
	date_vars = grep("(d|D)ate",colnames(dataset))		### Check which 'date' variables are surviving
	date_vars %in% idx
	
	idx2 = c()
	
	for (col in idx) {
		if (length(grep(pattern,dataset[,col])) != 0) {		### TODO Optimize
			dataset[,col] = sapply(dataset[,col], function(str) {ifelse(str != "", date_conv(str), NA)} )
			idx2 = append(idx2,col)
		}
	}
	date_vars %in% idx2 		### Check whether all surviving 'date' variables have been manipulated
	idx = idx[!idx %in% idx2]
	
		### Additional manipulation for 'Stage' column (transform "1a" ~ "4" => 1,2,3,4)
	stage.ind	= extract_ind(dataset,"Stage")
	if (length(grep("Stage",names(idx))) > 0) {
		ref 		= c(rep(1,4),rep(2,4),rep(3,4),4)
		names(ref)	= c("1","1a","1b","1c","2","2a","2b","2c","3","3a","3b","3c","4")
		dataset[,stage.ind] = ref[dataset[,stage.ind]]
		idx 		= idx[idx != stage.ind]
	}
	
		### Additional manipulation for 'Recur_Site' column
	# recur.ind = extract_ind(dataset,"Recur_Site")
	# dataset[,recur.ind] = sapply(dataset[,recur.ind],
	# 								function(v) {paste(sort(unlist(strsplit(v,split=","))),collapse="")})
	# dataset[,recur.ind]
	# idx = idx[idx != recur.ind]
	
		### Extract columns w/ multiple selections
		### (Assumption : multiple choices concatenated w/ commas and no spaces)
	pattern = "(\\d,)+\\d"
	mul_choice = c()
	for (col in idx) {
		if (length(grep(pattern,dataset[,col])) != 0) mul_choice = append(mul_choice, col)
	}
	if (!is.null(mul_choice)) idx = idx[idx != mul_choice]
	
		### Extract numeric values at the beginning
	pattern = "(\\d).+"
	idx2 = c()
	for (col in idx) {
		if (length(grep(pattern,dataset[,col])) != 0) {
			dataset[,col] = sub(pattern,"\\1",dataset[,col])
			idx2 = append(idx2,col)
		}
	}
	idx = idx[!idx %in% idx2]
	
		### Delete the remaining non-numeric columns
	if (length(idx) > 0) dataset = dataset[,-idx]
	
	
	
	
	
	#####################################################################################
	##### Conversion to numeric / factor (including dummy variable generation)
	#####################################################################################
	
	##### Check that the response variable is still included in the dataset
	if (length(extract_ind(dataset,resp.var)) == 0) stop("Response variable flawed.")
	
	##### Convert 'Refractory progression during CTx' to 'Yes'		### Obsolete
	# rec.ind = extract_ind(dataset,"Recurrence")
	# dataset[,rec.ind][!is.na(dataset[,rec.ind]) & dataset[,rec.ind]==2] <- 1
	
	##### Classify b/w categorical and numerical variables
	factor_check 		= sapply(strsplit(colnames(dataset),split="\n"), 
							function(v) ifelse(length(v) > 1,1,0))
	names(factor_check) = gsub(" ","_",sapply(strsplit(colnames(dataset),"\\n"), function(v) trimws(v[1])))
	
		### Add to 'factor_check' the variables that are categorical but not filtered	### Obsolete
	# manual_add			= c("Parity","Other_site")						# names of variables to add
	# factor_check[manual_add[manual_add %in% colnames(dataset)]] = 1
	
		### Convert the prespecified columns to factors
	for (i in 1:ncol(dataset)) { if (factor_check[i]==1) dataset[,i] = as.factor(dataset[,i]) }
	
	##### Choose main titles as column names (removing choices) and remove spaces within the main titles
	colnames(dataset) = names(factor_check)
	
	
	
	
	
	
	
	
	##### Dummy variable generation (including dealing w/ multiple choices)
	dummy.tmp = list()
	
	for (col in 1:ncol(dataset)) {
		if (factor_check[col] == 1 && nlevels(dataset[,col]) > 2) {
			cat("col : ", col, ", levels # : ", nlevels(dataset[,col]), "\n")
			
			### Generate dummy variables (w/ multiple-choices taking separate columns)
			dummy.mat = gen_dummy(dataset[,col,drop=F])
			#dummy.mat = model.matrix(~dataset[,col])[,-1]
			
			### Deal w/ multiple-choice cases
			if (col %in% mul_choice) {
				choices = levels(dataset[,col])[
							sapply(strsplit(levels(dataset[,col]),split=","), function(v) length(v) > 1)]
				
				### Check that all multiple-choice cases are not selected as the base factor for dummies
				if (all(match(choices,levels(dataset[,col])) > 1)) mul_choice = mul_choice[mul_choice != col]
				
				for (set in choices) {
					### items : vector of all values in the multiple choices
					items = unlist(strsplit(set,split=","))
					
					### Add column of 0's if any factor in the multiple-choice doesn't have a separate column
					for (miss in items[!items %in% levels(dataset[,col])]) {
						dummy.mat = cbind(dummy.mat, rep(0,nrow(dummy.mat)))
						colnames(dummy.mat)[ncol(dummy.mat)] = paste0("dataset[,col]",miss)
					}
					
					### index of the column falsely generated due to multiple choices
					false_column = grep(paste0(set,"$"),colnames(dummy.mat))
					
					### move each of the multiple choices to appropriate columns
					for (row in which(dummy.mat[,false_column] == 1)) {
						cat("Row w/ multiple choice : ", row, "\n")		### TODO DELETE
						for (item in items) {		### Add each item to corresponding column
							correct_column = grep(paste0(item,"$"),colnames(dummy.mat))
							dummy.mat[row,correct_column] = 1
						}
					}
					
					### remove the falsely-generated column
					dummy.mat = dummy.mat[,-false_column]
				}
			}
			
			### Final touches
				### Give appropriate names to each newly-generated columns ( (original_name)_(variable_value) )
			colnames(dummy.mat) = gsub("^dataset\\[, col\\](\\d+)$",
									   paste(colnames(dataset)[col],"\\1",sep="_"),colnames(dummy.mat))
			
				### Convert all dummy variable vectors into factors
			dummy.mat = apply(dummy.mat,2,as.character)
			dummy.mat = as.data.frame(dummy.mat,stringsAsFactors=T)
			
				### Store the dummy variables matrix to 'dummy.tmp' w/ 'col' as index
			dummy.tmp[[col]] = dummy.mat
		}
		else {		### Non-factors & factors w/ <= 2 levels
			dummy.tmp[[col]] = dataset[,col,drop=F]		### Just copy to 'dummy.tmp' for later merge
		}
	}
	
		### Check that 'dummy.tmp' has correctly been generated
	#lapply(dummy.tmp,head)
	
		### Check that the new dummy variables are all factors
	lapply(dummy.tmp, function(mat) apply(mat,2,class))
	
		
		### Check that all multiple-choice columns have been taken care of
	if (length(mul_choice) > 0) {
		stop(paste("Reorder factor levels so that multiple-choice level is not the base level",
				   " cols : ", mul_choice,sep=","))
	}
	
		### Merge 'dummy.tmp' into a data.frame that becomes the new 'dataset'
	dataset = as.data.frame(dummy.tmp)
		
		### Update 'factor_check'
	factor_check = ifelse(unlist(lapply(dataset,class)) == "factor",1,0)
	
	
	
	
	##### Convert character values in dataset into numeric values		### Obsolete
	# dataset <- as.data.frame(apply(dataset, 2, as.numeric))
	
	##### Check for flawed values in remaining columns
	# lv_data  = unlist(lapply(strsplit(colnames(dataset),split="\n"), length)) - 1
	# lv_count = unlist(lapply(dataset, function(v) length(unique(v[!is.na(v)]))))
	# names(lv_count) = NULL
	# lv_count[factor_check==0] = 0
	# 
	# names(factor_check)[which(lv_data - lv_count < 0)]
	# (lv_data - lv_count)[which(lv_data - lv_count < 0)]
	
	
	
	
	
	
	
	
	
	
	#####################################################################################
	##### Removal of unusable rows / columns
	#####################################################################################	
	
	##### Remove rows w/ rare values(< 0.01) in factor columns			### Obsolete
	# while(T) {
	# 	idx = c()
	# 	for (col in which(factor_check==1)) {
	# 		tmp = dataset[,col]
	# 		tbl = table(tmp) / length(tmp)
	# 		target = names(tbl)[which(tbl < 0.01)]
	# 		if (length(target) > 0) {
	# 			#cat(colnames(dataset)[col], '\n')
	# 			#cat(target, '---', which(tmp %in% target), '\n\n')
	# 			idx = c(idx,which(tmp %in% target))
	# 		}
	# 	}
	# 	idx = sort(unique(idx))
	# 	if (length(idx) == 0) break
	# 	dataset = dataset[-idx,] 
	# }

	##### Find and remove columns w/ NA/total ratio > 0.1
	idx = which(apply(dataset, 2, function(v) sum(is.na(v))/nrow(dataset)) > 0.1)
	if (length(idx) > 0) dataset <- dataset[,-idx]
	
	##### Find rows w/ NA/total ratio > 0.2
	idx = which(apply(dataset, 1, function(v) sum(is.na(v))/ncol(dataset)) > 0.2)
	if (length(idx) > 0) dataset = dataset[-idx,]
	
	return(dataset)
}