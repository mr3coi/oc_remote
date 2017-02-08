
##### ============================================================================================
##### ============================================================================================
##### =============== Utility Functions ===============
##### ============================================================================================
##### ============================================================================================

##### Date-to-int conversion function			##### TODO convert to # of days since 1900-01-01
date_conv = function(date_string) {
	date = strsplit(date_string,split=c("\\.|-"))
	date = as.numeric(unlist(date))
	
	### Return 'NA' if invalid 'date_string' given
	if (any(is.na(date))) return(NA)
	#date = sum(date*c(10000,100,1))
	
	leap = (date[1]-1901) %/% 4
	days = 365 * (date[1]-1901) + leap + date[3]
	this_leap = date[1] %% 4 == 0
	days_per_month = c(31,28,31,30,31,30,31,31,30,31,30,31)
	if (date[2] > 1) days = days + sum(days_per_month[1:(date[2]-1)])
	if (date[2] > 2 && this_leap) days = days + 1
	
	return(days)
}

##### Returns the column index in 'df' of the variable w/ specified 'name'
search_colname = function(df, name) { return(grep(name, colnames(df))) }

#####
extract_ind	= function(df,name) { 
	which(sapply(strsplit(colnames(df),"\\n"), function(v) trimws(v[1])) == name) }
