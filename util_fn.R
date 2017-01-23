
##### ============================================================================================
##### ============================================================================================
##### =============== Utility Functions ===============
##### ============================================================================================
##### ============================================================================================

##### Date-to-int conversion function 			##### TODO convert to # of days since 1900-01-01
date_conv = function(date_string) {
	date = as.vector(sapply(strsplit(date_string,split="\\."),strtoi))
	date = sum(date*c(10000,100,1))
	return(date)
}

##### Returns the column index in 'df' of the variable w/ specified 'name'
search_colname = function(df, name) { return(grep(name, colnames(df))) }

#####
extract_ind	= function(df,name) { 
	which(sapply(strsplit(colnames(df),"\\n"), function(v) trimws(v[1])) == name) }
