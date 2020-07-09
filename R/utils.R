# utility functions
# internal to the package

# Make a good guess is something is a biallelic SNP
isSNP <- function(a, b){
	a <- c(a, b)
	a <- unique(a[!is.na(a)])
	if(length(a) < 1) return(TRUE) # if all missing
	if(any(nchar(a) != 1)) return(FALSE)
	if(length(a) < 3) return(TRUE)
	return(FALSE)
}
