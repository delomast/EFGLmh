# data extraction from EFGLdata
# yes, these are simple to do be just referencing tibbles,
# but these functions make it easy for non-R afficionados to
# get data and perform analyses

#' get a vector of populations (pedigrees) present
#' @param x an EFGLdata object
#' @return a vector of the unique population names present
#' @export
getPops <- function(x){
	return(unique(x$genotypes$Pop))
}

#' get a vector of individuals present
#' @param x an EFGLdata object
#' @param pops a vector of pops that you want individual names for. If not
#'   specified, names for all pops are returned
#' @return a vector of the Individual names present
#' @export
getInds <- function(x, pops = NULL){
	if(is.null(pops)) return(x$genotypes$Ind)
	return(x$genotypes %>% filter(Pop %in% pops) %>% pull(Ind))
}

#' get the number of individuals present in each pop
#' @param x an EFGLdata object
#' @param pops a vector of pops that you want individual names for. If not
#'   specified, names for all pops are returned
#' @return a named vector with the number of individuals in each pop
#' @export
numInds <- function(x, pops = NULL){
	if(is.null(pops)) pops <- unique(x$genotypes$Pop)
	x <- x$genotypes %>% filter(Pop %in% pops) %>% count(Pop)
	# return as a vector
	y <- x$n
	names(y) <- x$Pop
	return(y)
}
