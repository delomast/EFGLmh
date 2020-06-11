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

#' get a vector of loci names present
#' @param x an EFGLdata object
getLoci <- function(x){
	if(ncol(x$genotypes) < 3) stop("no genotypes")
	return(gsub("\\.A1$", "", colnames(x$genotypes)[seq(3,ncol(x$genotypes) - 1, 2)]))
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

#' calculate allelic richness of loci
#' @param x an EFGLdata object
#' @return a tibble giving the allelic richness of each locus
#' @export
aRich <- function(x){
	if(ncol(x$genotypes) < 3) stop("No loci found")
	return(
	x$genotypes %>% select(-Pop, -Ind) %>% tidyr::gather(locus, allele, 1:ncol(.)) %>%
		filter(!is.na(allele)) %>%
		mutate(locus = gsub("\\.A[12]$", "", locus)) %>% group_by(locus) %>%
		summarise(aRich = n_distinct(allele), .groups = "drop")
	)
}

#' calculate genotyping success of loci (uses only allele 1 for each
#'   genotype - assumes if allele 1 is (is not) NA, so is (is not) allele 2)
#' @param x an EFGLdata object
#' @return a tibble giving the genotyping success of each locus as a proportion
#' @export
lociSuccess <- function(x){
	if(ncol(x$genotypes) < 3) stop("No loci found")
	return(
	x$genotypes %>% select(-Pop, -Ind) %>% tidyr::gather(locus, allele, 1:ncol(.)) %>%
		filter(grepl("\\.A1$", locus)) %>%
		mutate(locus = gsub("\\.A1$", "", locus)) %>% group_by(locus) %>%
		summarise(success = sum(!is.na(allele)) / length(allele), .groups = "drop")
	)
}

#' calculate genotyping success of individuals (uses only allele 1 for each
#'   genotype - assumes if allele 1 is (is not) NA, so is (is not) allele 2)
#' @param x an EFGLdata object
#' @return a tibble giving the genotyping success of each individual as a proportion
#'   and number of missing genotypes
#' @export
genoSuccess <- function(x){
	if(ncol(x$genotypes) < 3) stop("No loci found")
	return(
		x$genotypes %>% tidyr::gather(locus, allele, 3:ncol(.)) %>%
			filter(grepl("\\.A1$", locus)) %>% group_by(Pop, Ind) %>%
			summarise(success = sum(!is.na(allele)) / length(allele),
						 numFail = sum(is.na(allele)), .groups = "drop")
	)
}

#' wrapper for write table with commonly used options - carried over from IDFGEN
#' @param x object to write out
#' @param filename filename to write out as
#' @param row.names passed to write.table
#' @param sep passed to write.table
#' @export
#'
dumpTable <- function(x, filename, row.names = FALSE, sep = "\t") {
	write.table(x, file = filename, append = FALSE, quote = FALSE, sep = sep,
					row.names = row.names, col.names = if(row.names) NA else TRUE)
}
