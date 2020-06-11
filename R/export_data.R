# functions to export data from an EFGLdata object in particular formats

#' export a rubias baseline
#' @param x an EFGLdata object
#' @param pops a vector of pops to include in the baseline. If not specified, all pops are used.
#' @param repunit the column name of the metadata variable designating repunit.
#'   This can be Pop to use the population name. If not specified, NA is used for all samples.
#' @param collection the column name of the metadata variable designating collection.
#'   This can be Pop to use the population name. If not specified, NA is used for all samples.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @return
#' @export
#'
exportRubias_baseline <- function(x, pops = NULL, repunit = NULL,
											 collection = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")
	if(any(!pops %in% x$genotypes$Pop)) stop("not all pops are in this EFGLdata object")

	# check that genotypes and metadata contain same pop and inds in order
	if(any(x$genotypes$Pop != x$metadata$Pop)) stop("Pop does not match between genotypes and metadata")
	if(any(x$genotypes$Ind != x$metadata$Ind)) stop("Ind does not match between genotypes and metadata")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	loci <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci %in% l)) stop("one or more loci were not found in input")
	# put in order

	l <- l[l %in% loci]
	rm(loci) # defensive
	# columns in order:
	# sample_type, repunit, collection, indiv, genotypes
	g <- x$genotypes %>% filter(Pop %in% pops) %>% select(Ind, l) %>% rename(indiv = Ind)
	x$metadata <- x$metadata %>% filter(Pop %in% pops)
	repunit <- if(is.null(repunit)) rep(NA, nrow(g)) else x$metadata %>% pull(repunit)
	collection <- if(is.null(collection)) rep(NA, nrow(g)) else x$metadata %>% pull(collection)
	return(tibble::tibble(sample_type = "reference", repunit = repunit,
							  collection = collection) %>% bind_cols(g))

}

#' export a rubias mixture
#' @param x an EFGLdata object
#' @param pops a vector of pops to include in the baseline. If not specified, all pops are used.
#' @param collection the column name of the metadata variable designating collection.
#'   This can be Pop to use the population name. If not specified, NA is used for all samples.
#'   For mixtures, this variable indicates what samples come from the same "stratum" - to be analyzed together.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @return
#' @export
#'
exportRubias_mixture <- function(x, pops = NULL, collection = NULL, loci = NULL){
	return(exportRubias_baseline(x, pops = pops, repunit = NULL,
										collection = collection, loci = loci) %>%
		mutate(sample_type = "mixture"))
}


