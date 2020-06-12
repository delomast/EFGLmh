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
#' @return a tibble
#' @export
#'
exportRubias_baseline <- function(x, pops = NULL, repunit = NULL,
											 collection = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	# check that genotypes and metadata contain same pop and inds in order
	if(any(x$genotypes$Pop != x$metadata$Pop)) stop("Pop does not match between genotypes and metadata")
	if(any(x$genotypes$Ind != x$metadata$Ind)) stop("Ind does not match between genotypes and metadata")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% x$genotypes$Pop)) stop("not all pops are in this EFGLdata object")
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
#' @return a tibble
#' @export
#'
exportRubias_mixture <- function(x, pops = NULL, collection = NULL, loci = NULL){
	return(exportRubias_baseline(x, pops = pops, repunit = NULL,
										collection = collection, loci = loci) %>%
		mutate(sample_type = "mixture"))
}

#' write a genepop input file
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param header a string to use as the header line of the genepop file
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @param useIndNames TRUE to use individual names as sample identifiers. Otherwise,
#'   population names are used
#' @param return nothing, just writes a file
#' @export
#'
exportGenoPop <- function(x, filename, header = "genePop file", pops = NULL,
								  loci = NULL, useIndNames = FALSE){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% x$genotypes$Pop)) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	# select only pop and loci
	g <- x$genotypes %>% filter(Pop %in% pops) %>% select(Pop, Ind, l)
	if(useIndNames){
		gp <- g %>% select(Pop, Ind) %>% mutate(Ind = paste0(Ind, ","))
	} else {
		gp <- g %>% select(Pop) %>% mutate(Ind = paste0(Pop, ","))
	}
	# convert genotypes to genepop formatting
	to_remove <- c()
	for(i in loci){
		a <- g %>% select(paste0(i, ".A1"), paste0(i, ".A2"))
		alleleIndex <- a %>% tidyr::gather(locus, allele, 1:2) %>%
			filter(!is.na(allele)) %>% pull(allele) %>% unique
		if(length(alleleIndex) < 1){
			warning("all missing data for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		}
		if(length(alleleIndex) > 99) stop("More than 99 alleles at locus ", i)
		alleleIndex <- tibble(allele = alleleIndex,
									 index = gsub(" ", "0", format(1:length(alleleIndex), width = 2)))
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		a1[is.na(a1)] <- "00"
		a2[is.na(a2)] <- "00"
		gp <- gp %>% tibble::add_column(!!i := paste0(a1,a2))
	}
	loci <- loci[!loci %in% to_remove]

	#write genepop file
	cat(header, "\n", file = filename, append = FALSE, sep = "")
	for(i in loci) cat(i, "\n", sep = "", file = filename, append = TRUE)
	for(p in pops){
		cat("POP\n", file = filename, append = TRUE)
		write.table(filter(gp, Pop == p) %>% select(-Pop), file = filename,
						sep = " ", append = TRUE, col.names = FALSE,
						row.names = FALSE, quote = FALSE)
	}
}


