# functions relating to simulating data

#' Create F1 hybrids from two populations
#'
#' Creates hybrids by sampling one allele from one pop
#' and the other allele from the other pop
#' based on observed allele frequencies. Missing genotypes
#' are sampled based on the combined rate of missing genotypes
#' or data can be simulated wiht no missing genotypes. Alleles
#' at different loci are sampled independently. Note that if a sex marker exists
#' in the data set, it is treated like any other marker, resulting in XX, XY, and YY
#' genotypes in the hybrids. You may want to remove it prior to simulation using
#' the `removeLoci` function. Loci with all missing genotypes in one population
#' are simulated to have all missing genotypes in the hybrids REGARDLESS of
#' the input value for `missingGenos`.
#'
#' @param x an EFGLdata object
#' @param pop1 one parent population (must be in `x`)
#' @param pop2 the other parent population (must be in `x`)
#' @param newName a string giving the name of the population to put hybrids
#'   into. This MUST be a new pop.
#' @param n the number of hybrids to simulate
#' @param missingGenos TRUE to simulate missing genotypes, FALSE to not
#' @return an EFGLdata object with all the populations in `x`, plus the new simulated hybrid population
#' @export
createF1Hybrids <- function(x, pop1, pop2, newName, n = 50, missingGenos = TRUE){
	if(length(newName) != 1) stop("newName must be one string")
	if(any(!c(pop1, pop2) %in% x$genotypes$Pop)) stop("pop1 or pop2 are not in this EFGLdata object")
	if(newName %in% x$genotypes$Pop) stop(newName, " already exists")
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	# calculate allele frequencies
	af1 <- x %>% removePops(pops = getPops(x)[!getPops(x) %in% c(pop1, pop2)]) %>%
		calcAF()
	af2 <- af1 %>% filter(Pop == pop2)
	af1 <- af1 %>% filter(Pop == pop1)

	# calculate missing genotype frequency
	ls <- x %>% removePops(pops = getPops(x)[!getPops(x) %in% c(pop1, pop2)]) %>%
		lociSuccess() %>% mutate(miss = 1 - success)

	# set up for sampling
	af1_list <- list()
	af2_list <- list()
	# loci genotyped at least once in both pops
	loci <- intersect(unique(af1$locus), unique(af2$locus))
	for(i in 1:length(loci)){
		af1_list[[i]] <- af1 %>% filter(locus == loci[i])
		af2_list[[i]] <- af2 %>% filter(locus == loci[i])
	}
	miss <- data.frame(locus = loci) %>% left_join(ls, by = "locus") %>%
		pull(miss)

	# sample alleles
	g <- matrix(nrow = n, ncol = 0) # new genotypes
	for(i in 1:length(loci)){
		a1 <- sample(af1_list[[i]]$allele, size = n, prob = af1_list[[i]]$freq, replace = TRUE)
		a2 <- sample(af2_list[[i]]$allele, size = n, prob = af2_list[[i]]$freq, replace = TRUE)
		if(missingGenos){
			tempBool <- sample(c(TRUE, FALSE), size = n, prob = c(miss[i], 1 - miss[i]), replace = TRUE)
			a1[tempBool] <- NA
			a2[tempBool] <- NA
		}
		g <- cbind(g, a1, a2)
	}
	colnames(g) <- as.character(1:ncol(g)) # placeholder to initialize
	colnames(g)[seq(from = 1, to = (ncol(g) - 1), by = 2)] <- paste0(loci, ".A1")
	colnames(g)[seq(from = 2, to = ncol(g), by = 2)] <- paste0(loci, ".A2")
	g <- as.data.frame(g)
	g$Pop <- newName
	g$Ind <- paste0(newName, "_", 1:n)

	# add to EFGLdata object
	x$genotypes <- x$genotypes %>% bind_rows(g) # bind_rows should match column names and put NA for missing columns
	x$metadata <- x$metadata %>% bind_rows(g %>% select(Pop, Ind)) # bind_rows should put NA for missing columns

	return(construct_EFGLdata(x))
}
