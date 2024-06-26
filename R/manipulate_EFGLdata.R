# basic manipulations of EFGL objects

#' combine multiple EFGL objects into one
#' @param genoComb if the objects have different loci, whether to create
#'   a new object with the intersection or union of loci. If union, genotypes
#'   for missing loci are all NA.
#' @param metaComb if the objects have different metadata fields, whether to create
#'   a new object with the intersection or union of the fields. If union, missing
#'   fields are all NA.
#' @param ... multiple EFGLdata objects separated by commas
#' @export
#'
combineEFGLdata <- function(..., genoComb = c("intersect", "union"),
									 metaComb = c("intersect", "union")){
	genoComb <- match.arg(genoComb)
	metaComb <- match.arg(metaComb)
	toComb <- list(...)
	for(i in toComb) if(class(i) != "EFGLdata") stop("All input objects must be of class EFGLdata")
	if(length(toComb) < 2){
		warning("Less than two EFGLdata objects input")
		return(toComb[[1]])
	}

	gD <- toComb[[1]]$genotypes
	mD <- toComb[[1]]$metadata
	for(i in 2:length(toComb)){
		# intersect or union of loci
		m1 <- c() # in case no genotypes
		m2 <- c()
		if(ncol(gD) > 2) m1 <- colnames(gD)[3:ncol(gD)]
		if(ncol(toComb[[i]]$genotypes) > 2) m2 <- colnames(toComb[[i]]$genotypes)[3:ncol(toComb[[i]]$genotypes)]
		if(genoComb == "intersect"){
			m_new <- m1[m1 %in% m2]
			gD <- gD %>% select(Pop, Ind, m_new)
			toComb[[i]]$genotypes <- toComb[[i]]$genotypes %>% select(Pop, Ind, m_new)
		} else {
			m_new <- c(m1, m2[!(m2 %in% m1)]) # to make sure order is valid
		}
		# combine
		gD <- bind_rows(gD, toComb[[i]]$genotypes) %>% select(Pop, Ind, m_new) # make sure order is valid

		# intersect or union of metadata
		d1 <- c() # in case no metadata
		d2 <- c()
		if(ncol(mD) > 2) d1 <- colnames(mD)[3:ncol(mD)]
		if(ncol(toComb[[i]]$metadata) > 2) d2 <- colnames(toComb[[i]]$metadata)[3:ncol(toComb[[i]]$metadata)]
		d_new <- d1[d1 %in% d2]
		if(metaComb == "intersect"){
			mD <- mD %>% select(Pop, Ind, d_new)
			toComb[[i]]$metadata <- toComb[[i]]$metadata %>% select(Pop, Ind, d_new)
		}
		# make sure variable types are the same
		for(c in d_new){
			if(class(pull(mD, c)) == class(pull(toComb[[i]]$metadata, c))) next
			if(any(!is.na(pull(mD, c))) && any(!is.na(pull(toComb[[i]]$metadata, c)))) stop("metadata variable type mismatch in ", c)
		}
		# combine
		mD <- bind_rows(mD, toComb[[i]]$metadata) %>% select(Pop, Ind, everything()) # make sure order is valid
	}

	return(construct_EFGLdata(list(genotypes = gD, metadata = mD)))

}

#' combine populations into one AND REMOVE the old populations
#' @param x an EFGLdata object
#' @param pops a vector of populations to combine
#' @param newName a string giving the name of the population to combine pops
#'   into. This can be a new pop or an existing pop (a warning is issued if existing).
#' @return an EFGLdata object
#' @export
#'
movePops <- function(x, pops, newName){
	if(length(newName) != 1) stop("newName must be one string")
	if(any(!pops %in% x$genotypes$Pop)) stop("not all pops are in this EFGLdata object")
	if(newName %in% x$genotypes$Pop) warning(newName, " already exists. Adding pops to an existing population.")

	x$genotypes$Pop[x$genotypes$Pop %in% pops] <- newName
	x$metadata$Pop[x$metadata$Pop %in% pops] <- newName

	return(construct_EFGLdata(x))
}

#' combine individuals into one population AND REMOVE the previous entry
#'   for those individuals
#' @param x an EFGLdata object
#' @param inds a vector of individuals to put in the new pop
#' @param newName a string giving the name of population to add the individuals
#'   too. This can be a new pop or an existing pop (a warning is issued if existing).
#' @return an EFGLdata object
#' @export
#'
moveInds <- function(x, inds, newName){
	if(length(newName) != 1) stop("newName must be one string")
	if(any(!inds %in% x$genotypes$Ind)) stop("not all inds are in this EFGLdata object")
	if(newName %in% x$genotypes$Pop) warning(newName, " already exists. Adding pops to an existing population.")

	x$genotypes$Pop[x$genotypes$Ind %in% inds] <- newName
	x$metadata$Pop[x$metadata$Ind %in% inds] <- newName

	return(construct_EFGLdata(x))
}

#' remove loci from an EFGLdata object
#' @param x an EFGLdata object
#' @param lociRemove a vector of loci names to remove
#' @export
removeLoci <-function(x, lociRemove){
	if(length(lociRemove) < 1) return(x) # if pass an empty vector, don't remove anything
	lociRemove <- c(paste0(lociRemove, ".A1"), paste0(lociRemove, ".A2"))
	if(any(!lociRemove %in% colnames(x$genotypes)[3:ncol(x$genotypes)])) stop("one or more loci were not found in input")
	x$genotypes <- x$genotypes %>% select(-lociRemove)
	return(x)
}

#' remove individuals from an EFGLdata object
#' @param x an EFGLdata object
#' @param inds a vector of individuals to remove
#' @return an EFGLdata object
#' @export
#'
removeInds <- function(x, inds){
	if(any(!inds %in% x$genotypes$Ind)) stop("not all inds are in this EFGLdata object")

	x$genotypes <- x$genotypes %>% filter(!(Ind %in% inds))
	x$metadata <- x$metadata %>% filter(!(Ind %in% inds))

	return(construct_EFGLdata(x))
}

#' remove pops from an EFGLdata object
#' @param x an EFGLdata object
#' @param pops a vector of pops to remove
#' @return an EFGLdata object
#' @export
#'
removePops <- function(x, pops){
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")

	x$genotypes <- x$genotypes %>% filter(!(Pop %in% pops))
	x$metadata <- x$metadata %>% filter(!(Pop %in% pops))

	return(construct_EFGLdata(x))
}

#' Identify which individual out of duplicate pairs has the lower genotyping success
#' @param dupTable the output of close_matching_samples (from rubias)
#' @param geno_success the output of genoSuccess
#' @return a vector of unique individual names representing the individuals with lower
#'   genotyping success from each pair
#' @export
#'
whichLower <- function(dupTable, geno_success){
	if(nrow(dupTable) < 1) return(c())
	dupTable <- dupTable %>%
		mutate(s1 = geno_success$success[match(indiv_1, geno_success$Ind)],
				s2 = geno_success$success[match(indiv_2, geno_success$Ind)],
				remove1 = s1 < s2)
	return(unique(c(dupTable$indiv_1[dupTable$remove1], dupTable$indiv_2[!dupTable$remove1])))
}
