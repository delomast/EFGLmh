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
combineEFGLdata <- function(genoComb = c("intersect", "union"),
									 metaComb = c("intersect", "union"), ...){
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
		} else {
			m_new <- c(m1, m2[!(m2 %in% m1)])
			# add new columns if needed
			# first to existing
			toAdd <- m_new[!(m_new %in% m1)]
			for(j in toAdd) gD <- gD %>% tibble::add_column(!!j := NA)
			# then to new
			toAdd <- m_new[!(m_new %in% m2)]
			for(j in toAdd) toComb[[i]]$genotypes <- toComb[[i]]$genotypes %>% tibble::add_column(!!j := NA)
		}
		# combine
		gD <- gD %>% select(Pop, Ind, m_new)
		toComb[[i]]$genotypes <- toComb[[i]]$genotypes %>% select(Pop, Ind, m_new)
		gD <- bind_rows(gD, toComb[[i]]$genotypes)

		# intersect or union of metadata
		m1 <- c() # in case no metadata
		m2 <- c()
		if(ncol(mD) > 2) d1 <- colnames(mD)[3:ncol(mD)]
		if(ncol(toComb[[i]]$metadata) > 2) d2 <- colnames(toComb[[i]]$metadata)[3:ncol(toComb[[i]]$metadata)]
		if(genoComb == "intersect"){
			d_new <- d1[d1 %in% d2]
		} else {
			d_new <- c(d1, d2[!(d2 %in% d1)])
			# add new columns if needed
			# first to existing
			toAdd <- d_new[!(d_new %in% d1)]
			for(j in toAdd) mD <- mD %>% tibble::add_column(!!j := NA)
			# then to new
			toAdd <- d_new[!(d_new %in% d2)]
			for(j in toAdd) toComb[[i]]$metadata <- toComb[[i]]$metadata %>% tibble::add_column(!!j := NA)

		}
		# combine
		mD <- mD %>% select(Pop, Ind, d_new)
		toComb[[i]]$metadata <- toComb[[i]]$metadata %>% select(Pop, Ind, d_new)
		mD <- bind_rows(mD, toComb[[i]]$metadata)
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
removeLoci <-function(x, lociRemove){
	lociRemove <- c(paste0(lociRemove, ".A1"), paste0(lociRemove, ".A2"))
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
