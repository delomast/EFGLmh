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
		m1 <- colnames(gD)[3:ncol(gD)]
		m2 <- colnames(toComb[[i]]$genotypes)[3:ncol(toComb[[i]]$genotypes)]
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
		d1 <- colnames(mD)[3:ncol(mD)]
		d2 <- colnames(toComb[[i]]$metadata)[3:ncol(toComb[[i]]$metadata)]
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
		mD <- bind_rows(mD, toComb[[i]]$genotypes)
	}

	return(construct_EFGLdata(list(genotypes = gD, metadata = mD)))

}

