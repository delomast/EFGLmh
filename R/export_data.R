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
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
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
	repunit <- if(is.null(repunit)) rep(NA, nrow(g)) else x$metadata[[repunit]]
	collection <- if(is.null(collection)) rep(NA, nrow(g)) else x$metadata[[collection]]
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
#' @return nothing, just writes a file
#' @export
#'
exportGenePop <- function(x, filename, header = "genePop file", pops = NULL,
								  loci = NULL, useIndNames = FALSE){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops)
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

#' export a gRandma baseline or mixture
#' @param x an EFGLdata object
#' @param pops a vector of pops to include in the baseline. If not specified, all pops are used.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @param baseline TRUE to make a baseline input, FALSE to make a mixture input.
#' @return a tibble
#' @export
#'
exportGrandma <- function(x, pops = NULL, loci = NULL, baseline = TRUE){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci %in% l)) stop("one or more loci were not found in input")
	# put in order
	l <- l[l %in% loci]
	rm(loci) # defensive
	# columns in order:
	# Pop, Ind, Loci (no Pop for mixture)
	g <- x$genotypes %>% filter(Pop %in% pops) %>% select(Pop, Ind, l)
	if(!baseline) g <- g %>% select(-Pop)

	return(g)
}

#' convenience function to remove loci with all fails or no variation from
#' gRandma input
#' @param baseline a gRandma baseline input
#' @param mixture a gRandma mixture input
#' @return a list with two components, one is the baseline, one is the mixture
#' @export
cleanGrandma <- function(baseline, mixture = NULL){
	colnames(baseline)[1] <- "Pop"
	if (is.null(mixture)){
		g <- baseline
	} else {
		# make sure loci colnames are the same
		if(any(colnames(baseline)[3:ncol(baseline)] !=
				 colnames(mixture)[2:ncol(mixture)])) stop("loci colnames do not match")
		colnames(mixture)[1] <- colnames(baseline)[2]
		g <- tibble::tibble(Pop = rep(NA, nrow(mixture))) %>%
			bind_cols(mixture) %>% bind_rows(baseline)
	}
	loci <- gsub("\\.A1$", "", colnames(g)[seq(3, ncol(g) - 1, 2)])
	to_remove <- c()
	for(l in loci){
		a <- g %>% select(paste0(l, ".A1"), paste0(l, ".A2")) %>%
			tidyr::gather(locus, allele, 1:2) %>%
			filter(!is.na(allele)) %>% pull(allele) %>% unique
		if(length(a) < 1){
			cat("Removing locus", l, "for all missing genotypes.\n", sep = " ")
			to_remove <- c(to_remove, l)
		} else if(length(a) < 2){
			cat("Removing locus", l, "for no variation.\n", sep = " ")
			to_remove <- c(to_remove, l)
		}
	}
	if(length(to_remove) > 0){
		to_remove <- c(paste0(to_remove, ".A1"), paste0(to_remove, ".A2"))
		g <- g %>% select(-to_remove)
	}
	if(is.null(mixture)){
		g <- list(baseline = g %>% filter(!is.na(Pop)),
					 mixture = NULL)
	} else {
		g <- list(baseline = g %>% filter(!is.na(Pop)),
					 mixture = g %>% filter(is.na(Pop)) %>% select(-Pop))
	}
	return(g)
}

#' write a GenAlEx input file
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @param title a string to use as the "title" row
#' @param useNames TRUE to use samples names, FALSE to replace with unique numerical identifiers
#' @return nothing, just writes a file
#' @export
#'
exportGenAlEx <- function(x, filename, pops = NULL, loci = NULL, title = "", useNames = TRUE){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops) %>% arrange(Pop)
	gp <- g %>% select(Ind, Pop)
	if(!useNames) gp <- gp %>% mutate(Ind = 1:nrow(.))

	# convert genotypes to numerical formatting
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
		alleleIndex <- tibble(allele = alleleIndex,
									 index = 1:length(alleleIndex))
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		a1[is.na(a1)] <- 0
		a2[is.na(a2)] <- 0
		gp <- gp %>% tibble::add_column(!!i := a1,
												  !!paste0(i, ".A2") := a2)
	}
	loci <- loci[!loci %in% to_remove]

	# count number of inds in each pop, in order
	popCount <- c()
	for(p in pops) popCount <- c(popCount, sum(gp$Pop == p))

	#write GenAlEx file
	# num loci, num samples, num pops, num samps in each pop
	cat((ncol(gp) - 2)/2, nrow(gp), length(pops), popCount, file = filename, append = FALSE, sep = "\t")
	cat("\n", file = filename, append = TRUE)
	cat(title, "\t\t\t", file = filename, append = TRUE, sep = "")
	cat(pops, file = filename, append = TRUE, sep = "\t")
	cat("\n", file = filename, append = TRUE)
	lNames <- colnames(gp)[3:ncol(gp)]
	lNames[grepl("\\.A2$", lNames)] <- ""
	cat("SampNumber", "Pop", lNames, file = filename, append = TRUE, sep = "\t")
	cat("\n", file = filename, append = TRUE)
	for(p in pops){
		write.table(filter(gp, Pop == p), file = filename,
						sep = "\t", append = TRUE, col.names = FALSE,
						row.names = FALSE, quote = FALSE)
	}
}

#' write a Structure input file.
#'
#' Columns are Ind, Pop (converted to integers), then genotypes. Missing alleles
#' are coded as -9
#'
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @return nothing, just writes a file
#' @export
exportStructure <- function(x, filename, pops = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops) %>% arrange(Pop) %>%
		mutate(Pop = as.numeric(as.factor(Pop))) # turn pops into integers
	gp <- g %>% select(Ind, Pop)

	# convert genotypes to numerical formatting
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
		alleleIndex <- tibble(allele = alleleIndex,
									 index = 1:length(alleleIndex))
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		a1[is.na(a1)] <- -9
		a2[is.na(a2)] <- -9
		gp <- gp %>% tibble::add_column(!!i := a1,
												  !!paste0(i, ".A2") := a2)
	}
	loci <- loci[!loci %in% to_remove]

	cat(loci, file = filename, sep = "\t", append = FALSE) # row of locus names
	cat("\n", file = filename, append = TRUE)
	# Ind, Pop, genotypes
	write.table(gp, file = filename, sep = "\t", append = TRUE,
					col.names = FALSE, row.names = FALSE, quote = FALSE)
}

#' return a CKMRsim allele frequency input tibble.
#'
#' No chromosome and position information is used.
#' All populations specified are combined into one allele frequency output.
#' Any loci with only missing genotypes are omitted from the output. The output
#' should then be run through the CKMRsim function \code{reindex_markers()}.
#'
#' @param x an EFGLdata object
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @return a tibble
#' @export
exportCKMRsimAF <- function(x, pops = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops) %>% select(-Pop) %>%
		gather(Locus, Allele, 2:ncol(.)) %>%
		mutate(Locus = gsub("\\.A[12]$", "", Locus)) %>%
		filter(Locus %in% loci, !is.na(Allele)) %>%
		group_by(Locus) %>% count(Allele) %>%
		mutate(Freq = n / sum(n)) %>% ungroup() %>%
		mutate(Chrom = "Unk", Pos = 1, LocIdx = NA, AlleIdx = NA) %>%
		select(Chrom, Locus, Pos, Allele, LocIdx, AlleIdx, Freq)

	locRemoved <- loci[!loci %in% g$Locus]
	if(length(locRemoved) > 0) warning("Removed ", paste(locRemoved, collapse = " "), " for all missing data.")

	# remove spaces from marker names if applicable
	if(any(grepl(" ", loci))){
		warning("Removing spaces from locus names in output.")
		if(length(unique(gsub(" ", "", loci))) != length(unique(loci))) stop("Marker names are not unique after removal of spaces.")
		g$Locus <- gsub(" ", "", g$Locus)
	}

	return(g)
}

#' Return a long genotype tibble for input to CKMRsim for relationship inference.
#'
#' @param x an EFGLdata object
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @return a tibble
#' @export
exportCKMRsimLG <- function(x, pops = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops) %>% select(-Pop) %>%
		gather(Locus, Allele, 2:ncol(.)) %>%
		mutate(gene_copy = ifelse(grepl("\\.A1$", Locus), 1, 2),
				 Locus = gsub("\\.A[12]$", "", Locus)) %>%
		filter(Locus %in% loci) %>% rename(Indiv = Ind) %>%
	select(Indiv, Locus, gene_copy, Allele) %>% arrange(Indiv, Locus, gene_copy)

	# remove spaces from marker names if applicable
	if(any(grepl(" ", loci))){
		warning("Removing spaces from locus names in output.")
		if(length(unique(gsub(" ", "", loci))) != length(unique(loci))) stop("Marker names are not unique after removal of spaces.")
		g$Locus <- gsub(" ", "", g$Locus)
	}

	return(g)
}


#' write a SNPPIT input file. Will warn about skipping loci with > 2 alleles.
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param baseline a vector of pops to use as the baseline (potential parents).
#' @param mixture a vector of pops to use as the mixture (potential offspring).
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @param errorRate per allele error rate for all loci
#' @param POPCOLUMN_SEX metadata column with sex info (coded as M, F, and ?)
#' @param POPCOLUMN_REPRO_YEARS metadata column with repro years
#' @param POPCOLUMN_SPAWN_GROUP metadata column with spawn group
#' @param OFFSPRINGCOLUMN_BORN_YEAR metadata column with birth year
#' @param OFFSRPINGCOLUMN_SAMPLE_YEAR metadata column with sample year
#' @param OFFSPRINGCOLUMN_AGE_AT_SAMPLING metadata column with age at sampling
#' @return nothing, just writes a file
#' @export
exportSNPPIT <- function(x, filename,
								 baseline,
                      mixture,
                      loci = NULL,
                      errorRate = .005,
                      POPCOLUMN_SEX = NULL,
                      POPCOLUMN_REPRO_YEARS = NULL,
                      POPCOLUMN_SPAWN_GROUP = NULL,
                      OFFSPRINGCOLUMN_BORN_YEAR = NULL,
                      OFFSRPINGCOLUMN_SAMPLE_YEAR = NULL,
                      OFFSPRINGCOLUMN_AGE_AT_SAMPLING = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	pops <- c(baseline, mixture)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")

	g <- x$genotypes %>% filter(Pop %in% pops) %>% arrange(Pop)
	gp <- g %>% select(Pop, Ind)
	# convert genotypes to numerical formatting
	to_remove <- c()
	for(i in loci){
		a <- g %>% select(paste0(i, ".A1"), paste0(i, ".A2"))
		alleleIndex <- a %>% tidyr::gather(locus, allele, 1:2) %>%
			filter(!is.na(allele)) %>% pull(allele) %>% unique
		if(length(alleleIndex) < 1){
			warning("all missing data for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		} else if (length(alleleIndex) > 2){
			warning("More than two alleles found for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		}

		alleleIndex <- tibble(allele = alleleIndex,
									 index = 1:length(alleleIndex))
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		a1[is.na(a1)] <- 0
		a2[is.na(a2)] <- 0
		gp <- gp %>% tibble::add_column(!!i := a1,
												  !!paste0(i, ".A2") := a2)
	}
	loci <- loci[!loci %in% to_remove]

	# write header
	cat("NUMLOCI ", length(loci), "\n", file = filename, append = FALSE, sep = "")
	cat("MISSING_ALLELE 0", "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(POPCOLUMN_SEX)) cat("POPCOLUMN_SEX",  "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(POPCOLUMN_REPRO_YEARS)) cat("POPCOLUMN_REPRO_YEARS", "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(POPCOLUMN_SPAWN_GROUP)) cat("POPCOLUMN_SPAWN_GROUP", "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(OFFSPRINGCOLUMN_BORN_YEAR)) cat("OFFSPRINGCOLUMN_BORN_YEAR", "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(OFFSRPINGCOLUMN_SAMPLE_YEAR)) cat("OFFSRPINGCOLUMN_SAMPLE_YEAR", "\n", file = filename, append = TRUE, sep = "")
	if(!is.null(OFFSPRINGCOLUMN_AGE_AT_SAMPLING)) cat("OFFSPRINGCOLUMN_AGE_AT_SAMPLING", "\n", file = filename, append = TRUE, sep = "")

	# write error rates
	write.table(data.frame(l = loci, e = errorRate), file = filename, append = TRUE,
					quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

	# write baseline pops
	for(p in baseline){
		cat("POP ", p, "\n", file = filename, append = TRUE, sep = "")

		t <- gp %>% filter(Pop == p) %>% select(-Pop)
		# add metadata if needed
		metaCols <- c()
		if(!is.null(POPCOLUMN_SEX)) metaCols <- c(metaCols, POPCOLUMN_SEX)
		if(!is.null(POPCOLUMN_REPRO_YEARS)) metaCols <- c(metaCols, POPCOLUMN_REPRO_YEARS)
		if(!is.null(POPCOLUMN_SPAWN_GROUP)) metaCols <- c(metaCols, POPCOLUMN_SPAWN_GROUP)
		if(length(metaCols) > 0){
			m <- x$metadata %>% filter(Pop = p) %>% select(Ind, metaCols)
			t <- t %>% left_join(m, by = "Ind") %>% select(Ind, metaCols, everything())
		}
		write.table(t, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
	}
	warnMeta <- FALSE
	if(length(metaCols) > 0) warnMeta <- TRUE

	# write mixture pops
	for(p in mixture){
		cat("OFFSPRING ", p, " ?", "\n", file = filename, append = TRUE, sep = "")
		t <- gp %>% filter(Pop == p) %>% select(-Pop)
		# add metadata if needed
		metaCols <- c()
		if(!is.null(OFFSPRINGCOLUMN_BORN_YEAR)) metaCols <- c(metaCols, OFFSPRINGCOLUMN_BORN_YEAR)
		if(!is.null(OFFSRPINGCOLUMN_SAMPLE_YEAR)) metaCols <- c(metaCols, OFFSRPINGCOLUMN_SAMPLE_YEAR)
		if(!is.null(OFFSPRINGCOLUMN_AGE_AT_SAMPLING)) metaCols <- c(metaCols, OFFSPRINGCOLUMN_AGE_AT_SAMPLING)
		if(length(metaCols) > 0){
			m <- x$metadata %>% filter(Pop = p) %>% select(Ind, metaCols)
			t <- t %>% left_join(m, by = "Ind") %>% select(Ind, metaCols, everything())
		}
		write.table(t, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
	}
	if(length(metaCols) > 0) warnMeta <- TRUE
	if(warnMeta) warning("Options for use of metadata have not been thoroughly tested. Please report errors and use at your own risk.")
}


#' export a hierfstat input dataframe
#' @param x an EFGLdata object
#' @param pops a vector of pops to include. If not specified, all pops are used.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @return a dataframe coded to be used as input for hierfstat
#' @export
#'
exportHierFstat <- function(x, pops = NULL, loci = NULL){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	if(any(!(c(paste0(loci, ".A1"), paste0(loci, ".A2")) %in%
				colnames(x$genotypes)[3:ncol(x$genotypes)]))) stop("one or more loci were not found in input")
	# translate to one col per call
	g <- x$genotypes %>% filter(Pop %in% pops) %>% arrange(Pop)
	gp <- g %>% select(Pop) # columns: Pop, genotypes (one col per locus)

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
		if(length(alleleIndex) > 90) stop("More than 90 alleles at locus ", i)
		alleleIndex <- tibble(allele = alleleIndex,
									 index = 10:(9 + length(alleleIndex))
							)
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		aC <- paste0(a1,a2) # deal with NA's to avoid warnings from as.numeric
		aC[aC == "NANA"] <- NA
		gp <- gp %>% tibble::add_column(!!i := as.numeric(aC))
	}

	return(as.data.frame(gp))
}

#' export a "Progeny-style" export file for later reading into EFGLmh
#'
#' Columns in order are Pop, Ind, metadata, genotypes (2-column per call)
#' Missing genotypes are "0" for SNPs (biallelic or nonvariable with alleles represented by 1 character)
#' and "000" for others. If a locus is all missing, it is treated as a SNP. Pop column is called
#' "Pedigree" and Ind column is called "Individual.Name".
#'
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param pops a vector of pops to include. If not specified, all pops are used.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @param metadata a vector of metadata fields to include. If not specified,
#'   all fields are used.
#' @return nothing, just writes a file
#' @export
#'
exportProgenyStyle <- function(x, filename, pops = NULL, loci = NULL, metadata = NULL){
	x <- construct_EFGLdata(x) # basic checks on structure and order
	if(ncol(x$genotypes) < 3) stop("no genotypes")
	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(is.null(metadata)){
		metadata <- colnames(x$metadata)
		metadata <- metadata[!(metadata %in% c("Pop", "Ind"))]
	}
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	if(any(!(c(paste0(loci, ".A1"), paste0(loci, ".A2")) %in%
				colnames(x$genotypes)[3:ncol(x$genotypes)]))) stop("one or more loci were not found in input")

	# filter loci
	# remove not included
	toRem <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	toRem <- toRem[!(toRem %in% c(paste0(loci, ".A1"), paste0(loci, ".A2")))]
	x$genotypes <- x$genotypes %>% select(-toRem)
	# change missing genotype symbols
	snpNames <- c()
	otherNames <-c()
	for(l in loci){
		lnames <- c(paste0(l, ".A1"), paste0(l, ".A2"))
		if(isSNP(pull(x$genotypes, lnames[1]), pull(x$genotypes, lnames[2]))){
			snpNames <- c(snpNames, lnames)
		} else {
			otherNames <- c(otherNames, lnames)
		}
	}
	if(length(snpNames) > 0) x$genotypes <- x$genotypes %>% mutate_at(snpNames, tidyr::replace_na, replace = "0")
	if(length(otherNames) > 0) x$genotypes <- x$genotypes %>% mutate_at(otherNames, tidyr::replace_na, replace = "000")

	# combine, select metadata columns, order columns, and write out
	gNames <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	x$metadata %>% left_join(x$genotypes, by = c("Pop", "Ind")) %>%
		select(Pop, Ind, metadata, gNames) %>%
		rename(Pedigree = Pop, Individual.Name = Ind) %>% dumpTable(filename = filename)

}

#' export a Plink PED file and optionally a (not very useful except as a template) MAP file
#'
#' Only biallelic markers are used, exports as two column per call. This format can be used
#' as the input for "Admixture".
#'
#' @param x an EFGLdata object
#' @param filename the name of the file to write
#' @param pops a vector of pops to include. If not specified, all pops are used.
#' @param loci a vector of loci to use. If not specified,
#'   all loci are used.
#' @param FID Metadata column name to use as the family ID
#' @param IID Metadata column name to use as the within-family ID
#' @param pa Metadata column name to use as the within-family ID of the father.
#'   If NULL, all are listed as unknown (0).
#' @param ma Metadata column name to use as the within-family ID of the mother.
#'   If NULL, all are listed as unknown (0).
#' @param sex Metadata column to use as the sex. This should be coded as "M"
#'   for male, "F" for female, and any other values are treated as Unknown.
#'   If NULL, all are listed as unknown (0).
#' @param pheno Metadata column to use as the phenotype. This should be coded
#'   as "1" for control, "2" for case, and "-9" or "0" for Unknown. If NULL,
#'   all are listed as unknown (0).
#' @param map filename to write a MAP file. If NULL, no MAP file is written.
#'   This just writes a dummy MAP file with the loci names in order and all
#'   given the same chromosome code and positions of "0".
#'   If you need a valid MAP file, you will need to edit this.
#' @return nothing, just writes a file
#' @export
#'
exportPlink <- function(x, filename, pops = NULL, loci = NULL,
								FID = "Ind", IID = "Ind", pa = NULL, ma = NULL,
								sex = NULL, pheno = NULL,
								map = NULL){
	x <- construct_EFGLdata(x) # basic checks on structure and order
	if(ncol(x$genotypes) < 3) stop("no genotypes")
	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)

	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	if(any(!(c(paste0(loci, ".A1"), paste0(loci, ".A2")) %in%
				colnames(x$genotypes)[3:ncol(x$genotypes)]))) stop("one or more loci were not found in input")

	meta <- x$metadata %>% filter(Pop %in% pops) %>% arrange(Pop, Ind)
	g <- x$genotypes %>% filter(Pop %in% pops) %>% arrange(Pop, Ind)
	if(any(meta$Pop != g$Pop) || any(meta$Ind != g$Ind)) stop("Internal error in ordering") # just to make sure and ease my mind
	gp <- g %>% select(Ind)

	# filter out unusable loci
	to_remove <- c()
	for(i in loci){
		a <- g %>% select(paste0(i, ".A1"), paste0(i, ".A2"))
		alleleIndex <- a %>% tidyr::gather(locus, allele, 1:2) %>%
			filter(!is.na(allele)) %>% pull(allele) %>% unique
		if(length(alleleIndex) < 1){
			warning("all missing data for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		} else if (length(alleleIndex) > 2){
			warning("More than two alleles found for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		} else if (length(alleleIndex) < 2){
			warning("Fewer than two alleles found for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		}

		# convert to 1/2 format for consistency
		alleleIndex <- tibble(allele = alleleIndex,
									 index = 1:length(alleleIndex))
		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]

		a1[is.na(a1)] <- "0"
		a2[is.na(a2)] <- "0"
		gp <- gp %>% tibble::add_column(!!i := a1,
												  !!paste0(i, ".A2") := a2)
	}
	rm(g)
	loci <- loci[!loci %in% to_remove]

	# build ped file
	ped <- tibble::tibble(FID = meta[[FID]], IID = meta[[IID]])
	if(is.null(pa)){
		ped <- ped %>% tibble::add_column(paID = "0")
	} else {
		ped <- bind_cols(ped, meta[[pa]])
	}
	if(is.null(ma)){
		ped <- ped %>% tibble::add_column(maID = "0")
	} else {
		ped <- bind_cols(ped, meta[[ma]])
	}
	if(is.null(sex)){
		ped <- ped %>% tibble::add_column(sex = "0")
	} else {
		meta[[sex]] <- ifelse(meta[[sex]] == "M", "1", ifelse(meta[[sex]] == "F", "2", "0"))
		ped <- bind_cols(ped, meta[[sex]])
	}
	if(is.null(pheno)){
		ped <- ped %>% tibble::add_column(pheno = "0")
	} else {
		ped <- bind_cols(ped, meta[[pheno]])
	}
	gp <- gp %>% select(-Ind)
	gp <- bind_cols(ped, gp)

	# write out PED file
	write.table(gp, file = filename, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

	# MAP file if appropriate
	if(!is.null(map)){
		write.table(data.frame(chrom = 1,
									  locus = loci,
									  position = 0),
						file = map, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
	}

}

#' write a Colony input file.
#'
#'
#' @param x an EFGLdata object
#' @param filename the name of the file to write. If NULL, the colony input is returned
#'   as a character vector with each line a separate entry in the vector.
#' @param pops a vector of pops to include. If not specified,
#'   all pops are used.
#' @param loci a vector of loci to include. If not specified,
#'   all loci are used.
#' @param candMales candidate male names
#' @param candFemales candidate female names
#' @param candMaleProb probability father is included in candidate males. if NULL, either 0 or 0.5 depending on whether candMales is NULL.
#' @param candFemaleProb probability father is included in candidate Females. if NULL, either
#'   0 or 0.5 depending on whether candFemales is NULL.
#' @param projName project name in colony input
#' @param colOutFileName output file name to direct colony to use
#' @param alleleDropoutRate allelic dropout rate to use for all loci
#' @paramotherGenotypingErrorRate other genotyping error rate to use for all loci
#' @return nothing, just writes a file
#' @export
exportColony <- function(x, filename = "Colony2.dat", pops = NULL, loci = NULL,
								 candMales = NULL, candFemales = NULL, candMaleProb = NULL,
								 candFemaleProb = NULL,
								 projName = "efglmhToColony",
								 colOutFileName = "colonyResults",
								 alleleDropoutRate = 0.01,
								 otherGenotypingErrorRate = 0.01){
	if(ncol(x$genotypes) < 3) stop("no genotypes")

	if(is.null(loci)) loci <- getLoci(x)
	if(is.null(pops)) pops <- getPops(x)
	if(any(!pops %in% getPops(x))) stop("not all pops are in this EFGLdata object")
	loci2 <- c(paste0(loci, ".A1"), paste0(loci, ".A2"))
	l <- colnames(x$genotypes)[3:ncol(x$genotypes)]
	if(any(!loci2 %in% l)) stop("one or more loci were not found in input")
	if(!is.null(candMales) && any(!candMales %in% x$genotypes$Ind)) stop("not all candMales found")
	if(!is.null(candFemales) && any(!candFemales %in% x$genotypes$Ind)) stop("not all candFemales found")
	if(is.null(candMaleProb)){
		if(is.null(candMales)){
			candMaleProb <- 0
		} else {
			candMaleProb <- 0.5
		}
	}
	if(is.null(candFemaleProb)){
		if(is.null(candFemales)){
			candFemaleProb <- 0
		} else {
			candFemaleProb <- 0.5
		}
	}


	g <- x$genotypes %>% filter(Pop %in% pops) %>% ungroup() # shouldn't be grouped, just making sure here
	gp <- g %>% select(Ind)

	# convert genotypes to numerical formatting
	to_remove <- c()
	alleleFreq <- list()
	for(i in loci){
		a <- g %>% select(paste0(i, ".A1"), paste0(i, ".A2"))
		alleleIndex <- a %>% tidyr::gather(locus, allele, 1:2) %>%
			filter(!is.na(allele)) %>% count(allele)
		if(nrow(alleleIndex) < 1){
			warning("all missing data for locus", i, ". Skipping this locus.")
			to_remove <- c(to_remove, i)
			next
		}
		alleleIndex <- alleleIndex %>% mutate(prop = n / sum(n), index = 1:length(allele))
		alleleFreq[[i]] <- alleleIndex

		a1 <- a %>% pull(1)
		a2 <- a %>% pull(2)
		a1 <- alleleIndex$index[match(a1, alleleIndex$allele)]
		a2 <- alleleIndex$index[match(a2, alleleIndex$allele)]
		a1[is.na(a1)] <- 0
		a2[is.na(a2)] <- 0
		gp <- gp %>% tibble::add_column(!!i := a1,
												  !!paste0(i, ".A2") := a2)
	}
	loci <- loci[!loci %in% to_remove]

	textOut <- c( # text vector with each entry a line
	paste0(projName, " ! Project name"),
	paste0(colOutFileName, " ! Output file name"),
	paste0(nrow(gp) - length(candMales) - length(candFemales), " ! Number of offspring in the sample"),
	paste0((ncol(gp) - 1) / 2, " !  Number of loci"),
	"7 ! Seed for random number generation",
	"1 ! 0/1 no/yes update allele frequencies",
	"2 ! 2/1 dioecious/monoecious",
	"0 ! 0/1 no/yes inbreeding",
	"0 ! 0/1 diploid/haplodiploid",
	"0 0 ! 0/1 polygamous/monogamous males females",
	"0 ! 0/1 no/yes infer clones",
	"1 ! 0/1 no/yes scale sibship",
	"0 ! 0/1/2/3/4 no/weak/medium/strong/calcualteFromNeAndR sibship prior",
	"1 ! 0/1 unknown/known allele frequencies",
	paste(sapply(alleleFreq, nrow), collapse = " ") # number of alleles for each locus
	)
	# Alleles and frequencies for each locus
	for(i in loci){
		textOut <- c(textOut,
		paste(alleleFreq[[i]]$index, collapse = " "),
		paste(alleleFreq[[i]]$prop, collapse = " "))
	}
	textOut <- c(
		textOut,
		"1 ! Number of replicated runs",
		"1 ! 1/2/3/4 short/medium/long/very long run",
		"0 ! 0/1 for monitor method",
		"10000 ! monitor interval",
		"0 ! 0/1 non-GUI/GUI",
		"1 ! 0/1/2 PLS/FL/FPLS likelihood calculation",
		"1 ! 0/1/2/3 Low/Medium/High/Very High precision in calculating the full likelihood",
		paste(loci, collapse = " "), # marker names
		paste(rep("0", length(loci)), collapse = " "), # marker co-dominant or dominant
		paste(rep(alleleDropoutRate, length(loci)), collapse = " "), # allelic dropout rate
		paste(rep(otherGenotypingErrorRate, length(loci)), collapse = " ") # other genotyping error rate
	)
	# split genotypes
	gp_cm <- gp %>% filter(Ind %in% candMales)
	gp_cf <- gp %>% filter(Ind %in% candFemales)
	gp <- gp %>% filter(!Ind %in% c(candMales, candFemales))

	# offspring genotypes
	for(i in 1:nrow(gp)) textOut <- c(textOut, paste(unlist(gp[i,]), collapse = " "))

	textOut <- c(
		textOut,
		paste0(candMaleProb, " ", candFemaleProb, " ! prob sire is in candMale prob dam is in candFemale"),
		paste0(length(candMales), " ", length(candFemales), " ! number of candMale and candFemale")
	)

	# candidate male genotypes
	if(nrow(gp_cm) > 0) for(i in 1:nrow(gp_cm)) textOut <- c(textOut, paste(unlist(gp_cm[i,]), collapse = " "))

	# candidate female genotypes
	if(nrow(gp_cf) > 0) for(i in 1:nrow(gp_cf)) textOut <- c(textOut, paste(unlist(gp_cf[i,]), collapse = " "))

	textOut <- c(
		textOut,
		"0 0 ! number offspring with known sires and mismatches allowed for known sires",
		# line(s) for known offspring-sire dyads ommited
		"0 0 ! number offspring with known dams and mismatches allowed for known dams",
		# line(s) for known offspring-dam dyads ommited
		"0 ! known paternal sibships",
		# known paternal sibship lines skipped
		"0 ! known maternal sibships",
		# known maternal sibship lines skipped
		"0 ! offspring with excluded sires",
		# excluded sires lines skipped
		"0 ! offspring with excluded dams",
		# excluded dam lines skipped
		"0 ! excluded paternal sibships",
		# excluded paternal sibship lines skipped
		"0 ! excluded maternal sibships"
		# excluded maternal sibship lines skipped
	)
	if(is.null(filename)){
		return(textOut)
	}
	writeLines(textOut, con = filename)
}
