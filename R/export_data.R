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
