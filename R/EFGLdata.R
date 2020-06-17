# creating and methods for EFGL objects
#
#' read in data from a Progeny-style output file
#' or matrix/dataframe/tibble adn create an EFGLdata object
#'
#' @param input Either a character string to the tab-separated input file with
#'   a header row or a
#'   matrix/dataframe/tibble. Structure of the input: one row per individual.
#'   One column giving
#'   pedigree (population) names, one column giving individual names,
#'   optional additional metadata columns, genotype columns. Pedigree and
#'   individual name columns can be anywhere, if specified. Genotype columns
#'   MUST be consecutive and be the right most columns. Genotypes are given as
#'   two columns per call (diploidy assumed).
#' @param genotypeStart The column number that genotypes start at. If not specified,
#'   the first column with a column name ending in ".A1", ".a1", "-A1", or "-a1" is chosen.
#' @param pedigreeColumn The column number that contains pedigree (population) names.
#' @param nameColumn The column number that contains individual names. These MUST be unique.
#' @param convertNames TRUE to convert genotype and pedigree names in the same way
#'   that IDFGEN does (remove special
#'   characters from both and remove "." from genotype names).
#' @param convertMetaDataNames TRUE to remove special characters and spaces from metadata
#'   column names. This makes accessing them easier.
#' @param missingAlleles a vector of values (not NA) to treat as missing alleles. They will be converted to NA.
#' @param guess_max If input is a character string, this is the maximum number of lines to use when
#'   guessing input data types. Making this smaller results in quicker loading, making it larger
#'   can fix some parsing errors
#' @return An EFGLdata object, which is just a list with two elements. The first element
#' is a tibble with genotype data, the second is a tibble with metadata
#' @import dplyr
#' @export
#'
readInData <- function(input, genotypeStart = NULL, pedigreeColumn = 1, nameColumn = 2,
							  convertNames = TRUE, convertMetaDataNames = TRUE,
							  missingAlleles = c("0", "00", "000"), guess_max = 1e4){

	if(is.character(input)) input <- readr::read_tsv(input, guess_max = guess_max)
	if(!tibble::is_tibble(input)){
		warning("converting input to a tibble")
		input <- tibble::as_tibble(input)
	}

	# find first genotype column
	if(is.null(genotypeStart)) {
		genotypeStart <- which(grepl("[\\.-][Aa]1$", colnames(input)))
		if(length(genotypeStart) < 1) stop("no genotype columns found")
		genotypeStart <- genotypeStart[1]
	}

	# split into genotypes and metadata
	gD <- input %>% select(pedigreeColumn, nameColumn, genotypeStart:ncol(.))
	mD <- input %>% select(-(genotypeStart:ncol(.))) %>%
		select(pedigreeColumn, nameColumn, everything())
	rm(input)

	# convert genotype names
	mNames <- colnames(gD)[seq(3, ncol(gD) - 1, 2)]
	if(convertNames){
		mNames <- gsub("\\.", "", make.names(gsub("[\\.-][aA]1$", "", mNames)))
	}
	# standardize suffixes on allele 1 and 2
	colnames(gD)[seq(3, ncol(gD) - 1, 2)] <- paste0(mNames, ".A1")
	colnames(gD)[seq(4, ncol(gD), 2)] <- paste0(mNames, ".A2")

	# convert metadata names
	if(convertMetaDataNames) colnames(mD) <- make.names(colnames(mD))

	# standardize pedigree/population and individual name column names
	colnames(gD)[1:2] <- c("Pop", "Ind")
	colnames(mD)[1:2] <- c("Pop", "Ind")

	# convert pedigree names
	if(convertNames) {
		gD <- gD %>% mutate(Pop = make.names(Pop))
		mD <- mD %>% mutate(Pop = make.names(Pop))
	}

	# convert missing alleles into NA
	na_if_mult <- function(x){
		x[x %in% missingAlleles] <- NA
		return(x)
	}
	gD <- gD %>% mutate_at(vars(-Pop, -Ind), na_if_mult)

	return(construct_EFGLdata(list(genotypes = gD, metadata = mD)))
}



#' some basic checks on EFGLdata objects
#' @param x an EFGLdata object
#' @export
construct_EFGLdata <- function(x){

	# check for genotypes and metadata
	if(names(x)[1] != "genotypes") stop("First object is not genotypes in EFGLdata object")
	if(names(x)[2] != "metadata") stop("Second object is not metadata in EFGLdata object")

	# check pop and ind names
	if(colnames(x$genotypes)[1] != "Pop") stop("First column in genotypes must be Pop")
	if(colnames(x$genotypes)[2] != "Ind") stop("Second column in genotypes must be Ind")
	if(colnames(x$metadata)[1] != "Pop") stop("First column in metadata must be Pop")
	if(colnames(x$metadata)[2] != "Ind") stop("Second column in metadata must be Ind")

	# check for NAs in pop and ind
	if(any(is.na(x$genotypes$Pop))) stop("NA is not allowed for Pop")
	if(any(is.na(x$genotypes$Ind))) stop("NA is not allowed for Ind")

	# check that individual names are unique
	if(nrow(x$genotypes) != n_distinct(pull(x$genotypes, Ind))) stop("Individual names must be unique")

	# order by pop and ind
	x$genotypes <- x$genotypes %>% arrange(Pop, Ind)
	x$metadata <- x$metadata %>% arrange(Pop, Ind)

	# check that genotypes and metadata contain same pop and inds in order
	if(any(x$genotypes$Pop != x$metadata$Pop)) stop("Pop does not match between genotypes and metadata")
	if(any(x$genotypes$Ind != x$metadata$Ind)) stop("Ind does not match between genotypes and metadata")

	class(x) <- "EFGLdata"
	return(x)
}

#' print method for EFGLdata
#' @param x an EFGLdata object
#' @param ... ignored
#' @export
print.EFGLdata <- function(x, ...){
	cat("Populations:\n")
	print(sort(unique(x$genotypes$Pop)))
	cat("\n", (ncol(x$genotypes)/2) - 1, "loci\n")
}

#' An example input dataset used in the vignette
#' @format a tibble
"exampleData"
