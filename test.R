# testing functions
library(dplyr)
options(tibble.max_extra_cols = 10)

t <- readr::read_tsv("example_snp_mh.txt", guess_max = 1e4)
# data for vignettes
# exampleData <- t
# usethis::use_data(exampleData)
# usethis::use_vignette("How_to_use_EFGLmh")
which(grepl("[\\.-][Aa]1$", colnames(t)))
colnames(t) <- make.names(colnames(t))

t %>% select(2,1,-(3:10), everything())
t %>% select((1:3))

d <- readInData(t)
str(d)
d[[1]]
d[[2]]
d[[1]] %>% pull(Pop) %>% unique
combineEFGLdata(d,d,d)

d1 <- readInData(filter(t, Pedigree == "OmyOXBO19S"))
d2 <- readInData(filter(t, Pedigree != "OmyOXBO19S"))
d2$metadata <- d2$metadata %>% select(Pop, Ind, Gender)
d2$genotypes <- d2$genotypes %>% select(Pop, Ind)

combineEFGLdata(d1,d2, metaComb = "union")$metadata

exportProgenyStyle(d1, filename = "testProg.txt")

c <- combineEFGLdata(d1, d2, genoComb = "union", metaComb = "union")
c <- combineEFGLdata(d1, d2, genoComb = "intersect", metaComb = "intersect")
c[[1]]
tail(c[[1]])
c[[2]]
tail(c[[2]])

c <- movePops(readInData(t), c("OmyOXBO19S", "OmyLSCR19S"), "AAA")
numInds(c)
c[[1]]
c[[2]]

c <- moveInds(readInData(t), c("OmyOXBO19S_0001", "OmyLSCR19S_0005"), "OmyOXBO19S")

d1$genotypes %>% select(-Pop, -Ind) %>% tidyr::gather(locus, allele, 1:ncol(.)) %>%
	filter(!is.na(allele)) %>%
	mutate(locus = gsub("\\.A[12]$", "", locus)) %>% group_by(locus) %>%
	summarise(aRich = n_distinct(allele), .groups = "drop") %>%
	arrange(aRich) %>% as.data.frame()
aRich(d1)

a <- removeLoci(d1, lociRemove = unique(gsub("\\.A[12]$", "", colnames(d1$genotypes)[5:ncol(d1$genotypes)])))
a[[1]]
aRich(a)
removeInds(d1, getInds(d1))[[1]]

library(rubias)
d_temp <- removeLoci(d1, lociSuccess(d1) %>% filter(success < .05) %>% pull(locus))
# d_temp <- removeLoci(d_temp, getLoci(d_temp)[grepl("_mh$", getLoci(d_temp))])
rm <- exportRubias_mixture(d_temp)
rb <- exportRubias_baseline(d_temp, repunit = "Pop", collection = "Pop")
rb$indiv <- paste0(rb$indiv, "_2")
rm$collection <- "one"
dups <- close_matching_samples(D = bind_rows(rm, rb), gen_start_col = 5)
dups <- close_matching_samples(D = rb, gen_start_col = 5)
# note: rubias throws error when only mixture samples, but not when only
#  baseline samples for checking for duplicates
dups


exportGenoPop(d, "test.txt", header = "genePop file testing testing",
				  pops = "OmyOXBO19S", loci = NULL)
genepop::basic_info("test.txt", outputFile = "testout.txt", verbose = interactive())

exportGrandma(d, baseline = TRUE)
exportGrandma(d, baseline = FALSE)

exportGenAlEx(d, filename = "testGenAlEx.txt", title = "test title here")
exportGenePop(d, filename = "testGenePop.txt")

genoSuccess(d, loci = getLoci(d))

d_fil <- genoSuccess(d) %>% filter(success < .9) %>% pull(Ind) %>% removeInds(x = d, inds = .)

numInds(d)
numInds(d_fil)

exportSNPPIT(d_fil, filename = "testSNPPIT.txt", baseline = c("OmyDWOR19S", "OmyEFSW19S"), mixture = c("OmyLYON19S", "OmyOXBO19S"))

dupTable
