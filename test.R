# testing functions
library(dplyr)
options(tibble.max_extra_cols = 10)

t <- readr::read_tsv("example_snp_mh.txt", guess_max = 1e4)
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

c <- moveInds(readInData(t), c("OmyOXBO19S_0001", "OmyLSCR19S_0005"), "AAA")
