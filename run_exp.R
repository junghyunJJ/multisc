rm(list = ls())
library(tidyverse)
library(data.table)
library(gemma2) # https://github.com/fboehm/gemma2

library(jjutil) # devtools::install_github("junghyunJJ/jjuitl")
source("functions.R")

# load pheno data from gemma2 r package
pheno <- readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE)
Y <- as.matrix(pheno[, c(1, 6)])
Y %>% dim

# load pheno data from gemma2 r package
K <- readr::read_tsv(system.file("extdata", "mouse100.cXX.txt", package = "gemma2"), col_names = FALSE)[, 1:100]  %>% as.matrix

# load X data from gemma2 r package
geno <- readr::read_csv(system.file("extdata", "mouse100.geno.txt", package = "gemma2"), col_names = FALSE)

# We only focused on the 'rs8275764'
sel_snps <- "rs8275764"
sel_idx <- match(sel_snps, geno$X1)
g1 <- geno[c(sel_idx), -c(1:3)] # first 3 columns are SNP annotations!
X <- t(as.matrix(g1))
g1 %>% dim

# # test for multiple snp.
# NOTE!! the resutls of multiple snps and single snp is the same.
# sel_snps <- c("rs8275764", "rs6212654", "rs13477740")
# sel_idx <- match(sel_snps, geno$X1)
# g1 <- geno[c(sel_idx), -c(1:3)]
# X <- t(as.matrix(g1))

# 1. data transformation
dat <- transformation(Y, X, K)
dat %>% glimpse()

# 2. extimate vg and ve
res <- mvLMM(dat$Ystar, dat$Kva, dat$Kve)
res

# 3. data rotation
new_dat <- rotation(dat$Ystar, dat$Xstar, res$vg, res$ve, dat$Kva)
new_dat %>% glimpse()

# Ordinary least squares: (X'X)-1(X'y)
solve(t(new_dat$new_x) %*% new_dat$new_x) %*% t(new_dat$new_x) %*% new_dat$new_y

# 4. run single
res_single <- single_anlaysis(new_dat$new_y, new_dat$new_x)
res_single

# 4-1. run single using multiple linear
# NOTE!! the 4-1 is only for one SNP at this time
lm(new_dat$new_y ~ new_dat$new_x - 1) %>% summary

######################################################################
### multivariate analysis ############################################
######################################################################
# NOTE!! the multivariate analysis is only for one SNP at this time
# 5-1. run multivariate (gemma)
gemma <- fread("mouse100.assoc.txt") # results from gemma
gemma %>%
    filter(rs %in% sel_snps) %>%
    arrange(p_wald)

# 5-2. run reverse multivariate (Canonical Correlation Analysis; CCA)
vegan::CCorA(X = new_dat$new_y, Y = new_dat$new_x)

# 5-3.  run reverse multivariate (Permutational Multivariate Analysis of Variance Using Distance Matrices)
vegan::adonis2(new_dat$new_x ~ new_dat$new_y, method = "euclidean")

# 5-3.  run reverse multivariate (mCPC Aschard et al, AJHG, 2014)
mCPC(y = new_dat$new_x, g = new_dat$new_y)

# NOTE!! We need to check whether this approach is mathematically reasonable.
# 5-4. susieR
res_susie <- susieR::susie(y = new_dat$new_y, X = new_dat$new_x)
coef(res_susie)[-1]
res_susie$pip
res_susie$sets
