rm(list = ls())
library(tidyverse)
library(data.table)
library(gemma2) # https://github.com/fboehm/gemma2

library(jjutil) # devtools::install_github("junghyunJJ/jjuitl")
source("R/multisc.R")

# load pheno data from gemma2 r package
pheno <- readr::read_tsv(system.file("extdata", "mouse100.pheno.txt", package = "gemma2"), col_names = FALSE)
Y <- as.matrix(pheno[, c(1, 6)])
Y %>% dim
sum(diag(cov(Y)))
# 1.666482

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

######################################################################
### 1.GEMMA ##########################################################
######################################################################

# NOTE!! the multivariate analysis is only for one SNP at this time
# 5-1. run multivariate (gemma)
# /common/jungj2/miniconda3/bin/gemma \
#     -g /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt \
#     -p /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt \
#     -k /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.cXX.txt \
#     -lmm \
#     -n 1 6 \
#     -o mouse100

raw_res_gemma <- fread("data/mouse100.assoc.txt") # results from gemma
res_gemma <- raw_res_gemma %>%
  filter(rs %in% sel_snps) %>%
  arrange(p_wald)
res_gemma$p_wald
# [1] 3.282434e-06

######################################################################
### 2. Canonical Correlation Analysis (CCA) ##########################
######################################################################

res_CCA <- vegan::CCorA(X = X, Y = Y)
res_CCA
# Canonical Correlation Analysis

# Call:
# vegan::CCorA(Y = Y, X = X) 

#              Y X
# Matrix Ranks 2 1

# Pillai's trace:  0.2039568 

# Significance of Pillai's trace:
# from F-distribution:   1.5683e-05 
#                        CanAxis1
# Canonical Correlations   0.4516

#                       Y | X  X | Y
# RDA R squares      0.102155 0.2040
# adj. RDA R squares 0.092993 0.1875

###############################################################################
### 3. Permutational Multivariate Analysis of Variance Using Distance Matrices#
###############################################################################

res_adonis2 <- vegan::adonis2(Y ~ X,  method = "euclidean")
res_adonis2
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# vegan::adonis2(formula = Y ~ X, method = "euclidean")
#          Df SumOfSqs      R2     F Pr(>F)
# X         1   16.854 0.10216 11.15  0.001 ***
# Residual 98  148.128 0.89784
# Total    99  164.982 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

######################################################################
### 4. mCPC (Aschard et al, AJHG, 2014) ##############################
######################################################################
# (mCPC Aschard et al, AJHG, 2014)
res_mCPC <- mCPC(Y, X)
res_mCPC
# [1] 3.195272e-05

######################################################################
### 5.multisc ########################################################
######################################################################
# (mCPC Aschard et al, AJHG, 2014)
res_multisc <- multisc(Y, X, K)
res_multisc
# [1] 2.879274e-05


######################################################################
### 6. hsq ###########################################################
### 6-1. ori hsq #####################################################

Y_std <- scale(Y)
sum(diag(cor(Y_std)))

# e1 <- lmmlite::eigen_rotation(K, Y[, 1], cbind(1, X))
# res_lmm1 <- lmmlite::fitLMM(e1$Kva, e1$y, e1$X)
# vg_lmm1 <- res_lmm1$sigmasq_g
# ve_lmm1 <- res_lmm1$sigmasq_e
# res_lmm1$hsq
# # 0.8951924

# e2 <- lmmlite::eigen_rotation(K, Y[, 2], cbind(1, X))
# res_lmm2 <- lmmlite::fitLMM(e2$Kva, e2$y, e2$X)
# vg_lmm2 <- res_lmm1$sigmasq_g
# ve_lmm2 <- res_lmm1$sigmasq_e
# res_lmm2$hsq
# # 0.468617

# (vg_lmm1 + vg_lmm2) / (vg_lmm1 + ve_lmm1 + vg_lmm2 + ve_lmm2)
# # 0.8951924

# (vg_lmm1 + vg_lmm2)
# # 4.278105


# res_lmm1 <- qgg::greml(y = Y[, 1], X = cbind(1, X), GRM = list(K))
res_lmm1 <- qgg::greml(y = Y_std[, 1], X = cbind(1, X), GRM = list(K))
vg_lmm1 <- res_lmm1$theta[1]
ve_lmm1 <- res_lmm1$theta[2]
(hsq1 <- vg_lmm1 / (vg_lmm1 + ve_lmm1))
# 0.895193

# res_lmm2 <- qgg::greml(y = Y[, 2], X = cbind(1, X), GRM = list(K))
res_lmm2 <- qgg::greml(y = Y_std[, 2], X = cbind(1, X), GRM = list(K))
vg_lmm2 <- res_lmm2$theta[1]
ve_lmm2 <- res_lmm2$theta[2]
(hsq2 <- vg_lmm2 / (vg_lmm2 + ve_lmm2))
# 0.4686314

(vg_lmm1 + vg_lmm2) / (vg_lmm1 + ve_lmm1 + vg_lmm2 + ve_lmm2)
# 0.7634734

vg_lmm1 + vg_lmm2
# 2.86141

hsq1 + hsq2
# 1.363824

######################################################################
### 6. hsq ###########################################################
### 6-2. PCA hsq #####################################################

calpc <- function(y) {
  m <- ncol(y)
  eig <- eigen(cor(y))
  # eig <- eigen(cov(y))

  # calculate the eigenvalue and eigenvector of the correlation(or covariance) matrix of y
  eig.val <- eig$values
  eig.vec <- eig$vectors
  PC <- y %*% eig.vec
  return(list(PC = PC, eigval = eig.val))
}

res_PCA <- calpc(Y)
PC <- res_PCA$PC
# sum(diag(cov(PC)))
# # 1.666482

sum(diag(cor(PC)))
# 2

sum(res_PCA$eigval)
# 2

sum(diag(cov(Y)))
# 1.666482

# e_PC1 <- lmmlite::eigen_rotation(K, PC[, 1], cbind(1, X))
# res_PC1_lmm <- lmmlite::fitLMM(e_PC1$Kva, e_PC1$y, e_PC1$X)
 
# vg_PC1_lmm <- res_PC1_lmm$sigmasq_g
# ve_PC1_lmm <- res_PC1_lmm$sigmasq_e
# (hsq_PC1 <- res_PC1_lmm$hsq)
# # 0.71283

# e_PC2 <- lmmlite::eigen_rotation(K, PC[, 2], cbind(1, X))
# res_PC2_lmm <- lmmlite::fitLMM(e_PC2$Kva, e_PC2$y, e_PC2$X)
 
# vg_PC2_lmm <- res_PC2_lmm$sigmasq_g
# ve_PC2_lmm <- res_PC2_lmm$sigmasq_e
# (hsq_PC2 <- res_PC2_lmm$hsq)
# # 0.7108151


# ((hsq_PC1 * res_PCA$eigval[1]) + (hsq_PC2 * res_PCA$eigval[2])) / sum(res_PCA$eigval)
# # 0.7120807

# (vg_PC1_lmm + vg_PC2_lmm) / (vg_PC1_lmm + ve_PC1_lmm + vg_PC2_lmm + ve_PC2_lmm)
# # 0.7120687

# vg_PC1_lmm + vg_PC2_lmm
# # 2.053598


res_PC1_lmm <- qgg::greml(y = PC[, 1], X = cbind(1, X), GRM = list(K))
vg_PC1_lmm <- res_PC1_lmm$theta[1]
ve_PC1_lmm <- res_PC1_lmm$theta[2]
(hsq_PC1 <- vg_PC1_lmm / (vg_PC1_lmm + ve_PC1_lmm))
# 0.6441243

res_PC2_lmm <- qgg::greml(y = PC[, 2], X = cbind(1, X), GRM = list(K))
vg_PC2_lmm <- res_PC2_lmm$theta[1]
ve_PC2_lmm <- res_PC2_lmm$theta[2]
(hsq_PC2 <- vg_PC2_lmm / (vg_PC2_lmm + ve_PC2_lmm))
# 0.8865233

((hsq_PC1 * res_PCA$eigval[1]) + (hsq_PC2 * res_PCA$eigval[2])) / sum(res_PCA$eigval)
# 0.7369442

# (vg_PC1_lmm + vg_PC2_lmm) / (vg_PC1_lmm + ve_PC1_lmm + vg_PC2_lmm + ve_PC2_lmm)

vg_PC1_lmm + vg_PC2_lmm
# 2.444472

hsq_PC1 + hsq_PC2
# 1.530648