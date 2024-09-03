###############################################################################
#
#     ######      #      #######  #######
#    #           # #        #     #
#    #  ######  #    #      #     #####
#    #    #    # # # #      #     #
#    ######   #       #     #     #######
#
#   This script runs the GATE procedure, described in
# Wei Zhang, Liu Yang, Larry L Tang, Aiyi Liu, James L Mills, Yuanchang Sun, and Qizhai Li
# Hadling  with multiple correlated phenotypes association studies:
###################################################################################

################### Explanation of variables    ###################################
##### Input
## y: the phenotype matrix where each column represents a phenotype and each row represents an individual
## g: the genotype vector for a SNP
## x: the covariate matrix which needs to be adjusted for (the dimension is n by q, where n is the sample size and q is the number of covariate)
## B: the number of Monte Carlo replications usd to calculate the empirical distirbution of the test statistics
##### Output
## gate_pval: is the p value of the gate statistic
####################################################################################

multisc <- function(y, g, K, x = NULL) {
  ## calculate the principal component of the multiple trait and the chi-squared statistics using lm
  m <- ncol(y)
  eig <- eigen(cor(y))
  ### calculate the eigenvalue and eigenvector of the correlation matrix of y
  eig.val <- eig$values
  eig.vec <- eig$vectors
  PC <- y %*% eig.vec
  temp.Tstat <- rep(NA, m)
  for (j in 1:m) {
    if (is.null(x)) {
      # browser()
      # temp.Tstat[j] <- summary(lm(PC[, j] ~ g))$coef[2, 3]
      res <- qgg::greml(y = PC[, j], X = cbind(1, X), GRM = list(K))
      # browser()
      # NOTE!!!! we need to check the res$b and res$vb for T statistics
      temp.Tstat[j] <- res$b[2] / sqrt(res$vb[2, 2])
    } else {
      # temp.Tstat[j] <- summary(lm(PC[, j] ~ g + x))$coef[2, 3]
      res <- qgg::greml(y = PC[, j], X = cbind(1, X) + x, GRM = list(K))
      temp.Tstat[j] <- res$b[2] / sqrt(res$vb[2, 2])
    }
  }
  chi.PC <- temp.Tstat^2
  ### the test statistics based on all PC (put them in a group)
  T.QF <- sum(chi.PC)
  pv.QF <- 1 - pchisq(T.QF, df = m)

  ### the test statistic for the grouping with 2 groups
  com1.stat <- cumsum(chi.PC[1:(m - 1)])
  com2.stat <- rev(cumsum((rev(chi.PC))[1:(m - 1)]))
  pv1 <- 1 - pchisq(com1.stat, df = 1:(m - 1))
  pv2 <- 1 - pchisq(com2.stat, df = (m - 1):1)
  T.2g <- -2 * log(pv1) - 2 * log(pv2)
  pv.2g <- 1 - pchisq(T.2g, df = 4)
  return(min(c(pv.QF, pv.2g)))
  # ## the minimal p value of all possible groupings
  # gate_stat <- min(c(pv.QF, pv.2g))
  # reg_gate <- refGate(m, B=B)
  # gate_pval <- mean(reg_gate<gate_stat)
  # return(gate_pval)
}

mCPC <- function(y, g, x=NULL)
{
  ##calculate the principal component of the multiple trait and the chi-squared statistics using lm
  m <- ncol(y)
  eig <- eigen(cor(y))
  ### calculate the eigenvalue and eigenvector of the correlation matrix of y
  eig.val <- eig$values
  eig.vec <- eig$vectors
  PC <- y %*% eig.vec
  temp.Tstat <- rep(NA, m)
  for(j in 1:m){
    if(is.null(x)) {
      temp.Tstat[j] <- summary(lm(PC[,j]~g))$coef[2,3]
    } else {
      temp.Tstat[j] <- summary(lm(PC[,j]~g+x))$coef[2,3]
    }
  }
  chi.PC <- temp.Tstat^2
  ### the test statistics based on all PC (put them in a group)
  T.QF <- sum(chi.PC)
  pv.QF <- 1-pchisq(T.QF, df=m)

  ### the test statistic for the grouping with 2 groups
  com1.stat <- cumsum(chi.PC[1:(m-1)])
  com2.stat <- rev(cumsum((rev(chi.PC))[1:(m-1)]))
  pv1 <- 1-pchisq(com1.stat, df=1:(m-1))
  pv2 <- 1-pchisq(com2.stat, df=(m-1):1)
  T.2g <- -2*log(pv1)-2*log(pv2)
  pv.2g <- 1-pchisq(T.2g, df=4)
  return(min(c(pv.QF, pv.2g)))
  # ## the minimal p value of all possible groupings
  # gate_stat <- min(c(pv.QF, pv.2g))
  # reg_gate <- refGate(m, B=B)
  # gate_pval <- mean(reg_gate<gate_stat)
  # return(gate_pval)
}

GATE_fun <- function(y, g, x, B)
{
    ##calculate the principal component of the multiple trait and the chi-squared statistics using lm
    m <- ncol(y)
    eig <- eigen(cor(y))
    ### calculate the eigenvalue and eigenvector of the correlation matrix of y
    eig.val <- eig$values
    eig.vec <- eig$vectors
    PC <- y %*% eig.vec
    temp.Tstat <- rep(NA, m)
    for(j in 1:m)
    {
        if(is.null(x))  temp.Tstat[j] <- summary(lm(PC[,j]~g))$coef[2,3]
        else
        {
            temp.Tstat[j] <- summary(lm(PC[,j]~g+x))$coef[2,3]
        }
     }
     chi.PC <- temp.Tstat^2
     ### the test statistics based on all PC (put them in a group)
     T.QF <- sum(chi.PC)
     pv.QF <- 1-pchisq(T.QF, df=m)
     
     ### the test statistic for the grouping with 2 groups
     com1.stat <- cumsum(chi.PC[1:(m-1)])
     com2.stat <- rev(cumsum((rev(chi.PC))[1:(m-1)]))
     pv1 <- 1-pchisq(com1.stat, df=1:(m-1))
     pv2 <- 1-pchisq(com2.stat, df=(m-1):1)
     T.2g <- -2*log(pv1)-2*log(pv2)
     pv.2g <- 1-pchisq(T.2g, df=4)
     ## the minimal p value of all possible groupings
     gate_stat <- min(c(pv.QF, pv.2g))
     reg_gate <- refGate(m, B=B)
     gate_pval <- mean(reg_gate<gate_stat)
     return(gate_pval)
}

#############################################################################################
#### the refernce cumulative distirbution function of the gate statistic for m correlated phenotypes
refGate <- function(m, B) 
{
  samp <- matrix(rchisq(m*B, df=1, ncp=0), nrow=B, ncol=m)
  T.QF <- rowSums(samp)
  pv.QF <- 1 - pchisq(T.QF, df=m)
  com1.stat <- matrix(apply(samp[, -m, drop=F], 1, cumsum), nrow=m-1, ncol=B)
  com2.stat <- matrix(apply(samp[, -1, drop=F], 1, function(u) rev(cumsum(rev(u)))), nrow=m-1, ncol=B)
  pv1 <- apply(com1.stat, 2, function(u) 1-pchisq(u, df=1:(m-1)))
  pv2 <- apply(com2.stat, 2, function(u) 1-pchisq(u, df=(m-1):1))
  T.2g <- - 2*log(pv1) - 2*log(pv2)
  pv.2g <- 1 - pchisq(T.2g, df=4)
  ref.GateStat <- apply(rbind(pv.QF, pv.2g), 2, min)
  return(ref.GateStat)
}


