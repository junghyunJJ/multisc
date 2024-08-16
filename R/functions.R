grm <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  p_hat <- apply(x, 2, sum) / (2 * n)
  w <- apply(rbind(x, p_hat), 2, function(z) {
    (z - 2 * z[length(z)]) / sqrt(2 * z[length(z)] * (1 - z[length(z)]))
  })[1:n, ]
  return(w %*% t(w) / m)
}

# the function is exactly same as gemma2::eigen2()
eigh <- function(spd, decreasing = FALSE) {
  foo <- eigen(spd, symmetric = TRUE)
  bar <- foo
  bar$values <- foo$values[order(foo$values, decreasing = decreasing)]
  bar$vectors <- foo$vectors[, order(foo$values, decreasing = decreasing)]
  return(bar)
}

transformation <- function(Y, X, K) {
  # eigendecomposition of K
  K_eigen <- eigh(K)
  Kva <- K_eigen$values
  Kve <- K_eigen$vectors
  leftTransform <- Kve %>% t()

  # stadadization of Y
  Y <- scale(Y)
  Ystar <- leftTransform %*% Y

  # X centering and transform
  X <- scale(X, center = TRUE, scale = FALSE)
  # We dont need to add the interaction term because we perform centeing the X matrix
  # X0_T <- leftTransform %*% matrix(rep(1, n_indi))
  Xstar <- leftTransform %*% X
  return(list(Ystar = Ystar, Xstar = Xstar, Kva = Kva, Kve = Kve))
}

# NOTE!!!!!! we need to update this function using mgreml (https://github.com/devlaming/mgreml?tab=readme-ov-file#variance-components)
mvLMM <- function(Ystar, Kva, Kve) {
  n_indi <- nrow(Ystar)
  n_pheno <- ncol(Ystar)

  # RUN gemma2 using EM method
  # GEMMA (https://github.com/genetics-statistics/GEMMA) and gemma2 (https://cran.r-project.org/web/packages/gemma2/vignettes/compare-with-GEMMAv096.html)
  res <- gemma2::MphEM(
    eval = Kva,
    X = t(matrix(rep(1, n_indi))) %*% Kve,
    Y = t(Ystar),
    V_g = diag(n_pheno),
    V_e = diag(n_pheno)
  )[[1]]
  Psi <- res$Vg
  Phi <- res$Ve
  return(list(vg = Psi, ve = Phi))
}

# rotation function using chol (Joo et al, Genetics 2016) | GAMMA (https://academic.oup.com/genetics/article/204/4/1379/6046819?login=true)
chol_solve <- function(K) {
  a <- eigen(K)$vectors
  b <- eigen(K)$values
  b[b < 1e-13] <- 1e-13
  b <- 1 / sqrt(b)
  return(a %*% diag(b) %*% t(a))
}

rotate <- function(Y, sigma) {
  U <- chol_solve(sigma)
  tU <- t(U)
  UY <- tU %*% Y
  return(UY)
}

rotation <- function(Ystar, Xstar, Psi, Phi, Kva) {

  Psi_eigen <- eigh(Psi)
  Psi_Kva <- Psi_eigen$values
  Psi_Kve <- Psi_eigen$vectors

  Phi_eigen <- eigh(Phi)
  Phi_Kva <- Phi_eigen$values
  Phi_Kve <- Phi_eigen$vectors

  Psi_Kva[Psi_Kva == 0] <- 1e-13
  Phi_Kva[Phi_Kva == 0] <- 1e-13

  # Diagonalizing two matrices (Furlotte et al, Genetics 2015) | mvLMM (https://academic.oup.com/genetics/article/200/1/59/5936216)
  R <- Psi_Kve %*% sqrt(solve(diag(Psi_Kva))) # Now R %*% t(R) is equal to Psi^{-1}
  RR <- (t(R) %*% Phi) %*% R
  eigen_values <- eigh(RR)$values
  D <- eigen_values
  rightTransform <- t(eigh(RR)$vectors) %*% t(R)

  # P: the variance of the jth phenotype
  # P <- c(Kva + D[1], Kva + D[2], Kva + D[3])
  P <- lapply(seq_len(length(D)), function(d) {
    Kva + D[d]
  }) %>% unlist

  # The transformed vector of Y
  Yt <- as.vector(Ystar %*% t(rightTransform)) %>% as.matrix()
  #Yt <- as.vector(Ystar %*% rightTransform) %>% as.matrix()

  # The transformed vector of X
  Xt <- kronecker(rightTransform, Xstar)

  # data rotation 
  new_x <- rotate(Xt, diag(P))
  new_y <- rotate(Yt, diag(P))    

  return(list(new_y = new_y, new_x = new_x))
}

single_anlaysis <- function(new_y, new_x) {
  n_snp <- ncol(new_x) / 2

  res_lm <- lapply(seq_len(n_snp), function(xx) {
    tmpres <- summary(lm(new_y ~ new_x[, c(xx, (xx + n_snp))] - 1))$coefficients
    n_pheno <- nrow(tmpres)

    tmpres2 <- lapply(seq_len(n_pheno), function(pp) {
      sel_tmpres <- tmpres[pp, ]
      sel_tmpres <- t(data.frame(sel_tmpres))
      colnames(sel_tmpres) <- paste0(c("beta", "betastderr", "t", "p"), "_", pp)
      sel_tmpres %>% as.data.table()
    })
    do.call(cbind.data.frame, tmpres2)
  }) %>% rbindlist()

  return(res_lm)
}

# Aschard et al, AJHG, 20214 | mCPC (https://www.cell.com/ajhg/fulltext/S0002-9297(14)00118-9)
# please also check the following methods:
#   GATE (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3928-7)
#   APCT (https://link.springer.com/article/10.1007/s10255-022-1019-2)
mCPC <- function(y, g, x = NULL) {

  # calculate the principal component of the multiple trait and the chi-squared statistics using lm
  m <- ncol(y)
  eig <- eigen(cor(y))

  # calculate the eigenvalue and eigenvector of the correlation matrix of y
  eig.val <- eig$values
  eig.vec <- eig$vectors
  PC <- y %*% eig.vec
  temp.Tstat <- rep(NA, m)
  for (j in 1:m) {
    if (is.null(x)) {
      temp.Tstat[j] <- summary(lm(PC[, j] ~ g))$coef[2, 3]
    } else {
      temp.Tstat[j] <- summary(lm(PC[, j] ~ g + x))$coef[2, 3]
    }
  }

  chi.PC <- temp.Tstat^2
  # the test statistics based on all PC (put them in a group)
  T.QF <- sum(chi.PC)
  pv.QF <- 1 - pchisq(T.QF, df = m)

  # the test statistic for the grouping with 2 groups
  com1.stat <- cumsum(chi.PC[1:(m - 1)])
  com2.stat <- rev(cumsum((rev(chi.PC))[1:(m - 1)]))
  pv1 <- 1 - pchisq(com1.stat, df = 1:(m - 1))
  pv2 <- 1 - pchisq(com2.stat, df = (m - 1):1)
  T.2g <- -2 * log(pv1) - 2 * log(pv2)
  pv.2g <- 1 - pchisq(T.2g, df = 4)
  return(min(c(pv.QF, pv.2g)))
}
