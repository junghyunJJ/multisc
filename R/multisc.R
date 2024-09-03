grm <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  p_hat <- apply(x, 2, sum) / (2 * n)
  w <- apply(rbind(x, p_hat), 2, function(z) {
    (z - 2 * z[length(z)]) / sqrt(2 * z[length(z)] * (1 - z[length(z)]))
  })[1:n, ]
  return(w %*% t(w) / m)
}
