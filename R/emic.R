# This R script calculates the effective median estimate of a gene expression
# causal effect (eta_em), the standard error of eta_em (sigma_em) and the
# p-value of nonzero causal effect as described in the paper, "Powerful and
# robust inference of causal genes of complex phenotypes with a novel median-based
# Mendelian randomization framework based on dependent expression quantitative
# trait loci". Variable names in this R script are consistent with the paper.


library(MASS)
set.seed(1234)

calculate_eta_em <- function(eta, eta_cor, w, n) {
  # reorder SNPs by their eta values
  order <- order(eta)
  eta <- eta[order]
  eta_cor <- eta_cor[order,][, order]
  w <- w[order]

  # calculate independency weights
  ne <- c()
  for (i in seq(n)) {
    eigenvals <- eigen(eta_cor[seq(i), seq(i)])$values
    eigenvals[eigenvals > 1] <- 1
    ne <- c(ne, sum(eigenvals))
  }
  ne1 <- c(0, ne[seq(n - 1)])
  delta <- ne - ne1

  # combine independency weights and squared inverse variance weights
  tau <- delta * w
  tau <- tau / sum(tau)

  # calculate the two effective middle orders, t and k
  s <- c() # sum of weights up to the weight of eta
  for (i in seq(n)) {
    s <- c(s, sum(tau[seq(i)]))
  }
  p <- c() # weight percentiles the eta
  for (i in seq(n)) {
    p <- c(p, (s[i] - tau[i] / 2))
  }
  t <- length(p[p < 0.5])
  k <- t + 1

  # calculate the effective median estimate of gene expression causal effect
  extrapolate <- (0.5 - p[t]) / (p[k] - p[t])
  eta_em <- eta[t] + (eta[k] - eta[t]) * extrapolate
  return(eta_em)
}


emic <- function(beta_y, sigma_y, beta_x, sigma_x, genotype_cor) {
  # betaY, betaY_se, betaX, and betaX_se are n-dimensional vector.
  # genotype_cor is a real symmetric n by n matrix.
  boot_n <- 1000  # the bootstrap sample size for sigma_em estimate
  c <- 1.1  # the empirical factor for the sigma_em adjustment

  # calculate eta_em
  eta <- beta_y / beta_x  # the vector of gene expression causal effect estimates
  eta_cor <- 0.4428 * (genotype_cor^2) + 0.5665 * abs(genotype_cor)
  diag(eta_cor) <- 1
  w <- (beta_x / sigma_y)^4  # squared inverse variance weights
  n <- length(beta_y)
  eta_em <- calculate_eta_em(eta, eta_cor, w, n)

  # calculate the se of eta_em
  beta_y_cov <- sigma_y %o% sigma_y * genotype_cor
  beta_x_cov <- sigma_x %o% sigma_x * genotype_cor
  beta_y_boot <- mvrnorm(boot_n, beta_y, beta_y_cov)
  beta_x_boot <- mvrnorm(boot_n, beta_x, beta_x_cov)
  eta_boot <- beta_y_boot / beta_x_boot
  eta_em_boot <- c()
  for (i in seq(boot_n)) {
    eta_em_boot <- c(eta_em_boot, calculate_eta_em(eta_boot[i,], eta_cor, w, n))
  }
  sigma_em <- sd(eta_em_boot)
  pval <- 2 * pnorm(-abs((eta_em / sigma_em * c)))
  return(c(n, eta_em, sigma_em, pval))
}


##### read data
gwas_beta <- read.table('testdata/gwas_beta.txt')$V1
gwas_beta_se <- read.table('testdata/gwas_beta_se.txt')$V1
eqtl_beta <- read.table('testdata/eqtl_beta.txt')$V1
eqtl_beta_se <- read.table('testdata/eqtl_beta_se.txt')$V1
ld_r <- data.matrix(read.table('testdata/ld_r.txt'))

##### run emic
result <- emic(gwas_beta, gwas_beta_se, eqtl_beta, eqtl_beta_se, ld_r)
label <- c('Number of IVs:', 'Effect size:', 'SE of effect size:', 'P-value:')
print(cbind(label, result))
