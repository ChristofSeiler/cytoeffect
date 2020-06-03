#' Generate Dataset for Vignettes and Simulation Studies
#'
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>% %<>%
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom Matrix toeplitz
#' @importFrom stringr str_pad
#' @export
#'
#' @param n_markers number of markers
#' @param n_true number of differential markers
#' @param n_samples donors per group
#' @param paired paired experimental design
#' @param n_cells per sample
#' @param pi zero inflation probability
#' @param beta_treatment log of mean in treatment
#' @param beta_control log of mean in control
#' @param rho_b cell level correlation
#' @param rho_u donor level correlation
#' @param seed set random seed
#' @return \code{\link[tibble]{tibble}} data frame
#'
simulate_data = function(n_markers = 5,            # number of markers
                         n_true = 1,               # number of differential markers
                         n_samples = 8,            # donors per group
                         paired = TRUE,            # paired experimental design
                         n_cells = 100,            # per sample
                         pi = 0.05,                # zero inflation probability
                         beta_treatment = log(20), # log of mean in treatment
                         beta_control = log(5),    # log of mean in control
                         rho_b = 0.5,              # cell level correlation
                         rho_u = 0.5,              # donor level correlation
                         seed = 0xdada             # set random seed
                         ) {

  # simulation parameters
  set.seed(seed)

  if (paired) {
    n_donors = n_samples / 2
  } else {
    n_donors = n_samples
  }

  # patient information
  donor = rep(1:n_donors, each = n_cells)
  if (paired) {
    donor %<>% rep(2)
    condition = c(rep("treatment", length(donor) / 2),
                  rep("control", length(donor) / 2))
  } else {
    condition = c(rep("treatment", length(donor) / 2),
                  rep("control", length(donor) / 2))
  }
  df = tibble(donor, condition)

  # generate protein counts
  protein_names = paste0("m", stringr::str_pad(1:n_markers, width = 2, pad = "0"))
  rcorr = function(rho) {
    corr = rho ^ toeplitz(0:(n_markers - 1))
    # random positive and negative correlations
    lower_tri = lower.tri(corr)
    upper_tri = upper.tri(corr)
    signs = sample(c(1, -1), size = sum(lower_tri), replace = TRUE)
    corr_sign = diag(1, n_markers)
    corr_sign[lower_tri] = signs * corr[lower_tri]
    corr_sign[upper_tri] = t(corr_sign)[upper_tri]
    # project to nearest positive definite matrix
    as.matrix(Matrix::nearPD(corr_sign, corr = TRUE)$mat)
  }
  Sigma_b = rcorr(rho_b) # cell level variability
  Sigma_u = rcorr(rho_u) # donor level variability
  b = MASS::mvrnorm(n = nrow(df), mu = rep(0, n_markers), Sigma_b)
  u = MASS::mvrnorm(n = n_donors, mu = rep(0, n_markers), Sigma_u)
  u = u[donor,]
  beta = matrix(beta_control, nrow = nrow(b), ncol = n_markers)
  # collider confounding
  beta[, 1:n_true] = ifelse(condition == "treatment", beta_treatment, beta_control)
  log_lambda = beta + b + u
  lambda = exp(log_lambda)
  y = rpois(length(lambda), lambda)
  dim(y) = dim(lambda)
  I = rbinom(n = length(y), size = 1, prob = pi)
  I = matrix(I, ncol = n_markers)
  y = y * (1 - I)
  colnames(y) = protein_names
  df %<>% bind_cols(as_tibble(y))
  df$condition %<>% factor(levels = c("control",
                                      "treatment"))
  df %<>% mutate(donor = paste0("pid", stringr::str_pad(donor, width = 2, pad = "0")))
  df

}
