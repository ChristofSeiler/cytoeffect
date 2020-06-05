#' Composite maximum likelihood estimation for Poisson log-normal model
#'
#' \code{fit_poilog} uses composite maximum likelihood estimation to fit the mean
#' and the covariance of a Poisson log-normal mixed model.
#'
#' @import PLNmodels
#' @import dplyr
#' @importFrom statmod gauss.quad
#' @importFrom parallel mclapply
#' @importFrom Matrix nearPD
#' @importFrom Matrix cov2cor
#' @importFrom matrixStats logSumExp
#' @importFrom magrittr %>% %<>%
#' @export
#'
#' @param Y Matrix with proteins counts with dimensions (number of cells) x (number of markers)
#' @param ncores Number of CPU cores to use
#' @return A list containing
#'   \item{beta}{estimated mean parameter}
#'   \item{Sigma}{estimated covariance parameter}
#'
#' @examples
#' set.seed(1)
#' df = simulate_data(n_cells = 10)
#' str(df)
#' Y = dplyr::select(df, names(df)[3:ncol(df)])
#' fit_poilog(Y)
#'
fit_poilog = function(Y, ncores = 1) {

  if(nrow(Y) < 2)
    stop("Y has zero or one observations")

  p = ncol(Y)

  # variational inference
  df = PLNmodels::prepare_data(counts = Y, covariates = rep(1, nrow(Y)))
  pln_fit = PLNmodels::PLN(Abundance ~ 1, df)
  beta_vi = coef(pln_fit) %>% as.vector()
  cov_vi = sigma(pln_fit) %>% as.matrix()
  sigma_vi = diag(cov_vi) %>% as.vector()
  cor_vi = cov2cor(cov_vi)

  # multivariate Gauss-Hermite quadrature

  # prepare grid and weights
  n = 1001
  prune_thresold = 0.00001
  gherm = statmod::gauss.quad(n, kind = "hermite")
  gherm %<>% as_tibble
  gherm %<>% add_column(index = 1:n)
  grid = tidyr::crossing(i = 1:n, j = 1:n)
  grid %<>% left_join(gherm, by = c("i" = "index"))
  grid %<>% left_join(gherm, by = c("j" = "index"))
  grid %<>% mutate(weight = weights.x*weights.y)
  grid %<>% filter(weight > prune_thresold)
  # # debugging plot
  # ggplot(grid, aes(nodes.x, nodes.y, size = weight)) +
  #   geom_point() +
  #   coord_equal()

  # converting pars for optim
  theta2mu = function(theta) {
    c(theta[1], theta[2])
  }
  theta2T = function(theta) {
    T =  matrix(0, nrow = 2, ncol = 2)
    T[1,1] = theta[3]
    T[2,2] = theta[4]
    T[2,1] = theta[5]
    T
  }

  # # log-likelihood over all the data
  # neg_log_lik_n = function(theta, positions, weights, Y) {
  #
  #   # parse parameters
  #   mu = theta2mu(theta)
  #   T = theta2T(theta)
  #
  #   # function from Aitchison and Ho, 1989
  #   logH = function(x, mu, T, u) {
  #     mu_Tu = mu + T %*% u * sqrt(2)
  #     sum(x*mu_Tu) - sum(exp(mu_Tu)) - sum(lgamma(x+1))
  #   }
  #
  #   # evaluate and weight H
  #   dpoislognorm = function(x) {
  #     log_fs = apply(positions, 1, function(u) logH(x, mu, T, u))
  #     # weights*fs = exp(log(weights*fs)) = exp(log(weights)+log(fs))
  #     logSumExp(log(weights)+log_fs)
  #   }
  #   log_lik = apply(Y, 1, function(x) dpoislognorm(x))
  #   -sum(log_lik)
  #
  # }

  # log-likelihood over all the data (vectorized version)
  neg_log_lik_n_vec = function(theta, positions, weights, Y) {

    # parse parameters
    mu = theta2mu(theta)
    T = theta2T(theta)

    # function from Aitchison and Ho, 1989 (vectorized)
    # sum(x*mu_Tu) - sum(exp(mu_Tu)) - sum(lgamma(x+1)) = A - B - C
    sqrt_2 = sqrt(2)
    a = mu[1] + T[1,1]*positions$nodes.x*sqrt_2
    b = mu[2] + T[2,1]*positions$nodes.x*sqrt_2 + T[2,2]*positions$nodes.y*sqrt_2
    B = exp(a) + exp(b)
    dpoislognorm_vec = function(x) {
      A = x[1]*a + x[2]*b
      C = lgamma(x[1]+1) + lgamma(x[2]+1)
      log_fs = A - B - C
      matrixStats::logSumExp(log(weights)+log_fs)
    }

    log_lik = apply(Y, 1, function(x) dpoislognorm_vec(x))
    -sum(log_lik)

  }

  # pairwise fitting using composite likelihood
  fit_pairwise = function(Y, i, j) {

    # make lower triangular matrix T
    mu = c(beta_vi[i], beta_vi[j])
    Sigma = matrix(nrow = 2, ncol = 2)
    Sigma[1,1] = diag(cov_vi)[i]
    Sigma[2,2] = diag(cov_vi)[j]
    Sigma[1,2] = Sigma[2,1] = cov_vi[i,j]
    T = t(chol(Sigma))

    # find maximum
    theta = c(mu[1], mu[2], T[1,1], T[2,2], T[2,1])
    fit = optim(par = theta, fn = neg_log_lik_n_vec,
                positions = grid %>% dplyr::select(nodes.x, nodes.y),
                weights = grid$weight,
                Y = Y[,c(i,j)],
                method = "BFGS",
                control = list(maxit = 1000, trace = 1, REPORT = 1))

    # convert for aggregation
    mu_mle = theta2mu(fit$par)
    T_mle = theta2T(fit$par)
    cov_mle = T_mle %*% t(T_mle)
    sigma_mle = sqrt(diag(cov_mle))
    rho_mle = cov2cor(cov_mle)[2,1]
    tibble(
      i = i, j = j,
      mu.i = mu_mle[1], mu.j = mu_mle[2],
      sig.i = sigma_mle[1], sig.j = sigma_mle[2],
      rho = rho_mle
    )

  }

  tb_pairs = tidyr::crossing(i = 1:p, j = 1:p) %>% dplyr::filter(i < j)
  tb_mle = mclapply(
    1:nrow(tb_pairs),
    function(row_id) fit_pairwise(Y, tb_pairs$i[row_id], tb_pairs$j[row_id]),
    mc.cores = ncores
  ) %>% bind_rows()
  avearge_mu = function(index) {
    a = tb_mle %>% dplyr::filter(i == index) %>% pull(mu.i)
    b = tb_mle %>% dplyr::filter(j == index) %>% pull(mu.j)
    mean(c(a, b))
  }
  beta_mle = sapply(1:p, avearge_mu)
  avearge_sig = function(index) {
    a = tb_mle %>% dplyr::filter(i == index) %>% pull(sig.i)
    b = tb_mle %>% dplyr::filter(j == index) %>% pull(sig.j)
    mean(c(a, b))
  }
  sigma_mle = sapply(1:p, avearge_sig)
  cor_mle = diag(nrow = p)
  for(row_id in 1:nrow(tb_mle)) {
    cor_mle[tb_mle$i[row_id], tb_mle$j[row_id]] = tb_mle$rho[row_id]
    cor_mle[tb_mle$j[row_id], tb_mle$i[row_id]] = tb_mle$rho[row_id]
  }
  cov_mle = diag(sigma_mle) %*% cor_mle %*% diag(sigma_mle)

  # make positive definite
  cov_mle = as.matrix(Matrix::nearPD(cov_mle, corr = FALSE)$mat)

  list(
    beta = beta_mle,
    Sigma = cov_mle
  )

}
