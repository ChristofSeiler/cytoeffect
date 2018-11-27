#' Plot posterior summaries.
#'
#' @import rstan
#' @import ggplot2
#' @import ggcorrplot
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @export
#'
plot.cytoeffect_poisson = function(obj, type = "distribution") {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  fit_mcmc = obj$fit_mcmc
  protein_names = obj$protein_names
  warmup = fit_mcmc@stan_args[[1]]$warmup
  conditions = obj$conditions
  celltypes = obj$celltypes
  covariates = obj$covariates
  Y = obj$Y

  if(type == "beta") {

    tb_beta = summary(fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as.tibble(rownames = "name")
    tb_beta %<>% add_column(protein_name = rep(protein_names, each = length(covariates)))
    tb_beta %<>% add_column(covariate = rep(covariates, length(protein_names)))
    tb_beta %<>% dplyr::filter(covariate == covariates[2])
    ggplot(tb_beta, aes(x = `50%`, y = protein_name)) +
      geom_vline(xintercept = 0,color = "red") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
      ggtitle("Regression Coefficients") +
      xlab(covariates[2]) +
      theme(axis.title.y = element_blank())

  } else if (type == "sigma" || type == "sigma_donor") {

    tb_sigma = summary(fit_mcmc, pars = type, probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% add_column(protein_name = protein_names)
    ggplot(tb_sigma, aes(x = `50%`, y = protein_name)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
      ggtitle("Marker Standard Deviations") +
      xlab(type) +
      theme(axis.title.y = element_blank())

  } else if (type == "sigma_all") {

    tb_sigma = summary(fit_mcmc, pars = c("sigma","sigma_donor"),
                       probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% add_column(protein_name = rep(protein_names, 2))
    tb_sigma %<>% add_column(type = c(rep("cell", length(protein_names)),
                                      rep("donor", length(protein_names))))
    ggplot(tb_sigma, aes(x = `50%`, y = protein_name, color = type)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
      ggtitle("Marker Standard Deviations") +
      xlab("sigma") +
      theme(axis.title.y = element_blank())

  } else if (type == "Cor" || type == "Cor_donor") {

    cor = rstan::extract(fit_mcmc, pars = type)[[1]]
    cor_median = apply(X = cor, MARGIN = c(2,3), FUN = median)
    cor_pos = apply(X = cor, MARGIN = c(2,3), FUN = function(x) mean(x > 0))
    cor_neg = apply(X = cor, MARGIN = c(2,3), FUN = function(x) mean(x < 0))
    p_mat = 1-pmax(cor_pos, cor_neg)
    colnames(cor_median) = rownames(cor_median) = protein_names
    ggcorrplot(cor_median, hc.order = TRUE, type = "lower",
               outline.col = "lightgray",
               colors = c("#6D9EC1", "white", "#E46726")) +
      ggtitle(paste0("Marker Correlations (",type,")")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

  } else if (type == "Y_hat") {

    Y_hat = rstan::extract(fit_mcmc, pars = type)[[1]]
    Y_hat_mean = apply(X = Y_hat, MARGIN = c(2,3), FUN = mean)
    # Standard residuals
    #Y_residual = Y - Y_hat_mean
    # Pearson residuals
    #Y_residual = (Y - Y_hat_mean)/sd(Y_hat_mean)
    # Deviance residuals
    Y_residual = sign(Y - Y_hat_mean)*sqrt(2*(Y*log(Y/Y_hat_mean) - (Y-Y_hat_mean)))
    Y_residual %<>% as.tibble
    ggplot(Y_residual %>% gather(marker, expr, protein_names), aes(expr)) +
      geom_histogram(bins = 30) +
      facet_wrap(~marker) +
      ggtitle("Deviance Residuals")
    # plot residuals correlations
    # Y_residual = Y - Y_hat_mean
    # ggcorrplot(cor(Y_residual), hc.order = TRUE, type = "lower",
    #            outline.col = "lightgray",
    #            colors = c("#6D9EC1", "white", "#E46726")) +
    #   ggtitle(paste0("Marker Correlations (",type,")")) +
    #   theme(panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank())

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
