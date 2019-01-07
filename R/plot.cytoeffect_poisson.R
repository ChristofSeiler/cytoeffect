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

  protein_names = obj$protein_names
  warmup = obj$fit_mcmc@stan_args[[1]]$warmup
  conditions = obj$conditions
  celltypes = obj$celltypes
  covariates = obj$covariates

  if(type == "beta") {

    tb_beta = summary(obj$fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as.tibble(rownames = "name")
    tb_beta %<>% add_column(protein_name = rep(protein_names, each = length(covariates)))
    tb_beta %<>% add_column(covariate = rep(covariates, length(protein_names)))
    #tb_beta %<>% dplyr::filter(covariate == covariates[2])
    ggplot(tb_beta, aes(x = `50%`, y = protein_name)) +
      geom_vline(xintercept = 0,color = "red") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.5) +
      ggtitle("Regression Coefficients") +
      xlab("log expected count") +
      theme(axis.title.y = element_blank()) +
      facet_wrap(~covariate)

  } else if (type == "sigma" || type == "sigma_term" || type == "sigma_donor") {

    tb_sigma = summary(obj$fit_mcmc, pars = type, probs = c(0.025, 0.5, 0.975))
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

    tb_sigma = summary(obj$fit_mcmc, pars = c("sigma","sigma_term","sigma_donor"),
                       probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% add_column(protein_name = rep(protein_names, 3))
    tb_sigma %<>% add_column(type = c(rep("cell", length(protein_names)),
                                      rep("term", length(protein_names)),
                                      rep("donor", length(protein_names))))
    tb_sigma$type %<>% factor(levels = c("cell","term","donor"))
    dodge = position_dodge(width = 0.9)
    ggplot(tb_sigma, aes(x = protein_name, y = `50%`, color = type)) +
      geom_point(size = 2,
                 position = dodge) +
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                    position = dodge) +
      ggtitle("Marker Standard Deviations") +
      ylab("sigma") +
      theme(axis.title.y = element_blank()) +
      coord_flip()

  } else if (type == "Cor") {

    lapply(c("Cor","Cor_term","Cor_donor"), function(type) {
      cor = rstan::extract(obj$fit_mcmc, pars = type)[[1]]
      cor_median = apply(X = cor, MARGIN = c(2,3), FUN = median)
      colnames(cor_median) = rownames(cor_median) = protein_names
      ggcorrplot(cor_median, hc.order = TRUE, type = "lower",
                 outline.col = "lightgray",
                 colors = c("#6D9EC1", "white", "#E46726")) +
        ggtitle(paste0("Marker Correlations (",type,")")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    })

  } else if (type == "Y_hat") {

    Y = obj$Y
    Y_hat = rstan::extract(obj$fit_mcmc, pars = type)[[1]]
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
