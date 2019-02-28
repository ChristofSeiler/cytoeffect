#' Plot posterior summaries.
#'
#' @import rstan
#' @import ggplot2
#' @import ggcorrplot
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @export
#'
plot.cytoeffect_poisson = function(obj, type = "distribution") {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  warmup = obj$fit_mcmc@stan_args[[1]]$warmup
  protein_names = obj$protein_names
  conditions = obj$conditions

  if(type == "beta") {

    # extract index estimates
    tb_beta = summary(obj$fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as.tibble(rownames = "name")
    tb_beta %<>% add_column(protein_name = rep(protein_names, each = length(conditions)))
    tb_beta %<>% add_column(condition = rep(conditions, length(protein_names)))

    # add contrast
    beta = rstan::extract(obj$fit_mcmc, pars = "beta")[[1]]
    beta_contrast = beta[,,2] - beta[,,1]
    beta_diff = apply(beta_contrast, 2,
                      function(x)
                        quantile(x, probs = c(0.025, 0.5, 0.975))
    ) %>% t %>% as.tibble
    beta_diff %<>% add_column(name = "contrast")
    beta_diff %<>% add_column(protein_name = protein_names)
    beta_diff %<>% add_column(condition = paste(rev(conditions),collapse = " - "))
    tb_beta %<>% bind_rows(beta_diff)

    # combined plot
    dodge = position_dodge(width = 0.9)
    ggplot(tb_beta, aes(x = protein_name, y = `50%`, color = condition)) +
      geom_hline(yintercept = 0,color = "grey") +
      geom_point(size = 2, position = dodge) +
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                    position = dodge) +
      ggtitle("Regression Coefficients") +
      ylab("log expected count") +
      theme(axis.title.y = element_blank()) +
      coord_flip()

  } else if (type == "sigma") {

    tb_sigma = summary(obj$fit_mcmc, pars = c("sigma","sigma_term","sigma_donor"),
                       probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% add_column(protein_name = rep(protein_names, 3))
    tb_sigma %<>% add_column(type = c(rep(conditions[1], length(protein_names)),
                                      rep(conditions[2], length(protein_names)),
                                      rep("donor", length(protein_names))))
    tb_sigma$type %<>% factor(levels = c(conditions,"donor"))
    dodge = position_dodge(width = 0.9)
    ggplot(tb_sigma, aes(x = protein_name, y = `50%`, color = type)) +
      geom_point(size = 2, position = dodge) +
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                    position = dodge) +
      ggtitle("Marker Standard Deviations") +
      ylab("sigma") +
      theme(axis.title.y = element_blank()) +
      coord_flip()

  } else if (type == "Cor") {

    var_names = c("Cor","Cor_term","Cor_donor")
    display_names = c(conditions,"donor")
    lapply(1:length(var_names), function(i) {
      cor = rstan::extract(obj$fit_mcmc, pars = var_names[i])[[1]]
      cor_median = apply(X = cor, MARGIN = c(2,3), FUN = median)
      colnames(cor_median) = rownames(cor_median) = protein_names
      ggcorrplot(cor_median, hc.order = TRUE, type = "lower",
                 outline.col = "lightgray",
                 colors = c("#6D9EC1", "white", "#E46726")) +
        ggtitle(paste0("Marker Correlations (",display_names[i],")")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    })

  } else if (type == "Y_hat") {

    if(sum(str_detect(names(obj$fit_mcmc), type)) == 0)
      stop("Y_hat not available")

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
