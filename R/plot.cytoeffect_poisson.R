#' Plot Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @aliases plot.cytoeffect_poisson
#' @method plot cytoeffect_poisson
#'
#' @import ggplot2
#' @import ggcorrplot
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @import stringr
#' @import tibble
#' @importFrom rstan extract summary
#' @export
#'
#' @param x Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param type A string with the parameter to plot:
#'   \code{type = "beta"}, \code{type = "sigma"}, or \code{type = "Cor"}
#' @param selection A vector of strings with a selection of protein names to plot
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal(...)
#' # plot(fit)
plot.cytoeffect_poisson = function(x, type = "beta",
                                   selection = x$protein_names, ...) {

  fit_mcmc = x$fit_mcmc
  warmup = fit_mcmc@stan_args[[1]]$warmup
  protein_names = x$protein_names
  conditions = levels(pull(x$df_samples_subset, x$condition))
  donors = unique(pull(x$df_samples_subset, x$group))

  if(type == "beta") {

    # extract index estimates
    tb_beta = rstan::summary(fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as_tibble(rownames = "name")
    tb_beta %<>% add_column(protein_name = rep(protein_names, each = length(conditions)))
    tb_beta %<>% add_column(condition = rep(conditions, length(protein_names)))

    # add contrast
    beta = rstan::extract(fit_mcmc, pars = "beta")[[1]]
    beta_contrast = beta[,,2] - beta[,,1]
    beta_diff = apply(beta_contrast, 2,
                      function(x)
                        quantile(x, probs = c(0.025, 0.5, 0.975))
    ) %>% t %>% as_tibble
    beta_diff %<>% add_column(name = "contrast")
    beta_diff %<>% add_column(protein_name = protein_names)
    beta_diff %<>% add_column(condition = paste(rev(conditions),collapse = " - "))
    tb_beta %<>% bind_rows(beta_diff)

    # combined plot
    dodge = position_dodge(width = 0.9)
    ggplot(tb_beta %>% dplyr::filter(protein_name %in% selection),
           aes(x = protein_name, y = `50%`, color = condition)) +
      geom_hline(yintercept = 0,color = "grey") +
      geom_point(size = 2, position = dodge) +
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),
                    position = dodge) +
      ggtitle("Regression Coefficients") +
      ylab("log expected count") +
      theme(axis.title.y = element_blank()) +
      coord_flip()

  } else if (type == "sigma") {

    # summarize cell level effect
    tb_cond = rstan::summary(fit_mcmc, pars = c("sigma"),
                             probs = c(0.025, 0.5, 0.975))
    tb_cond = tb_cond$summary[,c("2.5%","50%","97.5%")]
    tb_cond %<>% as_tibble(rownames = "name")
    tb_cond %<>% add_column(protein_name = protein_names)
    tb_cond %<>% add_column(type = "cell")

    # summarize donor level effect
    to_tb = function(par_sigma, par_q, type_name) {
      sigma = rstan::extract(fit_mcmc, pars = par_sigma)[[1]]
      Q = rstan::extract(fit_mcmc, pars = par_q)[[1]]
      qvtq = lapply(1:nrow(Q), function(i) {
        cov = Q[i,,] %*% diag(sigma[i,]^2) %*% t(Q[i,,])
        sqrt(diag(cov))
      }) %>% bind_cols()
      tb = apply(
        qvtq, 1, function(x)
          quantile(x, probs = c(0.025, 0.5, 0.975))
        ) %>% t %>% as_tibble %>%
        add_column(protein_name = protein_names) %>%
        add_column(type = rep(type_name, length(protein_names)))
      tb
    }
    tb_donor = to_tb(par_sigma = "sigma_donor", par_q = "Q_donor", type_name = x$group)

    # combine and plot
    tb_sigma = bind_rows(tb_cond, tb_donor)
    tb_sigma$type %<>% factor(levels = c("cell",x$group))
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

    var_names = c("Cor","Cor_donor")
    display_names = c("cell",x$group)
    lapply(1:length(var_names), function(i) {
      cor = rstan::extract(fit_mcmc, pars = var_names[i])[[1]]
      cor_median = apply(X = cor, MARGIN = c(2,3), FUN = median)
      colnames(cor_median) = rownames(cor_median) = protein_names
      ggcorrplot(cor_median, hc.order = TRUE, type = "lower",
                 outline.color = "lightgray",
                 colors = c("#6D9EC1", "white", "#E46726")) +
        ggtitle(paste0("Marker Correlations (",display_names[i],")")) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    })

  } else if (type == "theta") {

    theta = rstan::extract(fit_mcmc, pars = "theta")[[1]]
    theta_median = apply(X = theta, MARGIN = c(2,3), FUN = median)
    colnames(theta_median) = donors
    theta_median %<>% as_tibble()
    theta_median %<>% mutate(protein_name = protein_names)
    theta_median %<>% tidyr::gather(donor, probability, -protein_name)
    ggplot(theta_median, aes(protein_name, donor)) +
      geom_tile(aes(fill = probability), color = "white") +
      scale_fill_gradient(low = "white", high = "steelblue") +
      coord_fixed() +
      theme(panel.background = element_rect(fill = "lightgray"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title = element_blank())

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
