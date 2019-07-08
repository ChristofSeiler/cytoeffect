#' Plot Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import ggcorrplot
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import tibble
#' @export
#'
#' @aliases plot.cytoeffect_poisson
#' @method plot cytoeffect_poisson
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param type A string with the parameter to plot:
#'   \code{type = "beta"}, \code{type = "sigma"}, or \code{type = "Cor"}
#' @param selection A vector of strings with a selection of protein names to plot
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal(...)
#' # plot(fit)
plot.cytoeffect_poisson = function(obj, type = "beta",
                                   selection = obj$protein_names) {

  warmup = obj$fit_mcmc@stan_args[[1]]$warmup
  protein_names = obj$protein_names
  conditions = levels(pull(obj$df_samples_subset, obj$condition))

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

    to_tb = function(par_sigma, par_q, type_name) {
      sigma = rstan::extract(obj$fit_mcmc, pars = par_sigma)[[1]]
      Q = rstan::extract(obj$fit_mcmc, pars = par_q)[[1]]
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
    tb_donor = to_tb(par_sigma = "sigma", par_q = "Q", type_name = conditions[1])
    tb_cond1 = to_tb(par_sigma = "sigma_term", par_q = "Q_term", type_name = conditions[2])
    tb_cond2 = to_tb(par_sigma = "sigma_donor", par_q = "Q_donor", type_name = "donor")
    tb_sigma = bind_rows(tb_donor, tb_cond1, tb_cond2)
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

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
