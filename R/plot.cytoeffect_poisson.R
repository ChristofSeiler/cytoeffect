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

  if(type == "beta") {

    tb_beta = summary(fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as.tibble(rownames = "name")
    #tb_beta = tb_beta[seq(2,nrow(tb_beta),2),]
    tb_beta %<>% add_column(protein_name = rep(protein_names, each = length(covariates)))
    tb_beta %<>% add_column(covariate = rep(covariates, length(protein_names)))
    #ind = sort.int(tb_beta$`50%`,index.return=TRUE)$ix
    #reordered_names = tb_beta$covariate[ind]
    #tb_beta$covariate %<>% factor(levels = reordered_names)
    #xlab_str = paste(conditions,collapse = " <-> ")
    ggplot(tb_beta, aes(x = `50%`, y = protein_name)) +
      geom_vline(xintercept = 0,color = "red") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
      ggtitle("Regression Coefficients") +
      #xlab(xlab_str) +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
      facet_wrap(~ covariate, ncol = 4)

  } else if (type == "sigma" || type == "sigma_donor") {

    tb_sigma = summary(fit_mcmc, pars = type, probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% add_column(protein_name = protein_names)
    ind = sort.int(tb_sigma$`50%`,index.return=TRUE)$ix
    reordered_names = tb_sigma$protein_name[ind]
    tb_sigma$protein_name %<>% factor(levels = reordered_names)
    ggplot(tb_sigma, aes(x = `50%`, y = protein_name)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`)) +
      ggtitle("Marker Standard Deviations") +
      xlab(type) +
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
               colors = c("#6D9EC1", "white", "#E46726"),
               p.mat = p_mat, insig = "blank") +
      ggtitle(paste0("Marker Correlations (",type,")")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
