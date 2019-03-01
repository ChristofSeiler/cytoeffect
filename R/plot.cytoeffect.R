#' Plot Posterior Summaries of Logistic Linear Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import ggcorrplot
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @param obj Object of class \code{cytoeffect} computed using \code{\link{glmm}}
#' @param type A string with the variable name to plot: \code{type = "beta"}, \code{type = "sigma_donor"}, or \code{type = "Cor_donor"}
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::glmm(...)
#' # plot(fit)
plot.cytoeffect = function(obj, type = "distribution") {

  if (class(obj) != "cytoeffect")
    stop("Not a cytoeffect object.")

  fit_mcmc = obj$fit_mcmc
  protein_names = obj$protein_names
  warmup = fit_mcmc@stan_args[[1]]$warmup
  conditions = obj$conditions

  if(type == "beta") {

    tb_beta = summary(fit_mcmc, pars = "beta", probs = c(0.025, 0.5, 0.975))
    tb_beta = tb_beta$summary[,c("2.5%","50%","97.5%")]
    tb_beta %<>% as.tibble(rownames = "name")
    tb_beta %<>% filter(name != "beta[1]") # remove intercept
    tb_beta %<>% add_column(protein_name = protein_names)
    ggplot(tb_beta, aes(x = `50%`, y = protein_name)) +
      geom_vline(xintercept = 0,color = "red") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, height = 0.3)) +
      ggtitle("Regression Coefficients") +
      xlab(conditions[2]) +
      theme(axis.title.y = element_blank())

  } else if (type == "sigma_donor") {

    tb_sigma = summary(fit_mcmc, pars = type, probs = c(0.025, 0.5, 0.975))
    tb_sigma = tb_sigma$summary[,c("2.5%","50%","97.5%")]
    tb_sigma %<>% as.tibble(rownames = "name")
    tb_sigma %<>% filter(name != "sigma_donor[1]") # remove intercept
    tb_sigma %<>% add_column(protein_name = protein_names)
    ggplot(tb_sigma, aes(x = `50%`, y = protein_name)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, height = 0.3)) +
      ggtitle("Marker Standard Deviations") +
      xlab(type) +
      theme(axis.title.y = element_blank())

  } else if (type == "Cor_donor") {

    cor = rstan::extract(fit_mcmc, pars = type)[[1]]
    cor = cor[,-1,-1] # remove intercept
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

  } else if (type == "Cor_donor_all") {

    type = strsplit(type, split = "_all")[[1]][1]
    cor = rstan::extract(fit_mcmc, pars = type)[[1]]
    cor = cor[,-1,-1] # remove intercept
    pairs = combn(protein_names, m = 2)
    cor_long = lapply(1:dim(cor)[1], function(j) {
      mcor = cor[j,,]
      rownames(mcor) = colnames(mcor) = protein_names
      tibble(
        name = sapply(1:ncol(pairs),
                      function(i) paste(pairs[,i], collapse = "-")),
        draw = as.integer(j),
        corr = sapply(1:ncol(pairs),
                      function(i) mcor[pairs[1,i], pairs[2,i]]))
    }) %>% bind_rows
    ggplot(cor_long, aes(corr, group = name)) +
      geom_line(color = "black", stat = "density", alpha = 0.3) +
      theme(legend.position="none") +
      ggtitle(paste0("Marker Correlations (",type,")"))

  # } else if (type == "uc") {

    # reordered_names = protein_names[order(summary(obj, "sigmac")$median)]
    # uc = rstan::extract(fit_mcmc,pars = "uc")[[1]]
    # uc = uc[,,-1] # remove intercept
    # # rehape posterior samples
    # convert = function(qt) {
    #   tb = as.tibble(qt)
    #   names(tb) = protein_names
    #   tb %>% mutate(celltype = celltypes)
    # }
    # tb_mid = convert(qt = apply(uc,c(2,3),median)) %>%
    #   gather(protein_name, mid, -celltype)
    # tb_low = convert(qt = apply(uc,c(2,3),
    #                             function(x) quantile(x,probs = 0.025))) %>%
    #   gather(protein_name, low, -celltype)
    # tb_high = convert(qt = apply(uc,c(2,3),
    #                              function(x) quantile(x,probs = 0.975))) %>%
    #   gather(protein_name, high, -celltype)
    # tb = bind_cols(tb_mid,tb_low,tb_high) %>%
    #   select(protein_name,celltype,mid,low,high)
    # tb$protein_name %<>% factor(levels = reordered_names)
    # tb$celltype %<>% as.factor
    # xlab_str = paste(conditions,collapse = " <-> ")
    # ggplot(tb, aes(x = mid, y = protein_name, color = protein_name)) +
    #   geom_vline(xintercept = 0,color = "black") +
    #   geom_point(size = 2) +
    #   geom_errorbarh(aes(xmin = low, xmax = high)) +
    #   ggtitle("Cell Type Random Effects") +
    #   xlab(xlab_str) +
    #   theme(legend.position="none",
    #         axis.title.y = element_blank()) +
    #   facet_wrap(~celltype,nrow = 2)

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
