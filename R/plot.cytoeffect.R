#' Plot posterior summaries.
#'
#' @import rstan
#' @import ggplot2
#' @import ggcorrplot
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import ggthemes
#' @export
#'
plot.cytoeffect = function(obj, type = "distribution") {

  if (class(obj) != "cytoeffect")
    stop("Not a cytoeffect object.")

  fit_mcmc = obj$fit_mcmc
  protein_names = obj$protein_names
  warmup = fit_mcmc@stan_args[[1]]$warmup
  conditions = obj$conditions
  celltypes = obj$celltypes

  if(type == "beta") {

    tb_beta = summary(obj, "beta")
    ind = sort.int(tb_beta$median,index.return=TRUE)$ix
    reordered_names = tb_beta$protein_name[ind]
    tb_beta$protein_name %<>% factor(levels = reordered_names)
    xlab_str = paste(conditions,collapse = " <-> ")
    ggplot(tb_beta, aes(x = median, y = protein_name)) +
      geom_vline(xintercept = 0,color = "red") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = low, xmax = high)) +
      ggtitle("Fixed Effects") +
      xlab(xlab_str) +
      theme(axis.title.y = element_blank()) +
      theme_few()

  } else if (type == "sigma") {

    tb_random_donor = summary(obj,"sigmad") %>% mutate(effect = "donor")
    tb_random_celltype = summary(obj,"sigmac") %>% mutate(effect = "celltype")
    ind = sort.int(tb_random_celltype$median,index.return=TRUE)$ix
    reordered_names = tb_random_celltype$protein_name[ind]
    tb_random = bind_rows(tb_random_celltype, tb_random_donor)
    tb_random$protein_name %<>% factor(levels = reordered_names)
    ggplot(tb_random, aes(x = median, y = protein_name, color = effect)) +
      geom_errorbarh(aes(xmin = low, xmax = high)) +
      geom_point(size = 2) +
      ggtitle("Standard Deviations") +
      theme(axis.title.y = element_blank()) +
      facet_wrap(~effect, ncol = 2) +
      scale_colour_few()

  } else if (type == "Corc" || type == "Cord") {

    if(type == "Corc") {
      title_str = "Cell Correlations"
    } else {
      title_str = "Donor Correlations"
    }
    # compute median correlations
    cor = rstan::extract(fit_mcmc, pars = type)[[1]]
    cor = cor[,-1,-1] # remove intercept
    cor_median = apply(X = cor, MARGIN = c(2,3), FUN = median)
    cor_pos = apply(X = cor, MARGIN = c(2,3), FUN = function(x) mean(x > 0))
    cor_neg = apply(X = cor, MARGIN = c(2,3), FUN = function(x) mean(x < 0))
    p_mat = 1-pmax(cor_pos, cor_neg)
    colnames(cor_median) = rownames(cor_median) = protein_names
    ggcorrplot(cor_median, hc.order = TRUE, type = "lower",
               outline.col = "lightgray",
               colors = c("#6D9EC1", "white", "#E46726"),
               p.mat = p_mat, insig = "blank") +
      ggtitle(title_str) +
      theme_few()

  } else if (type == "Corc_all" || type == "Cord_all") {

    if(type == "Corc_all") {
      title_str = "Cell Correlations"
    } else {
      title_str = "Donor Correlations"
    }
    type = strsplit(type, split = "_")[[1]][1]
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
      ggtitle(title_str) +
      theme_few()

  } else if (type == "uc") {

    reordered_names = protein_names[order(summary(obj, "sigmac")$median)]
    uc = rstan::extract(fit_mcmc,pars = "uc")[[1]]
    uc = uc[,,-1] # remove intercept
    # rehape posterior samples
    convert = function(qt) {
      tb = as.tibble(qt)
      names(tb) = protein_names
      tb %>% mutate(celltype = celltypes)
    }
    tb_mid = convert(qt = apply(uc,c(2,3),median)) %>%
      gather(protein_name, mid, -celltype)
    tb_low = convert(qt = apply(uc,c(2,3),
                                function(x) quantile(x,probs = 0.025))) %>%
      gather(protein_name, low, -celltype)
    tb_high = convert(qt = apply(uc,c(2,3),
                                 function(x) quantile(x,probs = 0.975))) %>%
      gather(protein_name, high, -celltype)
    tb = bind_cols(tb_mid,tb_low,tb_high) %>%
      select(protein_name,celltype,mid,low,high)
    tb$protein_name %<>% factor(levels = reordered_names)
    tb$celltype %<>% as.factor
    xlab_str = paste(conditions,collapse = " <-> ")
    ggplot(tb, aes(x = mid, y = protein_name, color = protein_name)) +
      geom_vline(xintercept = 0,color = "gray") +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = low, xmax = high)) +
      ggtitle("Cell Type Random Effects") +
      xlab(xlab_str) +
      theme(legend.position="none",
            axis.title.y = element_blank()) +
      facet_wrap(~celltype,nrow = 2) +
      scale_colour_hue()

  } else {

    stop("Plotting for this type is not yet implemented.")

  }

}
