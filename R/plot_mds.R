#' Plot Multivariate Posterior Summaries for Poisson Log-Normal Mixed Model
#'
#' @import rstan
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import parallel
#' @import MASS
#' @import ggrepel
#' @export
#'
#' @param obj Object of class \code{cytoeffect_poisson} computed
#'   using \code{\link{poisson_lognormal}}
#' @param asp Set \code{asp = FALSE} to avoid scaling aspect ratio by eigenvalues
#' @param ncores Number of cores
#' @param thinning Number of posterior draws after thinning
#' @param cor_scaling_factor Scaling factor for correlation arrows
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' # fit = cytoeffect::poisson_lognormal(...)
#' # cytoeffect::plot_mds(fit, asp = FALSE)
plot_mds = function(obj, asp = TRUE, ncores = parallel::detectCores(), thinning = 100,
                    cor_scaling_factor = 1, show_donors = TRUE) {

  if (class(obj) != "cytoeffect_poisson")
    stop("Not a cytoeffect_poisson object.")

  arrow_color = "darkgray"
  marker_color = "black"
  marker_size = 5
  seed = 0xdada
  stan_pars = rstan::extract(obj$fit_mcmc,
                             pars = c("beta",
                                      "sigma","sigma_term",
                                      "Q","Q_term",
                                      "b_donor"))
  condition_index = seq(obj$conditions)
  donor = obj$df_samples_subset %>%
    pull(obj$group) %>%
    as.factor %>%
    levels
  donor_index = seq(donor)
  # kth posterior draw
  sample_condition_donor = function(k, tb_info) {
    set.seed(seed)
    # fixed effects
    beta = stan_pars$beta[k,,tb_info$term_index]
    # cell random effect
    if(tb_info$term_index == 1) {
      sigma = stan_pars$sigma[k,]
      Q = stan_pars$Q[k,,]
    } else {
      sigma = stan_pars$sigma_term[k,]
      Q = stan_pars$Q_term[k,,]
    }
    Cov = Q %*% diag(sigma^2) %*% t(Q)
    # combine
    if(show_donors) {
      # donor random effect
      b_donor = stan_pars$b_donor[k,tb_info$donor_index,]
      mu = mvrnorm(n = tb_info$n, beta + b_donor, Cov)
    } else {
      mu = mvrnorm(n = tb_info$n, beta, Cov)
    }
    mu %<>% as.tibble
    names(mu) = obj$protein_names
    mu %<>% add_column(term  = tb_info$term)
    mu %<>% add_column(donor  = tb_info$donor)
    mu %<>% add_column(k  = k)
    mu
  }
  # count number of cells per term and donor
  subgroups = obj$df_samples_subset %>%
    group_by_(term = obj$condition, donor = obj$group) %>%
    tally %>%
    ungroup
  subgroups %<>% mutate(term_index = subgroups %>%
                          pull(term) %>%
                          as.integer)
  subgroups %<>% mutate(donor_index = subgroups %>%
                          pull(donor) %>%
                          as.factor %>%
                          as.integer)
  # sample one table
  sample_mu_hat = function(k) {
    lapply(seq(nrow(subgroups)), function(i) {
      sample_condition_donor(k = k, tb_info = subgroups[i,])
    }) %>% bind_rows()
  }
  # sample all tables
  sample_info_k = c("donor","term","k")
  set.seed(seed)
  expr_median = mclapply(round(seq(1, nrow(stan_pars$beta), length.out = thinning)),
                         function(i) {
                           sample_mu_hat(k = i) %>%
                             group_by(.dots = sample_info_k) %>%
                             summarize_at(obj$protein_names,median)
                         },
                         mc.cores = ncores
  ) %>% bind_rows
  # classical MDS on all posterior draws
  dist_matrix = dist(expr_median[,-seq(sample_info_k)])
  mds_res = cmdscale(dist_matrix,eig = TRUE, k = 2) # k is the number of dim
  explained_var = (100*mds_res$eig[1:2]/sum(mds_res$eig)) %>% round(digits = 1)
  expr_median %<>% bind_cols(tibble(MDS1 = mds_res$points[,1],
                                    MDS2 = mds_res$points[,2]))
  # plot MDS
  ggmds = ggplot(expr_median, aes(x = MDS1, y = MDS2, color = term)) +
    xlab(paste0("MDS1 (",explained_var[1],"%)")) +
    ylab(paste0("MDS2 (",explained_var[2],"%)")) +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A"),
                       name = obj$condition) +
    geom_density_2d()

  # make circle of correlation plot
  protein_sd = apply(as.data.frame(expr_median)[,obj$protein_names],2,sd)
  # only keep makers that have some variability
  protein_selection = obj$protein_names[protein_sd != 0]
  # correlations between variables and MDS axes
  expr_cor = cor(as.data.frame(expr_median)[,protein_selection],
                 expr_median[,c("MDS1","MDS2")]) %>% as.tibble
  # scaling factor (otherwise too crowded)
  expr_cor = expr_cor * cor_scaling_factor
  expr_cor %<>% add_column(protein_selection)
  # add arrows coordinates
  expr_cor %<>% add_column(x0 = rep(0,nrow(expr_cor)))
  expr_cor %<>% add_column(y0 = rep(0,nrow(expr_cor)))

  # add median donors
  if(show_donors) {
    expr_median_donor = expr_median %>%
      group_by(.dots = c("donor","term")) %>%
      summarize_at(c("MDS1","MDS2"), median)
    expr_median_donor %<>% add_column(
      color = sapply(expr_median_donor$term,
                     function(x) if(x == expr_median_donor$term[1]) "#5DA5DA" else "#FAA43A"))
    ggmds = ggmds +
      annotate("text",
               x = expr_median_donor$MDS1, y = expr_median_donor$MDS2,
               label = expr_median_donor$donor, color = expr_median_donor$color)
  }

  # add correlation arrows
  ggmds = ggmds + annotate("segment",
                           x = expr_cor$x0, xend = expr_cor$MDS1,
                           y = expr_cor$y0, yend = expr_cor$MDS2,
                           colour = arrow_color,
                           alpha = 1.0,
                           arrow = arrow(type = "open", length = unit(0.03, "npc")))

  ## add marker names labels
  ggmds = ggmds + geom_text_repel(data = expr_cor,
                                  aes(x = MDS1, y = MDS2,
                                      label = obj$protein_names),
                                  color = marker_color,
                                  alpha = 1.0) +
    ggtitle("Posterior MDS of Latent Variable"~mu~"(Aspect Ratio Unscaled)")

  if(asp) {
    # change aspect ratio according to explained variance
    ggmds = ggmds +
      coord_fixed(ratio = explained_var[2] / explained_var[1]) +
      ggtitle("Posterior MDS of Latent Variable"~mu)
  }
  ggmds

}
