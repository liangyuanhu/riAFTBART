#' Plot the propensity score by treatment
#'
#' This function estimates the propensity score for each treatment group and then plot the propensity score by each treatment to check covariate overlap.
#'
#' @param trt A numeric vector representing the treatment groups.
#' @param X A dataframe or matrix, including all the covariates but not treatments, with  rows corresponding to observations and columns to variables.
#' @param cluster.id A vector of integers representing the clustering id. The cluster id should be an integer and start from 1.
#' @param method A character indicating how to estimate the propensity score. The default is "Multinomial", which uses multinomial regression to estimate the propensity score.
#'
#' @return A plot
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' library(riAFTBART)
#' set.seed(20181223)
#' n = 5       # number of clusters
#' k = 50      # cluster size
#' N = n*k     # total sample size
#' cluster.id = rep(1:n, each=k)

#' tau.error = 0.8
#' b = stats::rnorm(n, 0, tau.error)

#' alpha = 2
#' beta1 = 1
#' beta2 = -1
#' sig.error = 0.5
#' censoring.rate = 0.02

#' x1 = stats::rnorm(N,0.5,1)
#' x2 = stats::rnorm(N,1.5,0.5)
#' trt.train = sample(c(1,2,3), N, prob = c(0.4,0.3,0.2), replace = TRUE)
#' plot_gps(trt = trt.train, X = cbind(x1, x2), cluster.id = cluster.id)
plot_gps <- function(trt, X, cluster.id, method = "Multinomial"){
  group <- NULL
  ps <- NULL
  es.max.ATE <- NULL
  if (method == "Multinomial"){
    multinom_result <- nnet::multinom(trt ~ ., data =   as.data.frame(cbind(trt, X, cluster.id)) %>% mutate(cluster.id = as.factor(cluster.id))) # Fit the multinomial regression using nnet::multinom
    pred_ps <- stats::fitted(multinom_result) # Get the fitted propensity score from the multinomial regression
    p_1 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 1) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 1|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # Figure for P(A = 1|X, V)
    p_2 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 2) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 2|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))# Figure for P(A = 2|X, V)
    p_3 <- pred_ps %>%
      cbind(trt) %>%
      as.data.frame() %>%
      tidyr::gather(group, ps,-trt) %>%
      dplyr::filter(group == 3) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x = trt, y = ps))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 3|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))# Figure for P(A = 3|X, V)
    (p <- cowplot::plot_grid(p_1, p_2, p_3,ncol = 3,align = "h")) # Combined figure
    return(p)
  } else if (method == "GBM"){
    X <- as.data.frame(X)
    temp<- noquote(names(X))
    strFormula  = sprintf("trt~%s", paste(temp, sep = "",collapse="+"))
    psmod <- twang::mnps(stats::as.formula(strFormula),
                         data=as.data.frame(cbind(trt, X, cluster.id)) %>% dplyr::mutate(trt = as.factor(trt), cluster.id = as.factor(cluster.id)), estimand = "ATE") # Fit the GBM using twang::mnps function
    p_1 <- psmod$psList$`1`$ps %>%
        cbind(trt) %>%
        dplyr::mutate(trt = as.factor(trt)) %>%
        ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
        ggplot2::geom_boxplot()+
        ggplot2::labs(x = "", y = "", title = "P(A = 1|X, V)")+
        ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
        ggplot2::theme_bw()+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))# Figure for P(A = 1|X, V)
    p_2 <- psmod$psList$`2`$ps %>%
      cbind(trt) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 2|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))# Figure for P(A = 2|X, V)
    p_3 <- psmod$psList$`3`$ps %>%
      cbind(trt) %>%
      dplyr::mutate(trt = as.factor(trt)) %>%
      ggplot2::ggplot(ggplot2::aes(x= trt,y = es.max.ATE)) +
      ggplot2::geom_boxplot()+
      ggplot2::labs(x = "", y = "", title = "P(A = 3|X, V)")+
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))# Figure for P(A = 3|X, V)
    (p <- cowplot::plot_grid(p_1, p_2, p_3,ncol = 3,align = "h")) # Combined figure
    return(p)
  }
}
