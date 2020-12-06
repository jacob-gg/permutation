#' permutation_mdiff
#'
#' Runs a non-parametric permutation test of the mean difference between two groups.
#'
#' \code{permutation_mdiff} will:
#' \itemize{
#'   \item Run a permutation test for the mean difference between two groups and report a two-sided p-value
#'   \item Report a symmetry test for the null distribution generated from permuted data (\code{symmetry.test} called from \code{lawstat})
#'   \item Produce a histogram of the generated null distribution and overlay the observed mean difference on it (using \code{ggplot2})
#' }
#' Requires: \code{ggplot2} and \code{lawstat}
#'
#' @param groups The vector/data frame column containing the (two) group identifiers
#' @param values The vector/data frame column containing the values to be compared (numeric or numeric-coercible)
#' @param iterations Number of permutations (default = 10000)
#' @param seed Set for reproducibility (default = 4)
#' @return Results reported in the console; plot generated in the viewer
#' @examples
#' set.seed(1969)
#' groups <- c(rep('a', 500), rep('b', 500))
#' values <- c(rnorm(500, 0, 1), rnorm(494, 0.2, 1), rep(NA, 6))
#' permutation_mdiff(groups, values)
#'
#' @export
permutation_mdiff <- function(groups, values, iterations = 10^4, seed = 4, one.sided = FALSE){
  
  if('lawstat' %in% installed.packages() == FALSE){
    stop('Error: lawstat not installed; use install.packages("lawstat")')
  }
  
  if('ggplot2' %in% installed.packages() == FALSE){
    stop('Error: ggplot2 not installed; use install.packages("ggplot2")')
  }
  
  if (length(groups) != length(values)) {
    stop('Error: Groups vector and values vector must be the same length')
  }
  
  if(length(unique(groups)) != 2){
    stop('Error: Ensure that there are only two unique groups; the grouping variable comes first in the list of arguments')
  }
  
  set.seed(seed)
  
  test_data <- data.frame(groups, values)
  test_data$groups <- as.factor(as.character(test_data$groups))
  test_data$values <- as.numeric(as.character(test_data$values))
  if (TRUE %in% is.na(test_data)) {
    test_data <- na.omit(test_data)
    cat('Data are incomplete: Observed mean difference and null distribution calculated from ', nrow(test_data), 'complete observations', '\n')
  }
  
  mean_table <- aggregate(test_data$values, by = list(test_data$groups), FUN = mean)
  colnames(mean_table) <- c('Group', 'Mean')
  if(mean_table$Mean[1] > mean_table$Mean[2]){
    observed_dif <- mean_table$Mean[1] - mean_table$Mean[2]
  } else {
    observed_dif <- mean_table$Mean[2] - mean_table$Mean[1]
  }
  null_mean_dif_dist <- rep(NA, iterations)
  for(i in 1:iterations) {
    index <- sample(nrow(test_data), size = nrow(test_data) / 2, replace = FALSE)
    null_mean_dif_dist[i] <- mean(test_data[index, c('values')]) - mean(test_data[-index, c('values')])
  }
  print(mean_table, row.names = FALSE)
  cat('The observed mean difference is ', observed_dif, '\n', sep = '')
  p_less <- (sum(null_mean_dif_dist <= observed_dif) + 1) / (length(null_mean_dif_dist) + 1)
  p_greater <- (sum(null_mean_dif_dist >= observed_dif) + 1) / (length(null_mean_dif_dist) + 1)
  if (p_less < p_greater) {
    p_val <- p_less*2
  } else {
    p_val <- p_greater*2
  }
  cat('The two-sided p-value for the permutation test is', p_val, '\n')
  invisible(symmetric <- lawstat::symmetry.test(null_mean_dif_dist, boot = FALSE))
  if (symmetric$p.value > .05) {
    cat('The null distribution for the permutation test is symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
  } else {
    cat('The null distribution for the permutation test is NOT symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
    cat('Change the seed and run the test again')
  }
  
  suppressMessages(print(
    ggplot2::ggplot(data.frame(null_mean_dif_dist)) +
      ggplot2::geom_histogram(ggplot2::aes(null_mean_dif_dist)) +
      ggplot2::geom_vline(xintercept = observed_dif, color = 'gold') +
      ggplot2::labs(title = 'Null mean-difference distribution',
                    subtitle = 'Gold line is observed mean difference',
                    x = 'Null distribution',
                    y = 'Frequency') +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(face = 'italic'))
  ))
  
}