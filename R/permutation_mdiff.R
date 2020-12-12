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
#' Note: If groups A and B are unbalanced, \code{permutation_mdiff}'s default is to generate permuted groups of sizes n_A and n_B; you can, however, use a downsampling procedure to generate permuted groups of size min(n_A, n_B); see above.
#'
#' @param groups The vector/data frame column containing the (two) group identifiers
#' @param values The vector/data frame column containing the values to be compared (numeric or numeric-coercible)
#' @param iterations Number of permutations (default = 10000)
#' @param seed Set for reproducibility (default = 4)
#' @param downsample If groups A and B are unbalanced, you can set the \code{downsample} argument to TRUE to use downsampling each permutation to generate permuted groups of size min(n_A, n_B); otherwise, \code{permutation_mdiff} will generate permuted groups of size n_A and n_B during the permutation procedure (default = FALSE)
#' @return Results reported in the console; plot generated in the viewer
#' @examples
#' set.seed(1969)
#' groups <- c(rep('control', 500), rep('experimental', 495))
#' values <- c(runif(500, 0, 5), runif(490, 0.25, 5), rep(NA, 5))
#' permutation_mdiff(groups, values)
#' permutation_mdiff(groups, values, downsample = T)
#'
#' @export
permutation_mdiff <- function(groups, values, iterations = 10^4, seed = 4, downsample = FALSE){
  # Checks
  if ('lawstat' %in% installed.packages() == FALSE) {
    stop('Error: lawstat not installed; use install.packages("lawstat")')
  }
  if ('ggplot2' %in% installed.packages() == FALSE) {
    stop('Error: ggplot2 not installed; use install.packages("ggplot2")')
  }
  if (length(groups) != length(values)) {
    stop('Error: Groups vector and values vector must be the same length')
  }
  if (length(unique(groups)) != 2) {
    stop('Error: Ensure that there are only two unique groups; the grouping variable comes first in the list of arguments')
  }
  
  set.seed(seed)
  
  # Structure data, omit NAs
  test_data <- data.frame(groups, values)
  test_data$groups <- as.factor(as.character(test_data$groups))
  test_data$values <- as.numeric(as.character(test_data$values))
  if (TRUE %in% is.na(test_data)) {
    test_data <- na.omit(test_data)
    if (downsample == F) {
      cat('Data are incomplete: Observed mean difference and null distribution calculated from', nrow(test_data), 'complete observations', '\n')
    }
    if (downsample == T) {
      cat('Data are incomplete: Observed mean difference calculated from', nrow(test_data), 'complete observations', '\n')
    }
  }
  
  # Calculate observed mean difference
  mean_table <- aggregate(test_data$values, by = list(test_data$groups), FUN = mean)
  colnames(mean_table) <- c('Group', 'Mean')
  if(mean_table$Mean[1] > mean_table$Mean[2]){
    observed_dif <- mean_table$Mean[1] - mean_table$Mean[2]
    group_one_greater <- T
  } else {
    observed_dif <- mean_table$Mean[2] - mean_table$Mean[1]
    group_one_greater <- F
  }
  
  # Permutation procedure begins
  null_mean_diff_dist <- rep(NA, iterations)
  
  # Note: In the permutation procedure below, using group_one_greater == T|F to determine whether null mean differences are calculated as a - b or b - a in the null-distribution generation isn't strictly necessary, as the two-sided p-value is calculated as the proportion of cases where the permuted absolute difference exceeds the observed absolute difference. However, given that the observed difference is set up in this function to always be positive, if the two-sided p-value was calculated as min(p_less, p_greater)*2, it would change (very) slightly depending on whether the null distribution was calculated with a - b or b - a (because the tails aren't perfectly symmetric).
  
  # Permutation procedure without downsampling (default; i.e., creating permuted groups of size n_a and n_b)
  if (downsample == F) {
    for (i in 1:iterations) {
      permuted_test_data <- test_data
      permuted_test_data$values <- sample(permuted_test_data$values, size = nrow(permuted_test_data), replace = F)
      if (group_one_greater == T) {
        null_mean_diff_dist[i] <- mean(permuted_test_data[permuted_test_data$groups == mean_table$Group[1], c('values')]) - mean(permuted_test_data[permuted_test_data$groups == mean_table$Group[2], c('values')])
      }
      if (group_one_greater == F) {
        null_mean_diff_dist[i] <- mean(permuted_test_data[permuted_test_data$groups == mean_table$Group[2], c('values')]) - mean(permuted_test_data[permuted_test_data$groups == mean_table$Group[1], c('values')]) 
      }
    }
  }
  
  # Permutation procedure with downsampling (i.e., creating permuted groups of size min(n_a, n_b) each permutation)
  if (downsample == T) {
    n_per_group <- table(test_data$groups)
    cat('Downsampling requested: Downsampling done each permutation to generate null-distribution groups of size n =', min(n_per_group), '\n')
    smaller_group_name <- names(which(n_per_group == min(n_per_group)))
    bigger_group_name <- names(which(n_per_group == max(n_per_group)))
    for (i in 1:iterations) {
      bigger_group_downsample <- sample(rownames(test_data[test_data$groups == bigger_group_name, ]), nrow(test_data[test_data$groups == smaller_group_name, ]), replace = F)
      smaller_group <- rownames(test_data[test_data$groups == smaller_group_name, ])
      downsampled_permuted_test_data <- test_data[c(bigger_group_downsample, smaller_group), ]
      downsampled_permuted_test_data$values <- sample(downsampled_permuted_test_data$values, size = nrow(downsampled_permuted_test_data), replace = F)
      if (group_one_greater == T) {
        null_mean_diff_dist[i] <- mean(downsampled_permuted_test_data[downsampled_permuted_test_data$groups == mean_table$Group[1], c('values')]) - mean(downsampled_permuted_test_data[downsampled_permuted_test_data$groups == mean_table$Group[2], c('values')])
      }
      if (group_one_greater == F) {
        null_mean_diff_dist[i] <- mean(downsampled_permuted_test_data[downsampled_permuted_test_data$groups == mean_table$Group[2], c('values')]) - mean(downsampled_permuted_test_data[downsampled_permuted_test_data$groups == mean_table$Group[1], c('values')])
      }
    }
  }
  
  # Calculate p-value: Proportion of permutations where abs(permuted mean difference) > abs(observed mean difference)
  print(mean_table, row.names = FALSE)
  cat('The observed mean difference is ', observed_dif, '\n', sep = '')
  p_val <- (sum(abs(null_mean_diff_dist) >= abs(observed_dif)) + 1) / (length(null_mean_diff_dist) + 1)
  cat('The two-sided p-value for the permutation test is', p_val, '\n')
  
  # Symmetry test
  invisible(symmetric <- lawstat::symmetry.test(null_mean_diff_dist, boot = FALSE))
  if (symmetric$p.value > .05) {
    cat('The null distribution for the permutation test is symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
  } else {
    cat('The null distribution for the permutation test is NOT symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
    cat('Change the seed and run the test again')
  }
  
  # Histogram of null distribution
  suppressMessages(print(
    ggplot2::ggplot(data.frame(null_mean_diff_dist)) +
      ggplot2::geom_histogram(ggplot2::aes(null_mean_diff_dist)) +
      ggplot2::geom_vline(xintercept = observed_dif, color = 'gold') +
      ggplot2::labs(title = 'Null mean-difference distribution',
                    subtitle = paste0('Gold line is observed mean difference, ', mean_table$Group[mean_table$Mean == max(mean_table$Mean)]
                                      , ' - ', mean_table$Group[mean_table$Mean == min(mean_table$Mean)]),
                    x = 'Null distribution',
                    y = 'Frequency') +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(face = 'italic'))
  ))
  
}