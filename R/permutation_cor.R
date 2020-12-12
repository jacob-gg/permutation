#' permutation_cor
#'
#' Runs a non-parametric permutation test of a correlation.
#'
#' \code{permutation_cor} will:
#' \itemize{
#'   \item Run a permutation test for a correlation (using pairwise-complete observations) and report a two-sided p-value
#'   \item Report a symmetry test for the null distribution generated from permuted data (\code{symmetry.test} called from \code{lawstat})
#'   \item Produce a histogram of the generated null distribution and overlay the observed correlation on it (using \code{ggplot2})
#' }
#' Requires: \code{ggplot2} and \code{lawstat}
#'
#' @param x One of the two vectors/data frame columns to be correlated
#' @param y The second of the two vectors/data frame columns to be correlated
#' @param iterations Number of permutations (default = 10000)
#' @param method Correlation method (default = 'pearson')
#' @param seed Set for reproducibility (default = 4)
#' @return Results reported in the console; plot generated in the viewer
#' @examples
#' set.seed(1969)
#' x <- rnorm(1000, 0, 1)
#' y <- x + rnorm(1000, 0, 30)
#' x[sample(length(x), 5, replace = F)] <- NA
#' permutation_cor(x, y)
#'
#' @export
permutation_cor <- function(x, y, iterations = 10^4, method = 'pearson', seed = 4){
  # Checks
  if ('lawstat' %in% installed.packages() == FALSE) {
    stop('Error: lawstat not installed; use install.packages("lawstat")')
  }
  if ('ggplot2' %in% installed.packages() == FALSE) {
    stop('Error: ggplot2 not installed; use install.packages("ggplot2")')
  }
  if (length(x) != length(y)) {
    stop('Error: x and y must be the same length')
  }
  if (FALSE %in% unlist(lapply(list(x, y), is.numeric))) {
    stop('Error: x and y must both be numeric')
  }
  
  set.seed(seed)
  
  # Structure data, omit NAs
  if (TRUE %in% unlist(lapply(list(x, y), is.na))) {
    dat_complete <- na.omit(data.frame(x, y))
    x <- dat_complete$x
    y <- dat_complete$y
    cat('Data are not pairwise complete: Observed correlation and null distribution calculated from ', length(x), 'pairwise-complete observations', '\n')
  }
  
  # Calculate observed correlation
  observed_cor <- cor(x, y, method = method)
  
  # Permutation procedure
  null_cor_dist <- rep(NA, iterations)
  for (i in 1:iterations) {
    index <- sample(length(x), replace = FALSE)
    permuted_x_vals <- x[index]
    null_cor_dist[i] <- cor(permuted_x_vals, y, method = method)
  }
  
  # Calculate p-value: Proportion of permutations where abs(permuted correlation) > abs(observed correlation)
  cat('The observed correlation is r =', observed_cor, '\n')
  p_val <- (sum(abs(null_cor_dist) >= abs(observed_cor)) + 1) / (length(null_cor_dist) + 1)
  cat('The two-sided p-value for the permutation test is', p_val, '\n')
  
  # Symmetry test
  invisible(symmetric <- lawstat::symmetry.test(null_cor_dist, boot = FALSE))
  if (symmetric$p.value > .05) {
    cat('The null distribution for the permutation test is symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
  } else {
    cat('The null distribution for the permutation test is NOT symmetric, p =', round(symmetric$p.value, digits = 4), '\n')
    cat('Change the seed and run the test again')
  }
  
  # Histogram of null distribution
  suppressMessages(print(
    ggplot2::ggplot(data.frame(null_cor_dist)) +
      ggplot2::geom_histogram(ggplot2::aes(null_cor_dist)) +
      ggplot2::geom_vline(xintercept = observed_cor, color = 'gold') +
      ggplot2::labs(title = 'Null correlation distribution',
                    subtitle = 'Gold line is observed correlation',
                    x = 'Null distribution',
                    y = 'Frequency') +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(face = 'italic'))
  ))
  
}