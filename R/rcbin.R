#' Generates correlated binary cluster data
#'
#' Generates correrlated binary cluster data given value of Intracluster Correlation, proportion of event and it's variation, number of clusters, cluster size and variation in cluster size
#' @param prop A numeric value between 0 and 1 denoting assumed proportion of event in interest, default value is 0.5. See Detail
#' @param prvar A numeric value between 0 and 1 denoting percent of variation in assumed proportion of event (\code{prvar}), default value is 0. See Detail
#' @param noc A numeric value telling the number of clusters to be generated
#' @param csize A numeric value denoting desired cluster size. See Deatil
#' @param csvar A numeric value between 0 and 1 denoting percent of variation in cluster sizes (\code{csize}), default value is 0. See Detail
#' @param rho A numeric value between 0 and 1 denoting desired level of Intracluster Correlation
#'
#' @details The minimum and maximum values of event proportion (\code{prop}) will be taken as 0 and 1 respectively in cases where it exceeds the valid limits (0, 1) due to larger value of percent variation (\code{prvar}) supplied
#' @details The minimum value of cluster size (\code{csize}) will be taken as 2 in cases where it goes below 2 due to larger value of percent variation (\code{csvar}) supplied
#'
#' @return A dataframe with two columns presenting cluster id (\code{cid}) and a binary response (\code{y}) variables
#'
#' @author Akhtar Hossain \email{mhossain@email.sc.edu}
#' @author Hrishikesh Chakraborty \email{rishi.c@duke.edu}
#'
#' @seealso \code{\link{iccbin}}
#'
#' @examples
#' rcbin(prop = .4, prvar = .2, noc = 30, csize = 20, csvar = .2, rho = .2)
#'
#' @importFrom stats rnorm rbinom
#'
#' @export
#'

rcbin <- function (prop = 0.5, prvar = 0, noc, csize, csvar = 0, rho)
{
  # Select individual cluster sizes
  cluster_sizes <- pmax(abs(round(csize + (csize * csvar * rnorm(noc)))),
                        2)
  # How many total units will we have?
  total_n <- sum(cluster_sizes)
  # Cluster level proportions
  cluster_proportions = pmin(pmax(abs(prop + (prop * prvar) * rnorm(noc)),
                                  0),
                             1)
  # Cluster ID generation
  cluster_id = rep(1:noc, times=cluster_sizes)
  # Individual proportion is just an expansion
  prop_individuals = rep(cluster_proportions, times = cluster_sizes)
  # Generate cluster level z_is and pad to be z_ijs
  z_ij = rep(rbinom(n = noc, size=1, prob = cluster_proportions),
             times = cluster_sizes)
  # Individual unit assignments
  y_ij = rbinom(n = total_n, size = 1, prob = prop_individuals)
  # Switching dummy
  u_ij = rbinom(n = total_n, size = 1, prob = sqrt(rho))
  # Assume x_ij based on switching dummy
  x_ij = ifelse(u_ij, z_ij, y_ij)

  # Return data
  data.frame(cid = as.factor(cluster_id), y = x_ij)
}
