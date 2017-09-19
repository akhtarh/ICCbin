#' Generates correlated binary cluster data
#'
#' Generates correrlated binary cluster data given value of Intracluster Correlation, proportion of event and it's variation, number of clusters, cluster size and variation in cluster size
#' @param prop A numeric value between 0 to 1 denoting assumed proportion of event in interest. Default value is 0.5
#' @param prvar A numeric value between 0 to 1 denoting percent of variation in assumed proportion of event (\code{prvar}). Default value is 0
#' @param noc A numeric value telling the number of clusters to be generated
#' @param csize A numeric value denoting desired cluster size
#' @param csvar A numeric value between 0 to 1 denoting percent of variation in cluster sizes (\code{csize}). Default value is 0
#' @param rho A numeric value between 0 to 1 denoting desired level of Intracluster Correlation
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


rcbin <- function(prop = .5, prvar = 0, noc, csize, csvar = 0, rho){
  cluster <- c(); x <- c()
  for(i in 1:noc){
    # Selecting individual cluster sizes
    min_csize <- csize - round(csize*csvar)
    max_csize <- csize + round(csize*csvar)
    csizen <- abs(round(csize + (csize*csvar)*rnorm(1)))
    while(csizen < min_csize | csizen > max_csize){
      csizen <- abs(round(csize + (csize*csvar)*rnorm(1)))
    }
    # Selecting individual cluster properties
    min_prop <- prop - prop*prvar
    max_prop <- prop + prop*prvar
    propn <- abs(prop + (prop*prvar)*rnorm(1))
    while(propn < min_prop | propn > max_prop){
      propn <- abs(prop + (prop*prvar)*rnorm(1))
    }
    # Generating binary data
    ri <- sqrt(rho)
    zi <- rbinom(n = 1, size = 1, prob = propn)
    for(j in 1:csizen){
      yij <- rbinom(n = 1, size = 1, prob = propn)
      uij <- rbinom(n = 1, size = 1, prob = ri)
      xij <- (1 - uij)*yij + uij*zi
      cluster <- c(cluster, i); x <- c(x, xij)
    }
  }
  cbcdata <- data.frame(cid = as.factor(cluster), y = x)
  return(cbcdata)
}




