#' Generates correlated binary cluster data
#'
#' Generates correrlated binary cluster data given value of Intracluster Correlation, proportion of event, perceent of variation in event proportion, number of clusters, cluster size and percent of variation in cluster size
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
#' @references Lunn, A.D. and Davies, S.J., 1998. A note on generating correlated binary variables. Biometrika, 85(2), pp.487-490.
#'
#' @author Akhtar Hossain \email{mhossain@email.sc.edu}
#' @author Hrishikesh Chakraborty \email{rishi.c@duke.edu}
#'
#' @seealso \code{\link{rcbin1}} \code{\link{iccbin}}
#'
#' @examples
#' rcbin(prop = .4, prvar = .2, noc = 30, csize = 20, csvar = .2, rho = .2)
#'
#' @importFrom stats rnorm rbinom
#'
#' @export
#'

rcbin <- function(prop = .5, prvar = 0, noc, csize, csvar = 0, rho){
  if(noc < 0) stop("The argument 'noc' should be > 0")
  if(prop < 0 || prop > 1) stop("The argument 'prob' should be in the range [0, 1]")
  if(prvar < 0 || prvar > 1) stop("The argument 'prvar' should be in the range [0, 1]")
  if(csize < 2) stop("The argument 'cssize' should be >= 2")
  if(csvar < 0) stop("The argument 'cssd' should be >= 0")
  if(rho < 0 || rho > 1) stop("The argument 'rho' should be in the range [0, 1]")
  cluster <- c(); x <- c()
  for(i in 1:noc){
    #Selecting individual cluster sizes
    min_csize <- ifelse((csize - round(csize*csvar)) >= 2, csize - round(csize*csvar), 2)
    #max_csize <- csize + round(csize*csvar)
    #Omitting to open cluster size to be as large it can be
    #Will benefit to generate skewed cluster size distribution
    csizen <- abs(round(csize + (csize*csvar)*rnorm(1)))
    while(csizen < min_csize){
      csizen <- abs(round(csize + (csize*csvar)*rnorm(1)))
    }
    #Selecting individual cluster properties
    min_prop <- ifelse((prop - prop*prvar) >= 0, prop - prop*prvar, 0)
    max_prop <- ifelse((prop + prop*prvar) <= 1, prop + prop*prvar, 1)
    propn <- abs(prop + (prop*prvar)*rnorm(1))
    while(propn < min_prop | propn > max_prop){
      propn <- abs(prop + (prop*prvar)*rnorm(1))
    }
    #Generating binary data
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
