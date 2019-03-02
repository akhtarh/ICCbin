
#' Generates correlated binary cluster data
#'
#' Generates correrlated binary cluster data given value of Intracluster Correlation, proportion of event and it's variance, number of clusters, cluster size and it's variance, and minimum cluster size
#' @param prop A numeric value between 0 and 1 denoting assumed proportion of event in interest, default value is 0.5. See Detail
#' @param prvar A numeric value between 0 and 1 denoting varince in assumed proportion of event (\code{prvar}), default value is 0. See Detail
#' @param noc A positive numeric value telling the number of clusters to be generated
#' @param csize A numeric value (\eqn{\ge 2}) denoting cluster size desired
#' @param csvar A positive numeric value denoting Variance of cluster size, default value is 0, see Detail
#' @param mincsize A numeric value (\eqn{\ge 2}) denoting the minimum cluster size desired, default value is 2, see Detail
#' @param rho A numeric value between 0 and 1 denoting desired level of Intracluster Correlation
#'
#' @details If supplied value of \code{prvar} is 0, the event proportion for all clusters is considered constant as supplied by \code{prop}.
#' If supplied \code{prvar} is > 0, cluster specific event proportions are generated from Beta distribution with
#' \code{shape1} and \code{shape2} parameters \eqn{a} and \eqn{b} respectively, see \code{\link{rbeta}}
#' The shape parameters are obtained using supplied values of \code{prop} and \code{prvar} by solving the equations
#'  \code{prop} \eqn{= a/(a + b)} and \code{prvar} \eqn{= ab/[(a + b)^2(1 + a + b)]}
#' @details If supplied value of \code{csvar} is 0, cluster of equal size (\code{csize}) will be generated. For \code{csvar} > 0, will be generated from
#' Normal or Negative Binomial dsitributions depending on relationship between \code{csize} and \code{csvar}.
#' If \code{csvar} < \code{csize}, the varying cluster sizes will be generated
#' from a Normal distribution with mean = \code{csize} and variacne = \code{csvar} (see \code{\link{rnorm}}).
#' If \code{csvar} \eqn{\ge} \code{csize} i.e. in the case of overdispersion,
#' cluster sizes will be generated from Negative Bionomial distribution using \code{mu} = \code{csize} and
#' \ifelse{latex}{\eqn{\code{size} = \code{csize}/[\code{csize}(\code{cscv}^2 - 1)]}}{\code{size} = \code{csize}/(\code{csize}*(\code{cscv}^2 - 1))}
#' (see \code{\link{rnbinom}}), where \code{cscv} is the coefficient of variation of cluster sizes defined as
#' \ifelse{latex}{\eqn{\sqrt{\code{csvar}}/\code{csize}}}{\code{sqrt(csvar)}/\code{csize}}. If the size of any cluster
#' is generated as less than 2, it will be replaced by the supplied value of minimum cluster size (\code{mincsize}) which has a default value
#' of 2
#'
#' @return A dataframe with two columns presenting cluster id (\code{cid}) and a binary response (\code{y}) variables
#'
#' @references Lunn, A.D. and Davies, S.J., 1998. A note on generating correlated binary variables. Biometrika, 85(2), pp.487-490.
#'
#' @author Akhtar Hossain \email{mhossain@email.sc.edu}
#'
#' @seealso \code{\link{rcbin}} \code{\link{iccbin}}
#'
#' @examples
#' rcbin1(prop = .6, prvar = .1, noc = 100, csize = 10, csvar = 12, rho = 0.2, mincsize = 2)
#'
#' @importFrom stats rnorm rbinom rbeta rnbinom
#'
#' @export
#'


rcbin1 <- function(prop = .5, prvar = 0, noc, csize, csvar = 0, mincsize = 2, rho){
  if(noc < 0) stop("The argument 'noc' should be > 0")
  if(prop < 0 || prop > 1) stop("The argument 'prob' should be in the range [0, 1]")
  if(prvar < 0 || prvar > 1) stop("The argument 'prvar' should be in the range [0, 1]")
  if(csize < 2) stop("The argument 'csize' should be >= 2")
  if(mincsize < 2) stop("The argument 'mincsize' should be >= 2")
  if(csvar < 0) stop("The argument 'cssd' should be >= 0")
  if(rho < 0 || rho > 1) stop("The argument 'rho' should be in the range [0, 1]")
  cssd <- sqrt(csvar)
  csmean <- csize
  #The case of equal cluster sizes
  if(csvar == 0){csizen <- rep(csmean, noc)}
  #The case of overdispersed cluster sizes
  if(csvar >= csmean && csvar > 0){
    csmean_n <- csmean
    #Getting CV for cluster sizes (csvar >= csmean)
    cscv <- cssd/csmean_n
    #Computing size for Negative Binomial Distribution
    r = csmean_n/(csmean_n*cscv^2 - 1)
    #Computing probability for Negative Bionomila Distribution
    #p <- r/(r + csmean_n)
    #Generating cluster sizes from Negative Bionomial Distribution
    csizen <- rnbinom(n = noc, size = r, mu = csmean_n)
    #mu can be replaced by p in above statement
    #Ensuring minimum cluster size to be mincsize (default is 2)
  }
  #The case when csvar < csmean
  if(csvar < csmean && csvar > 0){
    csmean_n <- csmean
    #Generating cluster sizes from Norrmal Distribution
    csizen <- round(rnorm(n = noc, mean = csmean_n, sd = cssd))
    #mu can be replaced by p in above statement
    #Ensuring minimum cluster size to be mincsize (default is 2)
  }
  csizen[csizen < mincsize] <- mincsize
  #Generating cluster wise event proportions
  #Case of no variation in prop (prvar = 0)
  if(prvar == 0){propn <- rep(prop, noc)}
  #Case of positive prop variation (prvar > 0)
  #Generating cluster specific proportion values from Beta distribution
  if(prvar > 0){
   a <- prop*(prop*(1 - prop)/prvar - 1)
   b <- (1 - prop)*(prop*(1 - prop)/prvar - 1)
   propn <- rbeta(noc, a, b)
  }
  #Generating cluster binary data
  cluster <- c(); x <- c()
  for(i in 1:noc){
    ri <- sqrt(rho)
    zi <- rbinom(n = 1, size = 1, prob = propn[i])
    for(j in 1:csizen[i]){
      yij <- rbinom(n = 1, size = 1, prob = propn[i])
      uij <- rbinom(n = 1, size = 1, prob = ri)
      xij <- (1 - uij)*yij + uij*zi
      cluster <- c(cluster, i); x <- c(x, xij)
    }
  }
  cbcdata <- data.frame(cid = as.factor(cluster), y = x)
  return(cbcdata)
}
