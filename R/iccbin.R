#' Estimate Intracluster Correlation coefficient (ICC) and it's confidence intervals
#'
#' Estimates Intracluster Correlation coefficients in 14 different methods and it's confidence intervals in 4 different methods given the data on cluster labels and outcomes
#'
#' @param cid Column name indicating cluster id in the dataframe \code{data}
#' @param y Column name indicating binary response in the dataframe \code{data}
#' @param data A dataframe containing \code{code} and \code{y}
#' @param method The method to be used to compute ICC. A single or multiple methods can be used together. By default, all 14 methods will be used. See Details section for more.
#' @param ci.type Type of confidence interval to be computed. By default all 4 types will be reported. See Detail section for more
#' @param alpha The alpha level to be used when computing confidence interval. Default is 0.05
#' @param kappa Value of Kappa to be used in computing Stabilized ICC when the method \code{stab} is chosen. Default is 0.45
#'
#' @details If in the dataframe, the cluster id (\code{cid}) is not a factor, it will be changed to a factor and a warning message will be given
#' @details If estimate of ICC in any method is outside the interval [0, 1], the estimate and corresponding confidence interval (if appropriate) and warning messages will be produced
#' @details If the lower limit of any confidence interval is below 0 and upper limit is above 1, they will be replaced by 0 and 1 respectively and warning message will be produced
#' @details The method \code{aov} computes the analysis of variance estimate of ICC. This estimator was originally proposed for continuous variables, but various authors (e.g. Elston, 1977) have suggested it's use for binary variables
#' @details The method \code{aovs} gives estimate of ICC using a modification of analysis of variance technique (see Fleiss, 1981)
#' @details The method \code{keq} computes moment estimate of ICC suggested by Kleinman (1973), uses equal weight, \eqn{w_{i} = 1/k} for each of \eqn{k} clusters
#' @details The method \code{kpr} computes moment estimate of ICC suggested by Kleinman (1973), uses weights proportionto cluster size \eqn{w_{i} = n_{i}/N}
#' @details The method \code{keqs} gives a modified moment estimate of ICC with equal weights (\code{keq}) (see Kleinman, 1973)
#' @details The method \code{kprs} gives a modified moment estimate of ICC with weights proportional to cluster size (\code{kpr}) (see Kleinman, 1973)
#' @details The method \code{stab} provides a stabilizestimate of ICC proposed by Tamura and Young (1987)
#' @details The method \code{ub} computes moment estimate of ICC from an unbiased estimating equation (see Yamamoto and Yanagimoto, 1992)
#' @details The method \code{fc} gives Fleiss-Cuzick estimate of ICC (see Fleiss and Cuzick, 1979)
#' @details The method \code{mak} computes Mak's estimate of ICC (see Mak, 1988)
#' @details The method \code{peq} computes weighted correlation estimate of ICC proposed by karlin, Cameron, and Williams (1981) using equal weight to every pair of observations
#' @details The method \code{pgp} computes weighted correlation estimate of ICC proposed by karlin, Cameron, and Williams (1981) using equal weight to each cluster irrespective of size
#' @details The method \code{ppr} computes weighted correlation estimate of ICC proposed by karlin, Cameron, and Williams (1981) by weighting each pair according to the total number of pairs in which the individuals appear
#' @details The method \code{rm} estimates ICC using resampling method proposed by Chakraborty and Sen (2016)
#' @details The CI type \code{aov} computes confidence interval for ICC using Simith's large sample approximation (see Smith, 1957)
#' @details The CI type \code{fc} gives Fleiss-Cuzick confidence interval for ICC (see Fleiss and Cuzick, 1979)
#' @details The CI type \code{peq} estimates confidenc einterval for ICC based on direct calculation of correlation between observations within clusters (see Wu, Crespi, and Wong, 2012)
#' @details The CI type \code{rm} gives confidence interval for ICC using resampling method by Chakraborty and Sen (2016)
#'
#' @author Akhtar Hossain \email{mhossain@email.sc.edu}
#'
#' @return
#' \item{estimates}{A dataframe containing the name of methods used and corresponding estimates of Intracluster Correlation coefficients}
#' \item{ci}{A dataframe containing names of confidence interval types and corresponding estimated confidence intervals}
#'
#' @references Chakraborty, H. and Sen, P.K., 2016. Resampling method to estimate intra-cluster correlation for clustered binary data. Communications in Statistics-Theory and Methods, 45(8), pp.2368-2377.
#' @references Elston, R.C., Hill, W.G. and Smith, C., 1977. Query: Estimating" Heritability" of a dichotomous trait. Biometrics, 33(1), pp.231-236.
#' @references Fleiss, J.L., Levin, B. and Paik, M.C., 2013. Statistical methods for rates and proportions. John Wiley & Sons.
#' @references Fleiss, J.L. and Cuzick, J., 1979. The reliability of dichotomous judgments: Unequal numbers of judges per subject. Applied Psychological Measurement, 3(4), pp.537-542.
#' @references Karlin, S., Cameron, E.C. and Williams, P.T., 1981. Sibling and parent--offspring correlation estimation with variable family size. Proceedings of the National Academy of Sciences, 78(5), pp.2664-2668.
#' @references Kleinman, J.C., 1973. Proportions with extraneous variance: single and independent samples. Journal of the American Statistical Association, 68(341), pp.46-54.
#' @references Mak, T.K., 1988. Analysing intraclass correlation for dichotomous variables. Applied Statistics, pp.344-352.
#' @references Smith, C.A.B., 1957. On the estimation of intraclass correlation. Annals of human genetics, 21(4), pp.363-373.
#' @references Tamura, R.N. and Young, S.S., 1987. A stabilized moment estimator for the beta-binomial distribution. Biometrics, pp.813-824.
#' @references Wu, S., Crespi, C.M. and Wong, W.K., 2012. Comparison of methods for estimating the intraclass correlation coefficient for binary responses in cancer prevention cluster randomized trials. Contemporary clinical trials, 33(5), pp.869-880.
#' @references Yamamoto, E. and Yanagimoto, T., 1992. Moment estimators for the beta-binomial distribution. Journal of applied statistics, 19(2), pp.273-283.
#'
#' @seealso \code{\link{rcbin}}
#'
#' @examples
#' bccdata <- rcbin(prop = .4, prvar = .2, noc = 30, csize = 20, csvar = .2, rho = .2)
#' iccbin(cid = cid, y = y, data = bccdata)
#' iccbin(cid = cid, y = y, data = bccdata, method = c("aov", "fc"), ci.type = "fc")
#'
#' @importFrom stats aggregate na.omit qnorm
#'
#' @export


iccbin <- function(cid, y, data = NULL, method = c("aov", "aovs", "keq", "kpr", "keqs", "kprs", "stab", "ub", "fc", "mak",
                                                   "peq", "pgp", "ppr", "rm"), ci.type = c("aov", "fc", "peq", "rm"), alpha = 0.05, kappa = 0.45){

  ic <- list(cid = substitute(cid), y = substitute(y))

  if(is.character(ic$y)){
    # warning("Please supply either an unquoted column name of 'data' or an object for 'y'")
    if(missing(data)) stop("Supply either the unqouted name of an object containing 'y' or supply both 'data' and then 'y' as an unquoted column name to 'data'")
    ic$y <- eval(as.name(y), data, parent.frame())
  }
  if(is.name(ic$y)) ic$y <- eval(ic$y, data, parent.frame())
  if(is.call(ic$y)) ic$y <- eval(ic$y, data, parent.frame())
  if(is.character(ic$y)) ic$y <- eval(as.name(ic$y), data, parent.frame())

  if(is.character(ic$cid)){
    # warning("Please supply either an unquoted column name of 'data' or an object for 'cid'")
    if(missing(data)) stop("Supply either the unquoted name of an object containing 'cid' or supply both 'data' and then 'cid' as an unquoted column name to 'data'")
    icall$cid <- eval(as.name(cid), data, parent.frame())
  }
  if(is.name(ic$cid)) ic$cid <- eval(ic$cid, data, parent.frame())
  if(is.call(ic$cid)) ic$cid <- eval(ic$cid, data, parent.frame())
  if(is.character(ic$cid) && length(ic$cid) == 1) ic$cid <- eval(as.name(ic$cid), data, parent.frame())

  dt <- data.frame(ic)
  dt <- na.omit(dt)

  # Number off clusters
  k <- length(unique(dt$cid))

  if(!is.null(attributes(dt)$na.action)){
    warning(cat("NAs removed from data rows:\n", unclass(attributes(dt)$na.action), "\n"))
  }

  if(!is.factor(dt$cid)){
    warning("'cid' has been coerced to a factor")
    dt$cid <- as.factor(dt$cid)
  } else{
    if(length(levels(dt$cid)) > k){
      dt$x <- factor(as.character(dt$cid), levels = unique(dt$cid))
      warning("Missing levels of 'cid' have been removed")
    }
  }

  zalpha <- qnorm(alpha/2, lower.tail = F)
  square <- function(z){z^2}

  # Cluster ID
  cid <- dt$cid
  # Response variable
  y <- dt$y
  # Number off clusters
  k <- length(unique(cid))
  # Number of observations in each cluster
  ni <- as.vector(table(cid))
  # Total number of observations
  N <- sum(ni)

  meth <- c(); est <- c()
  ci.typ <- c(); uci <- c(); lci <- c()

  # ::: ANOVA method :::

  n0 <- (1/(k - 1))*(N - sum((ni^2)/N))
  yi <- aggregate(y, by = list(cid), sum)[ , 2]
  yisq <- yi^2
  msb <- (1/(k - 1))*(sum(yisq/ni) - (1/N)*(sum(yi))^2)
  msw <- (1/(N - k))*(sum(yi) - sum(yisq/ni))
  rho.aov <- (msb - msw)/(msb + (n0 - 1)*msw)
  if("aov" %in% method){
    meth <- c(meth, "Analysis of Variance Estimator")
    if(rho.aov < 0 | rho.aov > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Analysis of Variance' method")
    } else{
      est <- c(est, rho.aov)
    }
  }

  if("aov" %in% ci.type){
    ci.typ <- c(ci.typ, "Smith's large sample confidence interval")
    if(rho.aov < 0 | rho.aov > 1){
      lci <- c(lci, "-")
      uci <- c(uci, "-")
      warning("Smith's large sample confidence interval for ICC is not estimable")
    } else{
      st0 <- 2*square(1 - rho.aov)/square(n0)
      st1 <- square(1 + rho.aov*(n0 - 1))/(N - k)
      st2 <- ((k - 1)*(1 - rho.aov)*(1 + rho.aov*(2*n0 - 1)) +
                square(rho.aov)*(sum(ni^2) - (2/N)*sum(ni^3) + (1/N^2)*square(sum(ni^2))))/square(k - 1)
      var.smith.rho.aov <- st0*(st1 + st2)
      ci.smith.rho.aov <- c(rho.aov - zalpha*sqrt(var.smith.rho.aov), rho.aov + zalpha*sqrt(var.smith.rho.aov))
      lci <- c(lci, ifelse(ci.smith.rho.aov[1] < 0, 0, ci.smith.rho.aov[1]))
      uci <- c(uci, ifelse(ci.smith.rho.aov[2] > 1, 1, ci.smith.rho.aov[2]))
      if(ci.smith.rho.aov[1] < 0 | ci.smith.rho.aov[2] > 1){
        warning("One or both of 'Smith's' confidence limits fell outside of [0, 1]")
      }
    }
  }

  # ::: Modified ANOVA (Fleiss, 1981) :::

  if("aovs" %in% method){
    meth <- c(meth, "Modified Analysis of Variance Estimator (Fleiss, 1981)")
    n0 <- (1/(k - 1))*(N - sum((ni^2)/N))
    yi <- aggregate(y, by = list(cid), sum)[ , 2]
    yisq <- yi^2
    msbs <- (1/(k))*(sum(yisq/ni) - (1/N)*(sum(yi))^2)
    msw <- (1/(N - k))*(sum(yi) - sum(yisq/ni))
    rho.aovs <- (msbs - msw)/(msbs + (n0 - 1)*msw)
    if(rho.aovs < 0 | rho.aovs > 1){
      est <- c(est, "-")
      warning("ICC Not Estimable by 'Modified Analysis of Variance' method")
    } else{
      est <- c(est, rho.aovs)
    }
  }

  # ::: Moment estimator :::
  pii <- yi/ni
  wi <- rep(1/k, k)
  piw <- sum(wi*pii)
  sw <- sum(wi*(pii - piw)^2)

  # First type; Kleinman (1973); Ridout, Dometrio et al. (1999)
  if("keq" %in% method){
    meth <- c(meth, "Moment Estimator with Equal Weights")
    pii <- yi/ni
    wi <- rep(1/k, k)
    piw <- sum(wi*pii)
    sw <- sum(wi*(pii - piw)^2)
    rho.keq <- (sw - piw*(1 - piw)*sum(wi*(1 - wi)/ni))/(piw*(1 - piw)*(sum(wi*(1 - wi))) - sum(wi*(1 - wi)/ni))
    if(rho.keq < 0 | rho.keq > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Moment with Equal Weights' method")
    } else{
      est <- c(est, rho.keq)
    }
  }

  # Second type; Kleinman (1973); Ridout, Dometrio et al. (1999)
  if("kpr" %in% method){
    meth <- c(meth, "Moment Estimator with Weights Proportional to Cluster Size")
    pii <- yi/ni
    wi <- ni/N
    piw <- sum(wi*pii)
    sw <- sum(wi*(pii - piw)^2)
    rho.kpr <- (sw - piw*(1 - piw)*sum(wi*(1 - wi)/ni))/(piw*(1 - piw)*(sum(wi*(1 - wi))) - sum(wi*(1 - wi)/ni))
    if(rho.kpr < 0 | rho.kpr > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Moment with Weights Proportional to Cluster Size' method")
    } else{
      est <- c(est, rho.kpr)
    }
  }

  # Third type; Kleinman (1973); Ridout, Dometrio et al. (1999)
  if("keqs" %in% method){
    meth <- c(meth, "Modified Moment Estimator with Equal Weights")
    pii <- yi/ni
    wi <- rep(1/k, k)
    piw <- sum(wi*pii)
    sw <- sum(wi*(pii - piw)^2)
    swn <- (k - 1)*sw/k
    rho.keqs <- (swn - piw*(1 - piw)*sum(wi*(1 - wi)/ni))/(piw*(1 - piw)*(sum(wi*(1 - wi))) - sum(wi*(1 - wi)/ni))
    if(rho.keqs < 0 | rho.keqs > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Modified Moment with Equal Weights' method")
    } else{
      est <- c(est, rho.keqs)
    }
  }

  # Fourth type; Kleinman (1973); Ridout, Dometrio et al. (1999)
  if("kprs" %in% method){
    meth <- c(meth, "Modified Moment Estimator with Weights Proportional to Cluster Size")
    pii <- yi/ni
    wi <- ni/N
    piw <- sum(wi*pii)
    sw <- sum(wi*(pii - piw)^2)
    swn <- (k - 1)*sw/k
    rho.kprs <- (swn - piw*(1 - piw)*sum(wi*(1 - wi)/ni))/(piw*(1 - piw)*(sum(wi*(1 - wi))) - sum(wi*(1 - wi)/ni))
    if(rho.kprs < 0 | rho.kprs > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Modified Moment with Weights Proportional to Cluster Size' method")
    } else{
      est <- c(est, rho.kprs)
    }
  }

  # Fifth type; Tamura & Young (1987)
  if("stab" %in% method){
    meth <- c(meth, "Stabilized Moment Estimator (Tamura & Young, 1987)")
    n0 <- (1/(k - 1))*(N - sum((ni^2)/N))
    kappa = 0.45
    p <- sum(yi)/sum(ni)
    wi <- ni/N
    sw <- sum(wi*(pii - piw)^2)
    rho.stab <- (1/(n0 - 1))*((N*sw)/((k - 1)*p*(1 - p)) + kappa - 1)
    if(rho.stab < 0 | rho.stab > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Stabilized Moment' method")
    } else{
      est <- c(est, rho.stab)
    }
  }

  # Sixth type; Yamamoto & Yanagimoto (1992)
  if("ub" %in% method){
    meth <- c(meth, "Moment Estimator from Unbiased Estimating Equation")
    n0 <- (1/(k - 1))*(N - sum((ni^2)/N))
    yi <- aggregate(y, by = list(cid), sum)[ , 2]
    yisq <- yi^2
    msw <- (1/(N - k))*(sum(yi) - sum(yisq/ni))
    rho.ub <- 1 - (N*n0*(k - 1)*msw)/(sum(yi)*(n0*(k - 1) - sum(yi)) + sum(yisq))
    if(rho.ub < 0 | rho.ub > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Unbiased Estimating Equation' method")
    } else{
      est <- c(est, rho.ub)
    }
  }

  # ::: Estimators based on a direct probabilistic method :::
  # Fleiss and Cuzick (1979) Method
  yi <- aggregate(y, by = list(cid), sum)[ , 2]
  ni <- as.vector(table(cid))
  piio <- sum(yi)/sum(ni)
  rho.fc <- 1 - (1/((N - k)*piio*(1 - piio)))*sum(yi*(ni - yi)/ni)
  if("fc" %in% method){
    meth <- c(meth, "Fleiss-Cuzick Kappa Type Estimator")
    if(rho.fc < 0 | rho.fc > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Fleiss-Cuzick Kappa' method")
    } else{
      est <- c(est, rho.fc)
    }
  }

  if("fc" %in% ci.type){
    ci.typ <- c(ci.typ, "Fleiss-Cuzick Confidence Interval")
    if(rho.fc < 0 | rho.fc > 1){
      lci <- c(lci, "-")
      uci <- c(uci, "-")
      warning("Fleiss-Cuzick confidence interval for ICC is not estimable")
    } else{
    t0.fc <- 1 - rho.fc
    t1.fc <- (1/(piio*(1 - piio)) - 6)*(sum(1/ni)/(N - k)^2) + (2*N + 4*k - (k/(piio*(1 - piio))))*(k/(N*(N - k)^2)) +
      (sum(ni^2)/(N^2*piio*(1 - piio)) - ((3*N - 2*k)*(N - 2*k)*sum(ni^2))/(N^2*(N - k)^2) - (2*N - k)/(N - k)^2)*rho.fc
    t2.fc <- ((4 - 1/(piio*(1 - piio)))*((sum(ni^2))/N^2))*rho.fc^2
    var.rho.fc <- t0.fc*(t1.fc + t2.fc)
    ci.rho.fc <- c(rho.fc - zalpha*sqrt(var.rho.fc), rho.fc + zalpha*sqrt(var.rho.fc))
    lci <- c(lci, ifelse(ci.rho.fc[1] < 0, 0, ci.rho.fc[1]))
    uci <- c(uci, ifelse(ci.rho.fc[2] > 1, 1, ci.rho.fc[2]))
    if(ci.rho.fc[1] < 0 | ci.rho.fc[2] > 1){
      warning("One or both of 'Fleiss-Cuzick' Confidence limits fell outside of [0, 1]")
    }
   }
  }

  # Mak (1988) method
  if("mak" %in% method){
    meth <- c(meth, "Mak's Unweighted Average Estimator (Mak, 1988)")
    yi <- aggregate(y, by = list(cid), sum)[ , 2]
    yisq <- yi^2
    ni <- as.vector(table(cid))
    rho.mak <- 1 - (k - 1)*sum((yi*(ni - yi))/(ni*(ni - 1)))/(sum(yisq/ni^2) + sum(yi/ni)*(k - 1 - sum(yi/ni)))
    if(rho.mak < 0 | rho.mak > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Mak's Unweighted' method")
    } else{
      est <- c(est, rho.mak)
    }
  }

  # ::: Estimators based on direct calculation of correlation within each group :::
  # wi = constant and equal weight for every pair of observations; Ridout, Demetrio, et al. (1999)

  yi <- aggregate(y, by = list(cid), sum)[ , 2]
  ni <- as.vector(table(cid))
  mu.peq <- sum((ni - 1)*yi)/sum((ni - 1)*ni)
  rho.peq <- (1/(mu.peq*(1 - mu.peq)))*(sum(yi*(yi - 1))/sum(ni*(ni - 1)) - mu.peq^2)
  if("peq" %in% method){
    meth <- c(meth, "Correlation Estimator with Equal Weight to Every Pair of Observations")
    if(rho.peq < 0 | rho.peq > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Correlation Method with Weight to Every Pair of Observations'")
    } else{
      est <- c(est, rho.peq)
    }
  }

  if("peq" %in% ci.type){
    ci.typ <- c(ci.typ, "Pearson Correlation Type Confidence Interval")
    if(rho.peq < 0 | rho.peq > 1){
      lci <- c(lci, "-")
      uci <- c(uci, "-")
      warning("Pearson Correlation type confidence interval for ICC is not estimable")
    } else{
    t0.peq <- (1 - rho.peq)/sum(ni*(ni - 1))^2
    t1.peq <- 2*sum(ni*(ni - 1)) + ((1/(piio*(1 - piio)) - 3)*sum(ni^2*(ni - 1)^2))*rho.peq
    t2.peq <- ((4 - 1/(piio*(1 - piio)))*(sum(ni*(ni - 1)^3)))*rho.peq^2
    var.rho.peq <- t0.peq*(t1.peq + t2.peq)
    ci.rho.peq <- c(rho.peq - zalpha*sqrt(var.rho.peq), rho.peq + zalpha*sqrt(var.rho.peq))
    lci <- c(lci, ifelse(ci.rho.peq[1] < 0, 0, ci.rho.peq[1]))
    uci <- c(uci, ifelse(ci.rho.peq[2] > 1, 1, ci.rho.peq[2]))
    if(ci.rho.peq[1] < 0 | ci.rho.peq[2] > 1){
      warning("One or both of 'Pearson Correlation Type' Confidence limits fell outside of [0, 1]")
    }
   }
  }

  # wi = 1/(k*ni*(ni - 1)) and equal weight for each group regardless of cluster size; Ridout, Demetrio, et al. (1999)
  if("pgp" %in% method){
    meth <- c(meth, "Correlation Estimator with Equal Weight to Each Cluster Irrespective of Size")
    yi <- aggregate(y, by = list(cid), sum)[ , 2]
    ni <- as.vector(table(cid))
    mu.pgp <- sum(yi/ni)/k
    rho.pgp <- (1/(mu.pgp*(1 - mu.pgp)))*(sum((yi*(yi - 1))/(ni*(ni - 1)))/k - mu.pgp^2)
    if(rho.pgp < 0 | rho.pgp > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Correlation Method with Equal Weight to Each Cluster Irrespective of Size'")
    } else{
      est <- c(est, rho.pgp)
    }
  }

  # wi = 1/(k*ni*(ni - 1)) and equal weight for each group regardless of cluster size; Ridout, Demetrio, et al. (1999)
  if("ppr" %in% method){
    meth <- c(meth, "Correlation Estimator with Weighting Each Pair According to Total Number of Pairs the Individuals Appear")
    yi <- aggregate(y, by = list(cid), sum)[ , 2]
    ni <- as.vector(table(cid))
    mu.ppr <- sum(yi)/N
    rho.ppr <- (1/(mu.ppr*(1 - mu.ppr)))*(sum(yi*(yi - 1)/(ni - 1))/N - mu.ppr^2)
    if(rho.ppr < 0 | rho.ppr > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Correlation Method with Weighting Each Pair According to Total Number of Pairs the Individuals Appear'")
    } else{
      est <- c(est, rho.ppr)
    }
  }

  # ::: Estimators using ressampling method; Chakraborty & Sen (2016) :::
  # U ststistics
  yi <- aggregate(y, by = list(cid), sum)[ , 2]
  ni <- as.vector(table(cid))
  n <- sum(ni)
  u1 <- sum(yi)/N
  alp <- u1
  # Within cluster paire probabilities
  ucid <- sort(unique(cid))
  nw11 <- 0; nw10 <- 0; nw01 <- 0; nw00 <- 0
  for(i in 1:k){
    dti <- dt[cid == ucid[i], ]
    wsamp1 <- c()
    wsamp2 <- c()
    for(m in 1:(nrow(dti) -1)){
      wsamp1 <- c(wsamp1, rep(dti$y[m], length((m + 1):nrow(dti))))
      wsamp2 <- c(wsamp2, dti$y[(m + 1):nrow(dti)])
    }
    wsamp <- rbind(wsamp1, wsamp2)
    for(j in 1:ncol(wsamp)){
      if(all(wsamp[ , j] == c(0, 0)) == TRUE){
        nw00 <- nw00 + 1}
      else if(all(wsamp[ , j] == c(0, 1)) == TRUE){
        nw01 <- nw01 + 1}
      else if(all(wsamp[ , j] == c(1, 0)) == TRUE){
        nw10 <- nw10 + 1}
      else{nw11 <- nw11 + 1}
    }
  }
  nw11n <- nw11*2; nw00n <- nw00*2
  nw10n <- nw10 + nw01; nw01n <- nw01 + nw10
  uw11 <- nw11n/sum(ni*(ni - 1))
  uw10 <- nw10n/sum(ni*(ni - 1))
  uw01 <- nw01n/sum(ni*(ni - 1))
  uw00 <- nw00n/sum(ni*(ni - 1))
  tw <- uw11 + uw00 - uw10 - uw01

  # Between cluster paire probabilities
  nb11 <- 0; nb10 <- 0; nb01 <- 0; nb00 <- 0
  for(i in 1:(k - 1)){
    dti <- dt[cid == ucid[i], ]
    for(m in (i + 1):k){
      dtm <- dt[cid == ucid[m], ]
      bsamp1 <- rep(dti$y, each = nrow(dtm))
      bsamp2 <-rep(dtm$y, times = nrow(dti))
      bsamp <- rbind(bsamp1, bsamp2)
      for(j in 1:ncol(bsamp)){
        if(all(bsamp[ , j] == c(0, 0)) == TRUE){
          nb00 <- nb00 + 1}
        else if(all(bsamp[ , j] == c(0, 1)) == TRUE){
          nb01 <- nb01 + 1}
        else if(all(bsamp[ , j] == c(1, 0)) == TRUE){
          nb10 <- nb10 + 1}
        else{nb11 <- nb11 + 1}
      }
    }
  }
  nb11n <- nb11*2; nb00n <- nb00*2
  nb10n <- nb10 + nb01; nb01n <- nb01 + nb10
  ub11 <- nb11n/(N*(N - 1) - sum(ni*(ni - 1)))
  ub10 <- nb10n/(N*(N - 1) - sum(ni*(ni - 1)))
  ub01 <- nb01n/(N*(N - 1) - sum(ni*(ni - 1)))
  ub00 <- nb00n/(N*(N - 1) - sum(ni*(ni - 1)))
  tb <- ub11 + ub00 - ub10 - ub01

  rho.rm <- (tw - tb)/(4*u1*(1 - u1))
  if("rm" %in% method){
    meth <- c(meth, "Resampling Estimator (Chakraborty & Sen, 2016)")
    if(rho.rm < 0 | rho.rm > 1){
      est <- c(est, "-")
      warning("ICC not estimable by 'Resampling' method")
    } else{
      est <- c(est, rho.rm)
    }
  }

  if("rm" %in% ci.type){
    ci.typ <- c(ci.typ, "Resampling based Confidence Interval")
    if(rho.rm < 0 | rho.rm > 1){
      lci <- c(lci, "-")
      uci <- c(uci, "-")
      warning("Resampling based confidence interval for ICC is not estimable")
    } else{
    t0.rm <- 1/(16*n*square(u1)*square(1 - u1))
    t1.rm <- 1/(n^2 - sum(ni^2)) + square(2*alp*(1 - alp) + 1)
    t2.rm <- (alp*(1  -alp)/(sum(ni^2) - n))*((1 + alp - rho.rm*alp)*(alp + rho.rm*(1 - alp)) +
                                                (1 - alp + rho.rm*alp)*(2 - alp - rho.rm*(1 - alp)) + 2*(1 + rho.rm)*(1 - alp*(1 - alp)*(1 + rho.rm)))
    t3.rm <- ((alp*(1 - alp)*square(tw - tb)*square(1 - 2*u1))/(square(u1*(1 - u1))))*(1/n + (rho.rm/n^2)*sum(ni*(1 - ni)))
    var.rho.rm <- t0.rm*(t1.rm + t2.rm + t3.rm)
    ci.rho.rm <- c(rho.rm - zalpha*sqrt(var.rho.rm), rho.rm + zalpha*sqrt(var.rho.rm))
    lci <- c(lci, ifelse(ci.rho.rm[1] < 0, 0, ci.rho.rm[1]))
    uci <- c(uci, ifelse(ci.rho.rm[2] > 1, 1, ci.rho.rm[2]))
    if(ci.rho.rm[1] < 0 | ci.rho.rm[2] > 1){
      warning("One or both of 'Resampling Based' Confidence limits fell outside of [0, 1]")
    }
   }
  }
  estimates <- data.frame(Methods = meth, ICC = est); row.names(estimates) <- NULL
  ci <- data.frame(Type = ci.typ, LowerCI = lci, UpperCI = uci); row.names(ci) <- NULL
  list(estimates = estimates, ci = ci)
} # End of function ICCbin
