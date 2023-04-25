#' Test series to apply the SNHT test.
#' 
#' @description Build a test series (Q) suitable to perform the SNHT test.
#' 
#' @param cand a vector containing the candidate data series.
#' @param aux a data frame containing the auxiliary data series.
#' @param method a character string specifying the method to use for creating
#' series of anomalies, either "ratio" or "difference". Defaults to "ratio".
#' @param min.cor the minimum correlation coefficient to use for selecting
#' appropriate reference time series. If NULL, all reference series with a
#' positive correlation coefficient will be used. Defaults to NULL.
#' @param ref.period a numeric vector of length 2 specifying the start and
#' end of the reference period used to calculate the anomalies. If NULL, the
#' entire time series will be used. Defaults to NULL.
#' @param n.min the minimum number of reference time series required for
#' consensus series to be built. If the number of reference time series is
#' less than this value, an error will be thrown. If NULL, no minimum is
#' required. Defaults to NULL.
#' @param n.max the maximum number of reference time series to use for the
#' SNHT test. If the number of reference time series is greater than this
#' value, only the top `n.max` correlated reference time series will be used.
#' If NULL, all available reference time series will be used. Defaults to NULL.
#' 
#' @details Given a candidate time series and a set of auxiliary series
#' (usually from nearby observatories), this functions creates a 'consensus'
#' series from the auxiliary series and then computes anomalies (differences
#' or ratios) between the candidate and the consensus. This latter time series,
#' termed `Q`, is later used to perform the SNHT test.
#' 
#' @return a list containing the input data, the method used, the Q series,
#' and the number of auxiliary series used.
#' 
#' #' @references
#' H. Alexandersson (1986), A homogeneity test applied to precipitation data,
#' \emph{Journal of Climatology} 6, 661--675.
#' 
#' @examples
#' set.seed(0124)
#' cand <- c(rnorm(20,2,0.15), rnorm(20,4,0.15))
#' r1 <- cand + rnorm(40,1,0.15)
#' r2 <- cand + rnorm(40,1,0.45)
#' r3 <- cand + rnorm(40,1,0.35)
#' aux <- cbind(r1,r2,r3)
#' Q <- snht.Q(cand, aux, method='ratio') # min.cor=0.8, ref.period=c(1:30), n.min = 3, n.max = 10
#' 
#' @export
snht.Q <- function(cand, aux, method=c('ratio', 'diff'),
                   min.cor=0, ref.period=1:length(cand),
                   n.min=1, n.max=Inf, robust=FALSE)
{
  
  require(MASS)
  
  # input checks
  if(length(cand) != nrow(aux)){
     stop("Vectors 'candidate' and 'reference' must have the same length.")
  }
  if(is.na(match(method,c('ratio', 'diff')))){
     stop("Parameter 'method' must be either 'ratio' or 'diff'.")
  }
  # TO DO: more checks:
  #   min.cor lies between (0,1)
  #   ref.period is a vector two, both positives, (min,max)
  #   n.min is positive
  #   n.max > n.min
  
  # create y1 and y2 objects
  y1 <- cand
  y2 <- aux
  
  # sort auxiliary series by correlation
  rho <- cor(y1, y2, use='complete.obs')
  y2 <- y2[,order(rho)]
  rho <- rho[,order(rho)]
  
  # select auxiliary series by min.cor, and sort
  y2 <- y2[, which(rho >= min.cor)]
  y2 <- as.data.frame(y2)

  # check that there are at least n.min auxiliary series
  if(ncol(y2) < n.min){
     stop("There are no sufficient auxiliary series available.")
  }
  
  # cut to n.max
  if (!is.infinite(n.max)) {
    n <- min(ncol(y2), n.max)
    y2 <- y2[, 1:n]
  }
  
  # compute correlation, cut to ref.period
  rho <- cor(y1[ref.period], y2[ref.period,], use='complete.obs')
  
  # compute means, cut to ref.period
  if (!robust) {
    y1_mean  <- mean(y1[ref.period], na.rm=TRUE)
    if (ncol(y2) == 1) {
      y2_mean <- mean(y2[ref.period,], na.rm=TRUE)
    } else {
      y2_mean <- colMeans(y2[ref.period,], na.rm=TRUE)
    }
  }
  if (robust) {
    y1_mean  <- huber(y1[ref.period])$m
    if (ncol(y2) == 1) {
      y2_mean <- huber(y2[ref.period,])$m
    } else {
      y2_mean <- apply(y2[ref.period,], 2, function(x) huber(x)$m)
    }
  }

  # compute weighted anomalies (ratios or differences w.r. to the mean,
  # times the correlation coefficient)
  y2_ano <- y2 * NA
  if(method == 'ratio') {
    for(j in 1:ncol(y2)) {
      y2_ano[,j] <- rho[j]^2 * y2[,j] * y1_mean / y2_mean[j]
    }
  }
  if(method == 'diff') {
    for(j in 1:ncol(y2)){
      y2_ano[,j] <- rho[j]^2 * (y2[,j] - y2_mean[j] + y1_mean)
    }
  }
  
  # consensus series (weighted average of auxiliary anomaly series)
  y2_cons <- apply(y2_ano, MARGIN=1, FUN=sum, na.rm=TRUE) / sum(rho^2)

  # Q series: ratio or difference of candidate to consensus
  if(method == 'diff') {
    Q <- y1 - y2_cons
  }
  if(method == 'ratio') {
    Q <- y1 / y2_cons
  }
  #Q <- na.omit(cbind(x1, Q))[,1:2]

  # return results
  res <- list(cand=cand, cons=y2_cons, method=method, Q=Q, n=ncol(y2))
  return(res)
}





#' Change-point detection in a time series.
#' 
#' @description Run the Standard Normal Homogeneity Test (SNHT) to detect a change point
#' in a time series.
#' 
#' @param x a vector containing the candidate data series.
#' @param clevel a numeric value representing the confidence level of the
#' test. Must be one of c(90, 92, 94, 95, 97.5, 99). Defaults to 95.
#' @param excl numeric, exclude the first and last elements from the analysis.
#' Defaults to 5.
#' 
#' @details Applies the Standard Normal Homogeneity Test (SNHT), in the
# version from Alexandersson (1986), that is the test for (one) shift without
# trend. The function is based on an R script by Pascal Haenggi, v20100219.
#' The test can be applied to a single time series, or to a time series of
#' anomalies with respect to a consensus series formed from neighbouring
#' auxiliary series (see function `snht.Q()`).
#' The function assumes that the time series is stationary (i.e., it contains
#' no trends or large-scale cycles, nor seasonality). Diverse strategies can be
#' followed to remove trends and cycles, such as computing anomalies with 
#' respect to a consensus series, and seasonality.
#' 
#' #' @return A list containing the input data, the snht statistic time series,
#' the loation of the highest change point, the means at both sides of the
#' change point, and the significance.
#' 
#' @references
#' H. Alexandersson (1986), A homogeneity test applied to precipitation data,
#' \emph{Journal of Climatology} 6, 661--675.
#' M. N. Khaliq, T. B. M. J. Ouarda (2007), On the critical values of the
#' standard normal homogeneity test (SNHT),
#' \emph{International Journal of Climatology} 27, 681--687.
#' 
#' @examples
#' set.seed(123)
#' Q <- rnorm(1000)
#' Q[201:500] <- Q[201:500] + 0.4
#' Q[501:600] <- Q[501:600] - 0.6
#' 
#' res <- snht(Q)
#' 
#' par(mfrow=c(2,1))
#' plot(baseData, type="p", pch='.')
#' abline(h=res$left_mean, col='red')
#' abline(h=res$right_mean, col='red')
#' plot(res$Tv, type="l")
#' abline(v=res$T0x, col="Red")
#' abline(h=res$Tc, lty="dashed")
#' 
#' @export
snht <- function(x, p=NULL, clevel=95, excl=5, cl.factor=1, robust=FALSE) {
  
  require(MASS)
  
  # input checks
  if(is.na(match(clevel,c(90,92,94,95,97.5,99)))){
     stop("Parameter 'clevel' must be 90, 92, 94, 95, 97.5 or 99")
  }
  
  # critical levels (from Khaliq and Ouarda, 2007)
  cval <- matrix(data = c(10,4.964,5.197,5.473,5.637,6.188,6.769,
    12,5.288,5.554,5.876,6.068,6.729,7.459,14,5.54,5.831,6.187,6.402,7.152,8.001,
    16,5.749,6.059,6.441,6.674,7.492,8.44,18,5.922,6.248,6.652,6.899,7.775,8.807,
    20,6.07,6.41,6.83,7.089,8.013,9.113,22,6.2,6.551,6.988,7.257,8.22,9.38,
    24,6.315,6.675,7.123,7.4,8.4,9.609,26,6.417,6.785,7.246,7.529,8.558,9.812,
    28,6.509,6.884,7.353,7.643,8.697,9.993,30,6.592,6.973,7.451,7.747,8.825,10.153,
    32,6.669,7.056,7.541,7.841,8.941,10.3,34,6.741,7.132,7.625,7.93,9.05,10.434,
    36,6.803,7.201,7.699,8.009,9.143,10.552,38,6.864,7.263,7.768,8.081,9.23,10.663,
    40,6.921,7.324,7.835,8.151,9.317,10.771,42,6.972,7.38,7.894,8.214,9.39,10.865,
    44,7.022,7.433,7.951,8.273,9.463,10.957,46,7.071,7.484,8.007,8.331,9.53,11.04,
    48,7.112,7.529,8.054,8.382,9.592,11.116,50,7.154,7.573,8.103,8.432,9.653,11.193,
    52,7.194,7.616,8.149,8.48,9.711,11.259,54,7.229,7.654,8.19,8.524,9.76,11.324,
    56,7.264,7.69,8.23,8.566,9.81,11.382,58,7.299,7.727,8.268,8.606,9.859,11.446,
    60,7.333,7.764,8.308,8.647,9.906,11.498,62,7.363,7.796,8.343,8.683,9.947,11.548,
    64,7.392,7.827,8.375,8.717,9.985,11.599,66,7.421,7.857,8.408,8.752,10.026,11.648,
    68,7.449,7.886,8.439,8.784,10.067,11.692,70,7.475,7.913,8.467,8.814,10.099,11.737,
    72,7.499,7.938,8.496,8.844,10.134,11.776,74,7.525,7.965,8.523,8.873,10.171,11.822,
    76,7.547,7.989,8.548,8.898,10.2,11.858,78,7.57,8.013,8.575,8.926,10.23,11.895,
    80,7.591,8.035,8.599,8.951,10.259,11.928,82,7.613,8.059,8.623,8.976,10.29,11.966,
    84,7.634,8.079,8.647,9.001,10.315,11.995,86,7.655,8.102,8.67,9.026,10.347,12.033,
    88,7.673,8.121,8.691,9.047,10.37,12.059,90,7.692,8.14,8.71,9.067,10.394,12.089,
    92,7.711,8.16,8.732,9.09,10.417,12.12,94,7.73,8.181,8.752,9.11,10.447,12.153,
    96,7.745,8.196,8.77,9.127,10.465,12.175,98,7.762,8.214,8.788,9.147,10.484,12.196,
    100,7.778,8.231,8.807,9.167,10.507,12.228,105,7.819,8.273,8.851,9.214,10.562,12.291,
    110,7.856,8.312,8.892,9.255,10.608,12.343,115,7.891,8.35,8.931,9.296,10.656,12.401,
    120,7.921,8.38,8.963,9.33,10.694,12.446,125,7.952,8.413,8.999,9.365,10.735,12.488,
    130,7.983,8.446,9.032,9.4,10.772,12.538,135,8.01,8.474,9.063,9.431,10.808,12.579,
    140,8.038,8.501,9.092,9.462,10.845,12.621,145,8.063,8.529,9.12,9.49,10.877,12.66,
    150,8.086,8.554,9.147,9.519,10.906,12.694,155,8.111,8.578,9.172,9.543,10.933,12.725,
    160,8.133,8.601,9.195,9.569,10.966,12.759,165,8.155,8.625,9.222,9.596,10.992,12.793,
    170,8.174,8.643,9.241,9.615,11.016,12.82,175,8.195,8.666,9.265,9.641,11.046,12.851,
    180,8.214,8.685,9.283,9.658,11.062,12.872,185,8.233,8.706,9.307,9.683,11.089,12.904,
    190,8.252,8.725,9.325,9.701,11.11,12.93,195,8.268,8.741,9.343,9.72,11.132,12.956,
    200,8.286,8.761,9.364,9.741,11.156,12.982,225,8.361,8.838,9.446,9.826,11.247,13.083,
    250,8.429,8.908,9.516,9.898,11.329,13.175,275,8.489,8.97,9.581,9.966,11.399,13.248,
    300,8.54,9.022,9.635,10.02,11.46,13.326,325,8.587,9.07,9.685,10.071,11.517,13.389,
    350,8.633,9.117,9.732,10.118,11.565,13.44,375,8.67,9.157,9.775,10.161,11.613,13.494,
    400,8.706,9.193,9.814,10.202,11.654,13.542,425,8.738,9.224,9.844,10.234,11.692,13.58,
    450,8.771,9.26,9.882,10.272,11.73,13.623,475,8.798,9.288,9.912,10.302,11.761,13.655,
    500,8.828,9.317,9.939,10.33,11.795,13.69,525,8.854,9.344,9.967,10.36,11.827,13.73,
    550,8.878,9.369,9.995,10.386,11.854,13.751,575,8.901,9.391,10.016,10.408,11.878,13.782,
    600,8.923,9.414,10.04,10.431,11.904,13.813,650,8.963,9.455,10.083,10.476,11.949,13.856,
    700,9.001,9.493,10.119,10.511,11.986,13.904,750,9.033,9.524,10.152,10.547,12.026,13.947,
    800,9.063,9.557,10.187,10.58,12.059,13.975,850,9.093,9.587,10.216,10.612,12.096,14.023,
    900,9.119,9.614,10.244,10.64,12.12,14.041,950,9.143,9.638,10.269,10.665,12.149,14.07,
    1000,9.168,9.664,10.295,10.692,12.176,14.105,1100,9.211,9.708,10.339,10.736,12.22,14.15,
    1200,9.246,9.745,10.377,10.775,12.263,14.197,1300,9.283,9.781,10.415,10.812,12.304,14.235,
    1400,9.313,9.812,10.446,10.845,12.34,14.271,1500,9.347,9.846,10.481,10.88,12.374,14.312,
    1600,9.372,9.871,10.506,10.904,12.396,14.339,2000,9.464,9.965,10.603,11.002,12.5,14.443,
    2500,9.551,10.052,10.69,11.089,12.591,14.54,3000,9.618,10.121,10.76,11.161,12.664,14.619,
    3500,9.675,10.178,10.818,11.219,12.727,14.683,4000,9.727,10.229,10.869,11.271,12.779,14.734,
    4500,9.766,10.269,10.911,11.313,12.82,14.777,5000,9.803,10.307,10.948,11.349,12.859,14.817,
    7500,9.938,10.442,11.085,11.487,12.997,14.959,10000,10.031,10.537,11.18,11.584,13.095,15.063,
    15000,10.152,10.658,11.302,11.707,13.221,15.186,20000,10.236,10.743,11.388,11.791,13.305,15.271,
    50000,10.48,10.988,11.634,12.039,13.556,15.523),
    nrow=108, ncol=7, byrow=TRUE,
    dimnames=list(NULL, c('n','clevel90','clevel92','clevel94','clevel95',
                          'clevel97.5','clevel99')))

  # some constants
  n <- length(x)
  y1 <- x
  Z <- scale(y1)
  
  # compute Tv
  Tv <- c()
  mean.1.2 <- array(data=NA, dim=c(0,2))
  # compute at fixed points, or ...
  if (!is.null(p)) {
    for (v in 1:length(p)) {
      if (!robust) {
        z1 <- mean(Z[1:v], na.rm=TRUE)
        z2 <- mean(Z[(v+1):n], na.rm=TRUE)
      }
      if (robust) {
        z1 <- huber(Z[1:v])$m
        z2 <- huber(Z[(v+1):n])$m
      }
      mean.1.2 <- rbind(mean.1.2, c(z1, z2))
      Tv[v] <- (v*(z1^2)) + ((n-v)*(z2^2))
    }
  } else {
  # ... iterate over the time series
    for (v in 1:(n-1)) {
      if (!robust) {
        z1 <- mean(Z[1:v], na.rm=TRUE)
        z2 <- mean(Z[(v+1):n], na.rm=TRUE)
      }
      if (robust) {
        if (v == 1) {
          z1 <- Z[1]
        } else {
          z1 <- huber(Z[1:v])$m
        }
        if (v == {n-1}) {
          z2 <- Z[n]
        } else {
          z2 <- huber(Z[(v+1):n])$m
        }
      }
      mean.1.2 <- rbind(mean.1.2, c(z1, z2))
      Tv[v] <- (v*(z1^2)) + ((n-v)*(z2^2))
    }
  }
  
  # test Statistic (omit tail ends +-excl elements)
  if (!is.null(p)) {
    T0 <- max(Tv)
    T0x <- p[Tv == T0]
  } else {
    T0 <- max(Tv[excl:(length(Tv)-excl)])
    T0x <- which(Tv == T0)
  }
  T0 <- T0 / cl.factor
  Tc <- as.numeric(cval[max(which(cval[,'n'] <= n)),
                        paste0('clevel', clevel)])
  Tc <- Tc
  if (T0 > Tc) {
    significant <- TRUE
  } else {
    significant <- FALSE
  }
  if (!robust) {
    mean.1 <- mean.1.2[T0x,1] * sd(y1, na.rm=TRUE) + mean(y1, na.rm=TRUE)
    mean.2 <- mean.1.2[T0x,2] * sd(y1, na.rm=TRUE) + mean(y1, na.rm=TRUE)
  } else {
    mean.1 <- mean.1.2[T0x,1] * huber(y1)$s + huber(y1)$m
    mean.2 <- mean.1.2[T0x,2] * huber(y1)$s + huber(y1)$m
  }
  
  res <- list(x=x, Tv=Tv, Tc=Tc, T0=T0, T0x=T0x,
              left_mean=mean.1, right_mean=mean.2,
              clevel=clevel, clf=cl.factor,
              sig=significant)
  return(res)
}




#' Changepoint detection in daily time series.
#' 
#' @description
#' Run the Standard Normal Homogeneity Test to detect a change point in a
#' daily time series, using a one-year moving filter to remove seasonality.
#' 
#' @param can a vector containing the candidate data series.
#' @param aux a matrix containing the auxiliary series.
#' @param tms optional, a vector of Date.
#' @param max.loops an integer indicating the maximum number of iterations of
#' the SNHT; defaults to length(can)-1.
#' @param method see `snht.Q()`.
#' @param min.cor see `snht.Q()`.
#' @param ref.period see `snht.Q()`.
#' @param n.min see `snht.Q()`.
#' @param n.max see `snht.Q()`.
#' @param clevel see `snht()`.
#' @param excl see `snht()`.
#' @param cl.factor a correction factor for the SNHT significance; defaults to
#' 365.
#' @param robust logical, indicating wether or not a robust method should be
#' used for computing the mean and standard deviations; defaults to FALSE.
#' 
#' @param clevel a numeric value representing the confidence level of the
#' test. Must be one of c(90, 92, 94, 95, 97.5, 99). Defaults to 95.
#' @param excl numeric, exclude the first and last elements from the analysis.
#' Defaults to 5.
#' 
#' @details The function computes a backward 365-days moving average time,
#' and then applies the Standard Normal Homogeneity Test recursively until
#' no change points are detected.
#' 
#' @return a list containing the corrected data series, the number of
#' breakpoints, the location of the breakpoints, and a list with the results
#' of each iteration of the SNHT test (see function `snht()`).
#' 
#' @examples
# another example
#' set.seed(123)
#' x <- rnorm(365*20)
#' x[2001:5000] <- x[2001:5000] + 0.4
#' x[5001:6000] <- x[5001:6000] - 0.6
#' 
#' y <- snht.daily(x)
#' y$br_pts
#' snht.plot(y$res[[1]])
#' 
#' @export
snht.daily <- function(can, aux, tms=NULL,
                       max.loops=length(can)-1,
                       method=c('diff', 'ratio'), min.cor=0,
                       ref.period=1:length(can), n.min=1,
                       n.max=Inf, clevel=95, excl=0, cl.factor=365,
                       corr=c('annual', 'monthly', 'daily'),
                       robust=FALSE) {

  if (robust) require(MASS)
  
  # main objects: main series, break points and significance
  x <- can
  br_pts <- NULL
  res <- NULL
  
  # time
  if (corr != 'annual') {
    if (is.null(tmes)) {
      stop('`tms` (a Date series of equal length than `can`) must be provided if `corr` is other than `annual`.')
    }
    mons <- as.character(tmes, format='%m')
    days <- as.character(tmes, format='%d')
  }
  
  # compute cumulative auxiliary series
  aux_acu <- apply(aux, 2, function(x) rowSums(embed(x, 365)) / 365)
  
  # # loop across provided change point candidates, or ...
  # if (!is.null(p)) {
  #   for (i in 1:length(p)) {
  #     
  #   }
  # } else {
  # ... loop the whole series until no significant break point is found
    for (i in 1:max.loops) {
      
      # compute cumulative candidate series
      x_acu <- rowSums(embed(x, 365)) / 365
      
      # compute test series
      Q <- snht.Q(x_acu, aux_acu, method, min.cor, ref.period, n.min, n.max,
                  robust)
      
      # SNHT test
      res[[i]] <- snht(Q$Q, clevel, excl, cl.factor, robust)
      
      # stop if no significant break point is found
      if (res[[i]]$sig == FALSE) {
        break()
      }
      
      # determine break point in the original series, left and right portions
      bp <- res[[i]]$T0x + 364
      n <- length(x)
      left <- 1:bp
      right <- (bp+1):n
      leftQ <- 1:res[[i]]$T0x
      rightQ <- (res[[i]]$T0x+1):n
      
      # stop if repeated break point
      if (bp %in% br_pts) {
        break()
      }
    }
    
    # compute correction coefficients
    QQ <- snht.Q(x, aux, method, min.cor, ref.period, n.min, n.max)
    if (method == 'diff') {
      if (corr == 'annual') {
        if (!robust) {
          lambda <- mean(QQ$Q[rightQ]) - mean(QQ$Q[leftQ])
        } else {
          lambda <- huber(QQ$Q[rightQ])$m - huber(QQ$Q[leftQ])$m
        }
        #lambda <- mean(x[right]) - mean(x[left])
        lambda <- rep(lambda, n)
      }
      if (corr == 'monthly' | corr == 'daily') {
        if (!robust) {
          lambda <- aggregate(QQ$Q[rightQ], list(mons[rightQ+365]), mean)[,2] -
            aggregate(QQ$Q[leftQ], list(mons[leftQ+365]), mean)[,2]
        } else {
          lambda <- aggregate(QQ$Q[rightQ], list(mons[rightQ+365]),
                      function(x) huber(x)$m)[,2] -
            aggregate(QQ$Q[leftQ], list(mons[leftQ+365]),
                      function(x) huber(x)$m)[,2]
        }
        #lambda <- aggregate(x[right], list(mons[right]), mean)[,2] -
        # aggregate(x[left], list(mons[left]), mean)[,2]
        lambda <- lambda[as.numeric(mons)]
      }
      if (corr == 'daily') {
        w <- which(days == '15')
        lambda <- ksmooth(w, lambda[w], 'normal', x.points=1:n, bandwidth=28)$y
      }
      x[left] <- x[left] + lambda[left]
    }
    if (method == 'ratio') {
      if (corr == 'annual') {
        if (!robust) {
          lambda <- mean(QQ$Q[rightQ]) / mean(QQ$Q[leftQ])
        } else {
          lambda <- huber(QQ$Q[rightQ])$m / huber(QQ$Q[leftQ])$m
        }
        #lambda <- mean(x[right]) / mean(x[left])
        lambda <- rep(lambda, n)
      }
      if (corr == 'monthly' | corr == 'daily') {
        if (!robust) {
          lambda <- aggregate(QQ$Q[rightQ], list(mons[rightQ+365]), mean)[,2] /
          aggregate(QQ$Q[leftQ], list(mons[leftQ+365]), mean)[,2]
        } else {
          lambda <- aggregate(QQ$Q[rightQ], list(mons[rightQ+365]),
                              function(x) huber(x)$m)[,2] /
                    aggregate(QQ$Q[leftQ], list(mons[leftQ+365]),
                              function(x) huber(x)$m)[,2]
        }
        #lambda <- aggregate(x[right], list(mons[right]), mean)[,2] /
        # aggregate(x[left], list(mons[left]), mean)[,2]
        lambda <- lambda[as.numeric(mons)]
      }
      if (corr == 'daily') {
        w <- which(days == '15')
        lambda <- ksmooth(w, lambda[w], 'normal', x.points=1:n, bandwidth=28)$y
      }
      x[left] <- x[left] * lambda[left]
    }
    
    # store break point
    br_pts[i] <- bp

  } # main loop
  
  return(list(x=x, n=length(br_pts), br_pts=br_pts, res=res))
}





#' Create a plot to check the results of the SNHT.
#' 
#' @param res the result of a call to `snht()`.
#' 
#' @return Two plots: a plot with the time series being checked with the means
#' corresponding to the two sections separated by the change point; and a time
#' series of the Tv statistic with indication of the confidence levels and
#' indication of the change point (if significant).
snht.plot <- function(res) {
  
  mean_left <- mean_right <- res$x * NA
  n <- length(res$x)
  mean_left[1:res$T0x] <- res$left_mean
  mean_right[res$T0x:n] <- res$right_mean
  par(mfrow=c(2,1), mar=c(2,4,1,1) + 0.1)
  
  # series
  plot(res$x, type="l", xlab=NA, ylab='Q')
  if (res$sig) {
    lines(mean_left)
    lines(mean_right)
  }
  
  # statistic
  plot(res$Tv / res$clf, type="l", xlab=NA, ylab='Tv',
       ylim=c(0, max(c(res$Tv / res$clf, res$Tc)))
       )
  if (res$sig) {
    abline(v=res$T0x, col="Red")
  }
  abline(h=res$Tc, lty="dashed")
}




#' Create a plot to check the results of the SNHT applied to daily series.
#' 
#' @param ori the original time series.
#' @param new the result of call to `snht.daily()`.
#' @param tms a vector of dates.
#' 
#' @return ´three plots: a time series of original and corrected data; a time
#' series of the original and corrected data with a one-year filter; the result
#' of all the SNHT tests applied recursively after no more change points are
#' found.
snht.daily.plot <- function(ori, new, tms) {
  
  # daily series, original and corrected
  par(mfrow=c(2,1), mar=c(2,4,1,1) + 0.1)
  plot(ori ~ tms, type='p', pch='.', ylab='Tmn (ºC)')
  abline(v=tms[new$br_pts], col='red')
  grid(nx=0)
  plot(new$x ~ tms, type='p', pch='.', ylab='Tmn_crr (ºC)')
  abline(v=tms[new$br_pts], col='red')
  grid(nx=0)
  
  # moving window, original and corrected
  par(mar=c(3,4,2,1) + 0.1)
  can_acu <- rowSums(embed(ori, 365)) / 365
  plot(can_acu ~ tms[-c(1:364)], type='p', pch='.',
       xlab=NA, ylab='Tmn 365_mean, ori. (ºC)', main=nam)
  abline(v=tms[new$br_pts], col='red')
  grid(nx=0)
  can_new_acu <- rowSums(embed(new$x, 365)) / 365
  par(mar=c(3,4,2,1) + 0.1)
  plot(can_new_acu ~ tms[-c(1:364)], type='p', pch='.',
       xlab=NA, ylab='Tmn_new 365_mean, corr. (ºC)', main=nam)
  abline(v=tms[new$br_pts], col='red')
  grid(nx=0)
  
  # SNHT tests
  nr <- length(new$res)
  if (nr > 0) {
    for (i in 1:nr) {
      snht.plot(new$res[[i]])
    }
  }
}



