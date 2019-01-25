# Verifies if needed packages are installed
pkgs_att <- c("tidyverse",
              "GA",
              "minpack.lm",
              "quadprog",
              "parallel")
new.pkgs <- pkgs_att[!(pkgs_att %in% installed.packages()[,"Package"])]

# Installs new packages
if (length(new.pkgs)) install.packages(new.pkgs)

# Loads and attaches all needed packages
for (s in pkgs_att) {
  if (!library(s, character.only = TRUE, logical.return = TRUE)) 
    stop("Error loading packages.")
}


#' Total variance for a given set of SVI RAW parameters.
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return Total variance evaluated at k points
#' @export
#'
#' @examples
svi_fun <- function(par, k){
  a <- par[1]
  b <- par[2]
  rho <- par[3]
  m <- par[4]
  sigma <- par[5]
  ans <- a + b*(rho*(k - m) + sqrt((k - m)^2 + sigma^2))
  return(ans)
}

#' First derivative of RAW SVI with respect to moneyness
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return First derivative evaluated at k points
#' @export
#'
#' @examples
diff_svi <- function(par, k) {
  a <- par[1]
  b <- par[2]
  rho <- par[3]
  m <- par[4]
  sigma <- par[5]
  ans <- b * (rho + (k - m) / (sqrt((k - m)^2 + sigma^2)))  
  return(ans)
}

#' Second derivative of RAW SVI with respect to moneyness
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return Second derivative evaluated at k points
#' @export
#'
#' @examples
diff2_svi <- function(par, k) {
  a <- par[1]
  b <- par[2]
  rho <- par[3]
  m <- par[4]
  sigma <- par[5]
  disc <- (k - m)^2 + sigma^2
  ans <- (b * sigma^2) / ((disc)^(3 / 2))
  return(ans)
}

#' Computes g(k) function. Auxiliary to retrieve density and essential to test for butterfly arbitrage.
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return Function g(k) evaluated at k points
#'
#' @examples
gfun <- function(par, k) {
  w <- svi_fun(par, k)
  w1 <- diff_svi(par, k)
  w2 <- diff2_svi(par, k)
  
  ans <- (1 - 0.5*(k*w1/w))^2 - (0.25*w1^2)*(0.25 + w^-1) + 0.5*w2
  return(ans)
}


#' Auxiliary function to compute d1 from BSM model, given a set of SVI RAW parameters and a vector of moneynesses.
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return
#'
#' @examples
d1_svi <- function(par, k) {
  v <- sqrt(svi_fun(par, k))
  return(-k/v + v/2)
}

#' Auxiliary function to compute d2 from BSM model, given a set of SVI RAW parameters and a vector of moneynesses.
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return
#'
#' @examples
d2_svi <- function(par, k) {
  v <- sqrt(svi_fun(par, k))
  return(-k/v - v/2)
}

#' Probability density implied by an SVI.
#'
#' @param par Set of raw parameters, (a, b, rho, m, sigma)
#' @param k Moneyness points to evaluate
#'
#' @return Implied risk neutral probability density from an SVI
#' @export
#'
#' @examples
svi_density <- function(par, k) {
  g <- gfun(par, k)
  w <- svi_fun(par, k)
  dtwo <- d2_svi(par, k)
  
  ans <- (g/sqrt(2*pi*w))*exp(-0.5*dtwo^2)
  return(ans)
}

#' Computes RMSE (root mean squared error) of a raw parameterization.
#' @param df Data frame
#' @param param Raw parameters in a vector (a, b, rho, m, sigma)
#' @param w Observed total variances
#' 
#' @return Root mean squared error
#' @export
#' 
#' @examples 
rmse <- function(par, k, w) {
  return(sqrt(mean((svi_fun(par, k) - w)^2)))
}

#' Inner optimization for parameter: c,d,a
#'
#' optimization of the inner parameters
#' let y(k) = (k - m)/sig
#' svi w = a + b*sig*(rho*y(k) + sqrt(y(x)^2+1))
#'        = a + dy(x) + cz(x)
#' d = roh*b*sig, c = b*sig, z = sqrt(y(x)^2+1)
#' Constrains:
#'     0 <= c <= 4sig
#'     |d| <= c
#'     |d| <= 4sig -c
#'     0 <= a <= max_i(w_i)
#'     
#' @param y (k-m)/sigma
#' @param w  Total iv, iv*sqrt(tau)
#' @param sig Sigma
#' @param method One of Nelder-Mead (default) or quadprog  
#' 
#' @return a, b, rho parameters
#' 
fit_inner <- function(y, w, sig, method = c("Nelder-Mead", "quadprog")){
  met = method[1] # If not set, then takes default method

  # constrain optimization for (a, d, c)
  fitfun <- function(x) {
    yhat <- (x[1] + x[2]*y + x[3]*sqrt(y^2 + 1))
    sum((yhat - w)^2)
  }
  
  # constrains (see above)
  ui <- rbind(c(0,0,1),c(0,0,-1),        # constrains on c
              c(0,-1,1), c(0, 1, 1),     # constrains on |d| <= c
              c(0, 1, -1), c(0, -1, -1), # constrains on |d| <= 4*sig - c
              c(1, 0, 0), c(-1, 0, 0))   # constrains on a
  ci <- c(0, -4*sig, 
          0, 0, 
          -4*sig, -4*sig, 
          0, -max(w))
  
  # initial guess for NM
  x0 <- c(max(w)/2, -0.9*sig, sig)
  
  # If quadprog method was set
  if (met == "quadprog") {
    # Min (-d^T b + 0.5 b^T D b) s.t A^T b >= b0
    # See quadprog documentation
    # Input values for D matrix
    n <- length(y)
    sum_y <- sum(y)
    sum_z <- sum(sqrt(y^2 + 1))
    sum_y2 <- sum(y^2)
    sum_yz <- sum(y * sqrt(y^2 + 1))
    sum_z2 <- sum(y^2 + 1)
    # Input values for d vector
    sum_w <- sum(w)
    sum_wy <- sum(w * y)
    sum_wz <- sum(w * (y^2 + 1))
    
    # Matrix D
    D <- matrix(c(n, sum_y, sum_z,
                  sum_y, sum_y2, sum_yz,
                  sum_z, sum_yz, sum_z2),
                nrow = 3,
                byrow = TRUE)
    # Vector d
    d <- c(sum_w, sum_wy, sum_wz)
    
    # b0 is equal to ci
    # A is equal to ui transposed
    op <- solve.QP(D, d, t(ui), ci)
    parameters <- op$solution
    
  } else {# Otherwise, constrOptim runs Nelder-Mead
    
    # constrained optimization, Nelder-Mead by default
    op <- constrOptim(x0, fitfun, grad = NULL, ui = ui, ci = ci)
    parameters <- op$par
  }
  
  a <- parameters[1]
  b <- parameters[3]/sig
  rho <- parameters[2]/parameters[3]
  return(c(a, b, rho))
}


#' Estimate SVI parameters with GA algorithm
#' 
#' This require ga package, GA gives robust result, though 
#' it is relative slow. Alternative, use constOptim
#' It takes two-step optimization:
#' 1. obtain inner parameters: a,b, rho
#' 2. estimate m and sigma, H using GA
#' @param k Forward log-moneyness
#' @param w Total implied variance
#' @param inner Method for inner optimization (i.e Reduced problem)
#' 
#' @return A Data Frame containing the following elements: 
#'     \item{method} String indicating method of optimization. "GA".
#'     \item{par}  List column of named vector of parameters (a, b, rho, m, sigma).
#'     \item{objective} Objective function. The inverse of squared errors mean.
#'     \item{convergence}  An integer code about reasons for optmization ending. 0 means convergence while 1 you should debug it.
#' @export
#' 
#' @examples 
svi_fit_ga <- function(k, w, inner = c("Nelder-Mead", "quadprog")){
  met = inner[1] # Method for inner optmization
  # Idea: weights can be passed to this function to weight the errors in 
  # ofitfun

  # optimization function: x (m, sigma)
  ofitfun <- function(x) { # m, sigma
    y <- (k - x[1])/x[2] # m, sigma
    # find inner parameters: a, b, rho
    par <- fit_inner(y, w, x[2], met) # a, b, rho
    param <- c(par[1],par[2], par[3], x[1], x[2]) #a, b, rho, m, sigma
    raw <- svi_fun(param, k)
    value <- 1/mean((raw - w)^2)
    return(value)
  }
  
  # genetic Algorith, set boxes for m, sigma
  op <- GA::ga("real-valued", fitness = ofitfun,
               lower = c(min(k), 0.00001),
               upper = c(max(k), 10),
               maxiter = 80, run = 80, optim = TRUE,
               monitor = FALSE)
  # possible multiple solutions, only select the first one
  opar <- op@solution[1,]
  
  y <- (k - opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2], met)
  # parameters: a, b, rho,m, sigma
  parameters <- c(par[1], par[2], par[3], opar[1], opar[2])
  names(parameters) <- c("a", "b", "rho", "m", "sigma")
  ans <- tibble::tibble(method = "GA",
                        par = list(parameters),
                        objective = op@fitnessValue,
                        convergence = 0) # GA always converge?
  return(ans)
}



#' Estimate SVI parameters with through quasi-explicit algorithm.
#'
#' @param k Forward log-moneyness
#' @param w Total implied variance
#' @param n_rst Number of complete problem restarts 
#' @param inner Method for inner optimization (i.e Reduced problem)
#' 
#' @return A Data Frame containing the following elements: 
#'     \item{method} String indicating method of optimization. "Quasi".
#'     \item{par}  List column of named vector of parameters (a, b, rho, m, sigma).
#'     \item{objective} The objective function. Sum of squared residuals.
#'     \item{convergence}  An integer code about reasons for optmization ending. 0 means convergence while 1 you should debug it. Note: it is complete problem convergence (m, sigma).
#'
#' @export
svi_fit_quasi <- function(k, w, n_rst = 50, inner = c("Nelder-Mead", "quadprog")) {
  met = inner[1] # Method for inner optmization
  # two-step optimization:
  # 1. obtain inner parametplers: a, b, rho
  # 2. estimate m and sigma
  ofitfun <- function(x) { # m, sigma
    y <- (k - x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2], met) # a, b, rho
    param <- c(par[1],par[2], par[3], x[1], x[2]) #a, b, rho, m, sigma
    raw <- svi_fun(param, k)
    return(sum((raw - w)^2))
  }
  
  # intial guess of m, sigma
  mi <- runif(n_rst, 2*min(k), 2*max(k)) 
  sigmai <- runif(n_rst, 0.00011, 1)
  
  # boundaries for m, sigma
  ui <- rbind(c(1,0),
              c(-1,0), 
              c(0,1), 
              c(0,-1))
  ci <- c( -999, 
           -999, 
           0.0001, 
           -999)
  # Holding list of results
  results <- vector("list", length = n_rst)
  
  # Proceed to optimizations. # of optmizations = restarts
  for (i in 1:n_rst) {
    op <- constrOptim(c(mi[i], sigmai[i]), ofitfun, NULL, ui, ci)
    opar <- op$par
    y <- (k - opar[1])/opar[2]
    par <- fit_inner(y, w, opar[2], met)
    # parameters: a, b, m, rho, sigma
    parameters <- c(par[1], par[2], par[3], opar[1], opar[2])
    names(parameters) <- c("a", "b", "rho", "m", "sigma")
    results[[i]] <- tibble::tibble(par = list(parameters),
                                   objective = op$value,
                                   convergence = op$convergence) 
  }
  
  # Select string for method
  if (met == "quadprog") 
    method <- "QuasiPQ"
  else 
    method <- "QuasiNM"
  
  ans <- do.call(rbind, results) %>% 
    filter(objective == min(objective)) %>% 
    mutate(method = method) %>% 
    select(method, par, objective, convergence)
  return(ans)
}

#' parallel computing a surface with a given H
#' @param chain surface table contains date, tau, k, iv
#' @param H H exponent
#' 
par_fit_svi <- function(chain, H){

  # use parallel computing
  # detect number of cores
  no_cores <- detectCores() - 1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr)})
  clusterExport(cl, c('fit_svih','fit_inner','svi_fun'))
  # fit function
  fit_curve <- function(c){
    param <- fit_svih(c$k,c$iv, c$tau[1])
    return(param)
  }
  
  # split chains to each dte
  lchain <- split(chain, chain[,'tau'])

  # use paralell computing
  paras <- parLapply(cl, lchain,fit_curve)
  # convert list to data.frame
  surface <- do.call(rbind, paras)
  # combine back date and dte
  ch <- chain %>% distinct(date, tau)
  surf <- left_join(surface, ch, by = c('tau'))
  # stop parallel
  stopCluster(cl)
  return(surf)
}

#' Direct procedure to fit an SVI slice. Uses a constrained Levenberg-Marquardt
#'     algorithm on the 5-dim problem. The algorithm is restarted a certain
#'     number of times with random initial guesses. The best round (i.e lower 
#'     objective function) is returned. Param df must be a data frame with at 
#'     least these columns: Spot Price (S), Risk-free rate (r), Forward 
#'     log-moneyness (k), implied volatility (iv) and time to expiration (tau).
#'     Depends on minpack.lm package
#' @param k Forward log-moneyness
#' @param w Total implied variance
#' @param n_rst Number of random restarts
#' 
#' @return A Data Frame containing the following elements: 
#'     \item{method} String indicating method of optimization. "Direct".
#'     \item{par}  List column of named vector of parameters (a, b, rho, m, sigma).
#'     \item{objective} Sum of squared residuals. The objective function.
#'     \item{convergence}  An integer code about reasons for optmization ending. 0 means convergence while 1 you should debug it.
#' @export
#' 
#' @examples 
svi_fit_direct <- function(k, w, n_rst = 10) {

  # Residual (error) function. It returns a vector of residues. The optmization 
  # will minimize the sum square of these residues.
  residFun <- function(par, x, observed){
    wi <- svi_fun(par, x)
    return(wi - observed)
  }
  
  # Initial guesses
  ai <- rep(min(w), n_rst)
  bi <- runif(n_rst, 0, 0.05)
  rhoi <- rep(ifelse(w[1] - w[length(w)] < 0, -0.5, 0.3), n_rst)
  mi <- runif(n_rst, 2 * min(k), 2 * max(k))
  sigmai <- runif(n_rst, 0, 1)
  
  # Holding list of results
  results <- vector("list", length = n_rst)
  
  # Proceed to optimizations. # of optmizations = restarts
  for (i in 1:n_rst) {
    init_pars <- c(ai[i], bi[i], rhoi[i], mi[i], sigmai[i])
    
    results[[i]] <- minpack.lm::nls.lm(init_pars, 
                                     lower = c(0, 0, -1, -Inf, 0),
                                     upper = c(max(w), Inf, 1, Inf, Inf),
                                     fn = residFun,
                                     x = k,
                                     observed = w,
                                     control = list(maxiter = 100,
                                                    ftol = .Machine$double.eps,
                                                    ptol = .Machine$double.eps))
  } # end of for loop
  
  ans <- enframe(results) %>% 
    mutate(method = "Direct",
           par = map(value, ~.x$par),
           objective = map_dbl(value, ~.x$deviance),
           convergence = map_int(value, ~.x$info)) %>% 
    filter(objective == min(objective)) %>% 
    mutate(convergence = if_else(convergence %in% 1:4, 0, 1)) %>% 
    select(method, par, objective, convergence)
  names(ans$par[[1]]) <- c("a", "b", "rho", "m", "sigma")
  return(ans)
}

#' Check data frame passed as argument to several functions
#' Param df must be a data frame with at least these columns:
#' Spot Price (S), Risk-free rate (r), Forward log-moneyness (k), 
#' implied volatility (iv) and time to expiration (tau)
#' @param df Data frame
#' 
check_df_argument <- function(df) {
  if (nrow(df) < 5) {
    stop("Data frame must have at least 5 observations (rows).")}
  
  col <- c("S", "r", "k", "iv", "tau")
  if (!prod(col %in% colnames(df))) {
    stop("Data frame must have at least 5 columns as in description.")
  }
}
