# Verifies if needed packages are installed
pkgs_att <- c("tidyverse",
              "GA",
              "minpack.lm",
              "parallel")
new.pkgs <- pkgs_att[!(pkgs_att %in% installed.packages()[,"Package"])]

# Installs new packages
if (length(new.pkgs)) install.packages(new.pkgs)

# Loads and attaches all needed packages
for (s in pkgs_att) {
  if (!library(s, character.only = TRUE, logical.return = TRUE)) 
    stop("Error loading packages.")
}

#' generic SVI function with parameters: 
#' a,b, rho, m, sigma
svi_fun <- function(param, k){
  a <- param[1]
  b <- param[2]
  rho <- param[3]
  m <- param[4]
  sigma <- param[5]
  ans <- a + b*(rho*(k - m) + sqrt((k - m)^2 + sigma^2))
  return(ans)
}

#' Estimate svi parameers with GA algorithm
#' 
#' This require ga package, GA gives robust result, though 
#' it is relative slow. Alternative, use constOptim
#' It takes two-step optimization:
#' 1. obtain inner parameters: a,b, rho
#' 2. estimate m and sigma, H using GA
#' @param k log(K/F) -a vector
#' @param iv implied volatility - a vector
#' @param tau maturity (scalar)
#' 
fit_ga_svi <- function(k, iv, tau){
  # Idea: weights can be passed to this function to weight the errors in 
  # ofitfun
  # total iv
  w <- iv*iv*tau
  
  # optimization function: x (m, sigma)
  ofitfun <- function(x) { # m, sigma
    y <- (k - x[1])/x[2] # m, sigma
    # find inner parameters: a, b, rho
    par <- fit_inner(y, w, x[2]) # a, b, rho
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
  par <- fit_inner(y, w, opar[2])
  # parameters: a, b, rho,m, sigma
  parameters <- c(a = par[1], b = par[2], rho = par[3],
                  m = opar[1], sigma = opar[2])
  ans <- tibble::tibble(par = list(parameters))
  return(ans)
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
#' @param y (k-m)/sigma
#' @param w  total iv, iv*sqrt(tau)
#' @param sig sigma
#' @param H fractional exponent
#' @return a,b,rho parameters
#' @export
fit_inner <- function(y, w, sig){
  
  # Idea: one can include a parameter method wich accepts "quadprog"
  # or "Nelder-Mead". 
  # If quadprog, then compute matrices D and A, and vector d from 
  # parameters y, w and sigma.
  # Then solves the quadratic problem with solve.QP
  
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
  
  # initial guess
  x0 <- c(max(w)/2, -sig, 2*sig)
  
  # constrained optimization, Nelder-Mead by default
  op <- constrOptim(x0,fitfun, NULL, ui, ci)
  parameters <- op$par
  
  a <- parameters[1]
  b <- parameters[3]/sig
  rho <- parameters[2]/parameters[3]
  return(c(a, b, rho))
}


#' Estimated svi parameters with a given H
#'
#' Estimate modified svi parameters with a given H.
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @param H exponent H for the modified SVI function
#' @export
fit_svih <- function(k, iv, tau) {
  # two-step optimization:
  # 1. obtain inner parametplers: a,b, rho
  # 2. estimate m and sigma
  w <- iv*iv*tau
  ofitfun <- function(x) { # m, sigma
    y <- (k - x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2]) # a, b, rho
    param <- c(par[1],par[2], par[3], x[1], x[2]) #a, b, m, rho, sigma
    raw <- svi_fun(param, k)
    return(sum((raw - w)^2))
  }
  
  # intial guess of m, sigma
  x0 <- c( 0.0, 0.02)
  
  ui <- rbind(c(1,0),c(-1,0), c(0,1), c(0,-1))
  ci <- c( min(k), -max(k), 0.0001, -5)
  
  op <- constrOptim(x0, ofitfun, NULL, ui, ci)
  opar <- op$par
  y <- (k - opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2])
  # parameters: a, b, m, rho, sigma
  parameters <- c(a = par[1], b = par[2], rho = par[3],
                  m = opar[1], sigma = opar[2])
  ans <- tibble::tibble(par = list(parameters))
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

#' direct procedure to fit an SVI slice. Uses a constrained Levenberg-Marquardt
#' algorithm on the 5-dim problem.
#' Depends on minpack.lm package
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @export
fit_direct_svi <- function(k, iv, tau) {
  w <- iv * iv * tau # Total implied variance observed on market
  restarts <- 10 # Number of random restarts
  
  # Residual (error) function. It returns a vector of residues. The optmization 
  # will minimize the sum square of these residues.
  residFun <- function(par, x, observed){
    wi <- svi_fun(par, x)
    return(wi - observed)
  }
  
  # Initial guesses
  ai <- rep(min(w), restarts)
  bi <- runif(restarts, 0, 0.05)
  rhoi <- rep(ifelse(w[1] - w[length(w)] < 0, -0.5, 0.3), restarts)
  mi <- runif(restarts, 2 * min(k), 2 * max(k))
  sigmai <- runif(restarts, 0, 1)
  
  # Holding list of results
  results <- vector("list", length = restarts)
  
  # Proceed to optimizations. # of optmizations = restarts
  for (i in 1:restarts) {
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
    mutate(par = map(value, ~.x$par),
           sumse = map_dbl(value, ~.x$deviance),
           info = map_int(value, ~.x$info)) %>% 
    filter(sumse == min(sumse))
  return(ans)
}