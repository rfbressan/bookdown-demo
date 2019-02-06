#' ---
#' title: "Functions to fit a SSVI"
#' author: "Rafael F. Bressan"
#' --- 
#' 
#' # Phi function
#'     
#' @param par Complete list of parameters for the SSVI function in that
#'     order, (rho, gamma, [eta]). For "heston" type only gamma will be
#'     selected, for "powerlaw", gamma and eta
#' @param theta ATM implied total variance
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#'     
#' @return Values for the chosen type of phi function
#' 
#' @details Auxiliary to compute a SSVI. Type can be either "heston" (
#'     default) or "powerlaw", in which case par must have an addtional third 
#'     value for eta.
#+ echo = TRUE
phi_fun <- function(par, theta, phitype){
  # modified powerlaw as given in Gatheral and Jacquier (2013) eq. 4.5, pg. 17
  if (phitype == "powerlaw") {
    gamma <- par[2]
    eta = par[3]
    return(eta/(theta^gamma * (1 + theta)^(1 - gamma)))
  }
  else {
    gamma = par[2]
    return(1. / (gamma*theta)*(1. - (1. - exp(-gamma*theta))/(gamma*theta)))
  }
}

#' # ATM implied total variance, interpolated through a spline.
#'
#' @param k A [KTx1] vector of forward log-moneynes
#' @param t A [KTx1] vector of times to maturity
#' @param iv A [KTx1] vector of implied volatilities
#'
#' @return A vector of ATM implied (interpolated) total variances of the same 
#'     length as k, t and iv.
#' @export
#' 
#' @examples
#+ echo = TRUE
theta_vec <- function(k, t, iv) {
  if (min(t) <= 0.0) stop("Parameter t must be greater than zero.")
  w <- iv^2 * t # Total variance
  df <- data.frame(k, t, w)
  groups <- split(df, df$t)
  # sp_fun: spline interpolation of smiles and returns ATM total variance
  sp_fun <- function(df) {
    spline(df$k, df$w, method = "natural", 
           xout = rep(0, length(df$k)))$y
  }
  # sapply sp_fun to groups to have a vector of ATM variances
  atm <- unlist(lapply(groups, sp_fun))
  names(atm) <- NULL  # No need to keep names 
  return(atm)
}

#' # Evalaute ATM total variance
#'
#' @param theta_vec Vector of ATM total variances
#' @param t_vec Corresponding vector of time to maturity
#' @param t_eval Vector of times to evaluate ATM total variance
#'
#' @return A vector of ATM total variances corresponding to t_eval. A spline 
#'     interpolation is performed over theta_vec and t_vec and then evaluated
#'     at t_eval.
#' @export
#'
#' @examples
#+ echo = TRUE
theta_fun <- function(theta_vec, t_vec, t_eval) {
  if (min(theta_vec) <= 0 | min(t_vec) <= 0 | min(t_eval) <= 0)
    stop("All values in theta_vec, t_vec or t_eval must be greater than zero.")
  
  return(spline(unique(t_vec), unique(theta_vec), method = "natural", 
                xout = t_eval)$y)
}

#' # SSVI function
#'
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k A [KTx1] vector of forward log-moneynes
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#'
#' @return Values of SSVI function for the pair (theta, k)
#' @export
#'
#' @examples
#+ echo = TRUE
ssvi_fun <- function(par, theta, k, phitype = "powerlaw") {
  p <- phi_fun(par, theta, phitype)
  return(0.5*theta*(1. + par[1]*p*k + sqrt((p*k + par[1])^2 + 1 - par[1]^2)))
}


#' # SSVI function 2. Constant sigma.
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param sigma Constant implied variance, such that $\theta := \sigma^2 t$
#' @param t Time to maturity
#' @param k Forward log-moneyness
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return: Values of SSVI function for the pair (theta, k), considering
#'     that $\theta := \sigma^2 t$.
#+ echo = TRUE
ssvi_fun2 <- function(par, sigma, t, k, phitype = "powerlaw"){
  theta <- sigma^2 * t
  
  p <- phi_fun(par, theta, phitype)
  
  return(0.5*theta*(1. + par[1]*p*k + sqrt((p*k + par[1])^2 + 1 - par[1]^2)))
}


#' # First derivative of SSVI with respect to k
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k A [KTx1] vector of forward log-moneynes
#' @param phitype Text string representing the type of phi function. One of
#' "heston" or "powerlaw".
#' @return: Values of diff(SSVI, k) for the pair (theta, k), considering
#' that $\theta := \sigma^2 t$.
#+ echo = TRUE
ssvi_diff <- function(par, theta, k, phitype){
  p <- phi_fun(par, theta, phitype)
  rho <- par[1]
  pkr <- p*k + rho
  
  return(0.5*theta*p*(rho + pkr/sqrt(pkr^2 + 1 - rho^2)))
}

#' # Second derivative of SSVI with respect to k
#' 
#' @param par Parameters of the function in that order, (rho, sigma, gamma,
#'     [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k A [KTx1] vector of forward log-moneynes
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return Values of diff(SSVI, k, k) for the pair (theta, k), considering
#'     that $\theta := \sigma^2 t$.
#+ echo = TRUE
ssvi_diff2 <- function(par, theta, k, phitype){
  p <- phi_fun(par, theta, phitype)
  rho <- par[1]
  pkr <- p * k + rho
  
  return(0.5*p^2*theta*(1 - rho^2)/(pkr^2 + 1 - rho^2)^(3/2))
}


#' # First derivative of SSVI with respect to t
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param t A [KTx1] vector of times to maturity
#' @param k A [KTx1] vector of forward log-moneynes
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return Values of diff(SSVI, t) for the pair (theta, k).
#+ echo = TRUE
ssvi_difft <- function(par, t, k, theta, phitype){
  # Finite difference method
  eps <- 1e-6
  tp <- theta_fun(theta, t, t + eps)
  tm <- theta_fun(theta, t, t - eps)
  return((ssvi_fun(par, tp, k, phitype) - ssvi_fun(par, tm, k, phitype))/(2*eps))
}


#' # Computes the g(k) function from an SSVI parameters
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k Forward log-moneyness
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return The g(theta, k) function, considering that $\theta := \sigma^2 t$.
#+ echo = TRUE
ssvi_g <- function(par, theta, k, phitype){
  w <- ssvi_fun(par, theta, k, phitype)
  w1 <- ssvi_diff(par, theta, k, phitype)
  w2 <- ssvi_diff2(par, theta, k, phitype)
  
  return((1 - 0.5*(k*w1/w))^2 - (0.25*w1^2)*(w^-1 + 0.25) + 0.5*w2)
}


#' # Auxiliary function to compute d1 from BSM model.
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k Forward log-moneyness
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return Returns d1 from BSM model evaluated at (k, w) points
#+ echo = TRUE
ssvi_d1 <- function(par, theta, k, phitype){
  v <- sqrt(ssvi_fun(par, theta, k, phitype))
  return(-k/v + 0.5 * v)
}


#' # Auxiliary function to compute d2 from BSM model.
#' 
#' @param par Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k Forward log-moneyness
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return Returns d2 from BSM model evaluated at (k, w) points
#+ echo = TRUE
ssvi_d2 <- function(par, theta, k, phitype){
  v <- sqrt(ssvi_fun(par, theta, k, phitype))
  return(-k/v - 0.5 * v)
}


#' # Probability density implied by an SSVI.
#' 
#' @param par: Parameters of the function in that order, (rho, gamma, [eta])
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k: Forward log-moneyness
#' @param phitype: Text string representing the type of phi function. One of
#' "heston" or "powerlaw".
#' @return: Implied risk neutral probability density from a SSVI
#+ echo = TRUE
ssvi_density <- function(par, theta, k, phitype){
  g <- ssvi_g(par, theta, k, phitype)
  w <- ssvi_fun(par, theta, k, phitype)
  dtwo <- ssvi_d2(par, theta, k, phitype)
  
  dens <- (g / sqrt(2 * pi * w)) * exp(-0.5 * dtwo^2)
  return(dens)
}


#' # Local Volatility function through Dupire's equation for a SSVI
#'     
#' @param par Parameters of the function in that order, (rho, sigma, gamma,
#'     [eta])
#' @param t Time to maturity
#' @param theta A [KTx1] vector of ATM total implied variances
#' @param k Forward log-moneyness
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @return Local volatility function evaluated at (k, t) points
#+ echo = TRUE
ssvi_local_vol <- function(par, t, theta, k, phitype){
  return(sqrt(ssvi_difft(par, t, k, theta, phitype) / ssvi_g(par, theta, k, phitype)))
}


#' # Fit a SSVI to market data
#'
#' @param t A [KTx1] vector of times to maturity
#' @param k A [KTx1] vector of forward log-moneyness
#' @param iv A [KTx1] vector of implied volatilities
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#' @param iter Number of iterations to optmize
#'
#' @return Parameters defining the SSVI function, (rho, gamma, [eta])
#' @export
#'
#' @examples
#+ echo = TRUE
fit_ssvi <- function(t, k, iv, phitype = "powerlaw", iter = 50) {
  theta <- theta_vec(k, t, iv)
  wmkt <- iv^2 * t
  
  # Bounds on parameters (rho, gamma, [eta])
  bounds <- ssvi_bounds(phitype)
  # Object function for GA
  ga_fun <- function(gapar) {
    w <- ssvi_fun(gapar, theta, k, phitype)
    obj <- 1 / mean((w - wmkt)^2)
    return(obj)
  }
  # Objective function for LM
  lm_fun <- function(lmpar) {
    w <- ssvi_fun(lmpar, theta, k, phitype)
    # Insert penalty terms to avoid arbitrage
    # add penalty for the parameter constrains
    pen <- sqrt(.Machine$double.xmax)/4  # penalty term
    penalty <- max(ssvi_butterfly_cons(lmpar, phitype), 0.0)*pen
    
    # Remember: nls.lm already tries to minimize the sum square
    ans <- w - wmkt + penalty
    return(ans)
  }  
  # genetic Algorith, set boxes for m, sigma
  # Slightly modify bounds in order to avoid initial guesses AT boundary for LM
  op <- GA::ga("real-valued", fitness = ga_fun,
               lower = bounds$lower + sqrt(.Machine$double.eps),
               upper = bounds$upper - sqrt(.Machine$double.eps),
               maxiter = iter, run = iter, optim = TRUE,
               monitor = FALSE)
  # possible multiple solutions, only select the first one
  gapar <- op@solution[1,]
  
  # Levenber-Marquardt to refine solution
  lmop <- minpack.lm::nls.lm(gapar, 
                             lower = bounds$lower,
                             upper = bounds$upper,
                             fn = lm_fun,
                             control = list(maxiter = iter,
                                            ftol = .Machine$double.eps,
                                            ptol = .Machine$double.eps))
  
  lmpar <- lmop$par
  names(lmpar) <- switch(phitype,
                         powerlaw = c("rho", "gamma", "eta"),
                         heston = c("rho", "gamma")
  )
  return(lmpar)
}


#' # Bounds on SSVI parameters for "powerlaw" and "heston" phi types
#'
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#'
#' @return A list with lower and upper elements
#'
#' @examples
#+ echo = TRUE
ssvi_bounds <- function(phitype) {
  
  # rho, gamma, [eta]
  inf <- 1e6
  switch(phitype,
    powerlaw = list(lower = c(-1., 0.0, 0.0),
                    upper = c(1., 1., inf)),
    heston = list(lower = c(-1., 0.0),
                  upper = c(1., inf))
  )
}


#' # Butterfly constraints on a SSVI
#'
#' @param x Parameters
#' @param phitype Text string representing the type of phi function. One of
#'     "heston" or "powerlaw".
#'
#' @return Value of constraint, which should be lower than or equal to zero to 
#'     guarantee butterfly arbitrage absence.
#'     
#' @details In a Heston parameterization constraint is read as: 
#'     $1 + |\rho| \leq 4 \gamma$, while in Power-law parameterization it is 
#'     $\eta(1+|\rho|)\leq 2$
#' @export
#'
#' @examples
#+ echo = TRUE
ssvi_butterfly_cons <- function(x, phitype) {
  # rho, gama, [eta]
  #  constraints: heston (1 + |rho|) <= 4 gamma, in power-law: eta(1+|rho|) <= 2
  names(x) <- NULL  # No need for names returned
  switch(phitype,
         heston = 1 + abs(x[1]) - 4*x[2],
         powerlaw =  x[3] * (1 + abs(x[1])) - 2
  )
}