
# SVI auxiliary functions -------------------------------------------------

# Author: Rafael F. Bressan
# Date: 2018-12-25
# Description: Auxiliary functions used to fit an (S)SVI implied volatility 
# model to market data.


# Libraries ---------------------------------------------------------------

library(Rsolnp)
library(nloptr)
library(ROI)
library(ROI.plugin.nloptr)

# svi_raw_fit -------------------------------------------------------------

# Inputs: [num]k, [num]w, [num]lambda, [list]init_pars
# Outputs: [list]params (a, b, rho, m, sigma)
svi_raw_fit <- function(k, w, init_pars = c(0.1, 0.1, -0.1, 0.0, 0.20), 
                        weights = NULL, solver = c("nloptr", "auto")) {
  # Selects defautl solver if not provided
  solver <- solver[1]
  # Initial checks
  if (is.null(weights)) {
    lambda <- 1
  } else if (length(lambda) != length(k)) {
    stop("Lambda must be either NULL or a numeric vector of the same size as k.")
  } else {
    lambda <- weights
  }
  # Number of parameters
  n <- length(init_pars)
  
  if (n != 5L) stop("Number of parameters must be 5.")
  
  w_mkt <- w
  
  # Objective function
  # first parameter must be a vector of parameters to be minimized on
  obj_fun <- function(par) {
    a <- par[1]
    b <- par[2]
    rho <- par[3]
    m <- par[4]
    sigma <- par[5]
    
    # SVI w
    w_svi <- a + b*(rho*(k - m) + sqrt((k - m)^2 + sigma^2))
    
    obj <- (1/sum(lambda)*sum(lambda*(w_svi - w_mkt)^2))
    return(obj)
  }
  
  # Gradient function
  grad_fun <- function(par) {
    # Finite difference to approximate Gradient
    a <- par[1]
    b <- par[2]
    rho <- par[3]
    m <- par[4]
    sigma <- par[5]
    eps <- 1.0e-8
    
    # Partials with respect to each parameter
    g1 <- (obj_fun(c(a + eps, b, rho, m, sigma)) - obj_fun(par)) / obj_fun(par)
    g2 <- (obj_fun(c(a, b + eps, rho, m, sigma)) - obj_fun(par)) / obj_fun(par)
    g3 <- (obj_fun(c(a, b, rho + eps, m, sigma)) - obj_fun(par)) / obj_fun(par)
    g4 <- (obj_fun(c(a, b, rho, m + eps, sigma)) - obj_fun(par)) / obj_fun(par)
    g5 <- (obj_fun(c(a, b, rho, m, sigma + eps)) - obj_fun(par)) / obj_fun(par)
    
    return(rbind(c(g1, g2, g3, g4, g5)))
  }
  
  # Constrain function
  const_fun <- function(par){
    return(-par[1] - par[2]*par[5]*sqrt(1 - par[3]^2))
  }
  
  # Jacobian of contrain function
  jacob_fun <- function(par) {
    j1 <- -1
    j2 <- -par[5]*sqrt(1 - par[3]^2)
    j3 <- par[2]*par[3]*par[4]/sqrt(1 - par[3]^2)
    j4 <- 0
    j5 <- -par[2]*sqrt(1 - par[3]^2)
    return(rbind(c(j1, j2, j3, j4, j5)))
  }
  
  # Parameters bounds
  par_lower = c(-Inf, 0, -1, -Inf, 0)
  par_upper = c(Inf, Inf, 1, Inf, Inf)
  
  # Otimizadores ######################################################
  # stats::optim() stats::nlminb()
  # Pacotes: optimx, lbfgsb3c, Rsolnp, nlsr, DEoptim(R), GA
  #####################################################################
  NM <- optim(init_pars, obj_fun, method = "Nelder-Mead")
  LBFGSB <- optim(init_pars, obj_fun, method = "L-BFGS-B",
              lower = par_lower,
              upper = par_upper)
  
  # Using nloptr
  # nl_optr <- nloptr(x0 = init_pars, 
  #                   eval_f = obj_fun, 
  #                   eval_grad_f = grad_fun,
  #                   eval_g_ineq = const_fun,
  #                   eval_jac_g_eq = jacob_fun,
  #                   lb = par_lower, 
  #                   ub = par_upper, 
  #                   opts = list("algorithm" = "NLOPT_LD_MMA",
  #                               "xtol_rel" = 1.0e-8,
  #                               "print_level" = 2,
  #                               "check_derivatives" = TRUE,
  #                               "check_derivatives_print" = "all"))
  nl_optr <- nloptr(x0 = init_pars, 
                    eval_f = obj_fun, 
                    eval_g_ineq = const_fun,
                    lb = par_lower, 
                    ub = par_upper, 
                    opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                "xtol_rel" = 1.0e-8))
  
  # Using Rsolnp
  go_nlp <- gosolnp(fun = obj_fun,
                    ineqfun = const_fun,
                    ineqLB = -Inf,
                    ineqUB = 0,
                    LB = par_lower,
                    UB = par_upper,
                    n.restarts = 4,
                    distr = c(2, 2, 2, 2, 2),
                    distr.opt = list(list(mean = init_pars[1], sd = 1),
                                     list(mean = init_pars[2], sd = 1),
                                     list(mean = init_pars[3], sd = 1),
                                     list(mean = init_pars[4], sd = 1),
                                     list(mean = init_pars[5], sd = 1)),
                    n.sim = 100000)
  Rnlp <- solnp(pars = init_pars,
                fun = obj_fun,
                ineqfun = const_fun,
                ineqLB = -Inf,
                ineqUB = 0,
                LB = par_lower,
                UB = par_upper)
  
  method_names <- c("NM", "LBFGSB", "NLOPTR", "GONLP", "RNLP")
  pars <- list(NM$par, LBFGSB$par, nl_optr$solution, go_nlp$pars, Rnlp$pars)
  const_value <- sapply(pars, const_fun)
  return(tibble(opt_method = method_names,
                pars = pars,
                obj_value = c(NM$value, LBFGSB$value, nl_optr$objective, 
                              go_nlp$values[go_nlp$outer.iter + 1], 
                              Rnlp$values[Rnlp$outer.iter + 1]), 
                convergence = c(NM$convergence, LBFGSB$convergence, 
                                nl_optr$status, 
                                go_nlp$convergence, Rnlp$convergence),
                const_value = const_value
  )
  )
  # return(list(NM = NM, LBFGSB = LBFGSB, GONLP = go_nlp, RNLP = Rnlp))
}

# svi_raw_simulate --------------------------------------------------------

# Inputs: [num]k, [list]params (a, b, rho, m, sigma)
# Outputs: [num]w
svi_raw_simulate <- function(par, k) {
  a <- par[1]
  b <- par[2]
  rho <- par[3]
  m <- par[4]
  sigma <- par[5]
  # SVI w
  w_svi <- a + b*(rho*(k - m) + sqrt((k - m)^2 + sigma^2))
  return(w_svi)
}

# svi_invert_w ------------------------------------------------------------

# Inputs: [num]w, [num]t
# Outputs: [num]sigma_imp
svi_invert_w <- function(w, t) {
  # t is a vector the same size as w 
  if (length(t) != length(w)) stop("t must be a vector of the same size as w.")
  sigma_imp <- sqrt(w / t)
  return(sigma_imp)
}


# Algumas variaveis para testes
nomes <- c("a", "b", "rho", "m", "sigma")
a <- 0.04
b <- 0.13
rho <- -1
m <- 0.06
sigma <- 0.20
restricao <- a + b*sigma*sqrt(1 - rho^2)
parametros <- c(a, b, rho, m, sigma)
names(parametros) <- nomes
k <- seq(-0.5, 0.2, by = 0.1)
w_mkt <- a + b*(rho*(k - m) + sqrt((k - m)^2 + sigma^2))
t <- rep(1, length(k))
sigma_mkt <- sqrt(w_mkt / t)
lambda <- c(0.4, 0.5, 0.7, 0.9, 1.2, 1.6, 1.2, 0.7)
mkt_df <- tibble(name = rep("MKT", length(k)),
                 w_mkt = w_mkt, 
                 k = k,
                 t = t,
                 sigma_mkt = sigma_mkt,
                 weights = lambda)

# init_pars <- list(a = 0.0,
#                   b = 0.10,
#                   rho = 0.0,
#                   m = 0.0,
#                   sigma = 0.15)

# Fitting -----------------------------------------------------------------

fit <- svi_raw_fit(mkt_df$k, mkt_df$w_mkt, weights = mkt_df$weights)

# Ploting the Fit ---------------------------------------------------------

fit_df <- fit %>% 
  mutate(w_sim = map(pars, ~svi_raw_simulate(.x, mkt_df$k)),
         sigma_sim = map(w_sim, ~svi_invert_w(.x, mkt_df$t)),
         k = list(mkt_df$k)) %>% 
  select(opt_method, sigma_sim, k) %>% 
  unnest() %>% 
  gather(key = measure, value = value, -c(opt_method, k)) 

fit_df %>% 
  filter(measure == "sigma_sim") %>% 
ggplot(aes(x = k, y = value, color = opt_method)) +
  geom_line() +
  geom_point(data = mkt_df, aes(x = k, y = sigma_mkt, color = name))
