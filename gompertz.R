survfit_gompertz = function(y, d, dados = NULL) {
  log_like = function(par) {
    a = par[1]
    b = par[2]
    return(sum(d*dgompertz(y, a, b, T) + (1-d)*(pgompertz(y, a, b, T, lower.tail = F))))
  }
  
  start = c(1,1)
  
  param = optim(start, log_like, control = list(fnscale = -1, maxit = 500),
                method="L-BFGS-B", lower = c(0.0001, 0.0001), upper = c(Inf,Inf))
  
  print(param$value)
  return(param$par)
}

#Créditos: função retirada do pacote flexsurv
.hess_to_cov <- function(hessian, tol.solve = 1e-9, tol.evalues = 1e-5, ...) {
  if(is.null(tol.solve)) tol.solve <- .Machine$double.eps
  if(is.null(tol.evalues)) tol.evalues <- 1e-5 
  # use solve(.) over chol2inv(chol(.)) to get an inverse even if not PD
  # less efficient but more stable
  inv_hessian <- solve(hessian, tol = tol.solve)
  evalues <- eigen(inv_hessian, symmetric = TRUE, only.values = TRUE)$values
  if (min(evalues) < -tol.evalues)
    warning(sprintf(
      "Hessian not positive definite: smallest eigenvalue is %.1e (threshold: %.1e). This might indicate that the optimization did not converge to the maximum likelihood, so that the results are invalid. Continuing with the nearest positive definite approximation of the covariance matrix.",
      min(evalues), -tol.evalues
    ))
  # make sure we return a plain positive definite symmetric matrix
  as.matrix(Matrix::nearPD(inv_hessian, ensureSymmetry = TRUE, ...)$mat)
}


## Funções densidade, sobrevivência, acumulada, quantil e de geração aleatória
dgompertz = function(x, a, b, ln = F) {
  if(min(a) <= 0){
    stop("'a' must be a positive value.")
  }
  else if (min(b) <= 0){
    stop("'b' must be a positive value.")
  }
  else if (ln == T){
      out = log(a) - log(b) + x/b -a*(expm1(x/b))
  }
  else {
      out = (a/b)*exp(x/b)*exp(-a*(expm1(x/b)))
  }
  return(out)
}

pgompertz = function(x, a, b, ln = F, lower.tail = T) {
  if (min(a) <= 0){
    stop("'a' must be a positive value.");
  }
  else if (min(b) <= 0){
    stop("'b' must be a positive value.");
  }
  else if (lower.tail == T){
    if (ln == T){
        out = log(1 - exp(-a*expm1(x/b)));
    }
    else {
        out = 1 - exp(-a*(expm1(x/b)));
    }
  }
  else {
    if (ln == T){
        out = -a*expm1(x/b);
      }
    else {
        out = exp(-a*(expm1(x/b)));
    }
  }
  return(out);
}

qgompertz = function(x, a, b) {
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else {
    out = ifelse((x <= 1) & ( x > 0), b*log(1 - (1/a)*log(1-x[i])),NaN)
  }
  return (out);
}

rgompertz = function(n, a, b) {
  p = runif(n, 0, 1)
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
    stop("'b' must be a positive value.");
  }
  else {
      out =  b*log1p(-(1/a)*log1p(-p));
  }
  return (out);
}
