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

dgompertz = function(x, a, b, ln = F) {
  if (a <= 0){
    stop("'a' must be a positive value.")
  }
  else if (b <= 0){
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
  if (a <= 0){
    stop("'a' must be a positive value.");
  }
  else if (b <= 0){
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
