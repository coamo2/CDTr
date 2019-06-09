getbias = function(Theta_real, Theta_est){
  diag_bias = sqrt(mean((diag(Theta_real) - diag(Theta_est))^2))
  diag(Theta_real) = 0
  idx = Theta_real != 0
  offdiag_bias = sqrt(mean((Theta_real[idx]-Theta_est[idx])^2))
  return(data.frame(diag_bias = diag_bias, offdiag_bias = offdiag_bias))
}
geteval = function(Theta_real, Theta_est, p, data){
  Thetabias = norm(Theta_real - Theta_est, type = "F")
  
  bias = getbias(Theta_real, Theta_est)
  diag_bias = bias$diag_bias
  offdiag_bias = bias$offdiag_bias
  
  ind = upper.tri(Theta_real,  diag = FALSE)
  TPR = sum((Theta_real != 0 & Theta_est != 0)[ind]) / sum(Theta_real[ind] != 0)
  TNR = sum((Theta_real == 0 & Theta_est == 0)[ind]) / sum(Theta_real[ind] == 0)
  
  return(data.frame(Thetabias, TPR, TNR, diag_bias, offdiag_bias))
}

pathadd = function(x){
  len = length(x)
  p = nrow(x[[1]])
  x[[len+1]] = matrix(0,p,p)
  x[[len+2]] = matrix(1,p,p)
  return(x)
}
getpath = function(cormat, cutnum = 100){
  corvec = cormat[upper.tri(cormat, diag = FALSE)]
  upper = max(abs(corvec))
  lower = min(abs(corvec))
  epsilon = (upper - lower)/10000
  upper= upper + epsilon
  lower = max(0, lower - epsilon)
  cutseq = seq(lower, upper, length.out = cutnum+1)
  diag(cormat) = 0
  path = lapply(cutseq, function(x){
    return(abs(cormat) >= x)
  })
  return(path)
}

onesim = function(Sigma, Theta, mu, p, n, seed, model){
  set.seed(seed)
  data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
  data = data / rowSums(data)
  clr.data = log(data) %*% (diag(p) - 1/p * matrix(1, p, p))
  
  Theta_offdiag = Theta
  diag(Theta_offdiag) = 0
  
  ##se(mb)
  mb.stime = Sys.time()
  mb.res = sparseiCov(clr.data, method = 'mb', lambda.min.ratio = 0.0001, nlambda = 100, sym = "or")
  mb.path = pathadd(mb.res$path)
  mb.roc = huge.roc(mb.path, Theta_offdiag, verbose = FALSE)
  mb.roc
  
  mb.res = sparseiCov(clr.data, method = 'mb', lambda.min.ratio = 0.1, nlambda = 20, sym = "or")
  mb.res = icov.select(mb.res, criterion = 'stars')
  mb.icov = symBeta(getOptBeta(mb.res), mode='ave')
  mbval = geteval(Theta_real = Theta, Theta_est = as.matrix(mb.icov), 
                  p = p, data = data)
  mb.etime = Sys.time()
  
  ##se(gl)
  gl.stime = Sys.time()
  gl.res = sparseiCov(clr.data, method = 'glasso', lambda.min.ratio = 0.0001, nlambda = 100)
  gl.path = pathadd(gl.res$path)
  gl.roc = huge.roc(gl.path, Theta_offdiag, verbose = FALSE)
  gl.roc
  
  gl.res = sparseiCov(clr.data, method = 'glasso', lambda.min.ratio = 0.001, nlambda = 15)
  gl.res = icov.select(gl.res, criterion = 'stars')
  gl.icov = gl.res$opt.icov
  glval = geteval(Theta_real = Theta, Theta_est = as.matrix(gl.icov), 
                  p = p, data = data)
  gl.etime = Sys.time()
  
  
  ###result of fang
  fang.stime = Sys.time()
  fang = gcoda(x = data, counts = F, lambda.min.ratio = 1e-3, 
               nlambda = 50, ebic.gamma = 0.5)  
  fang.path = pathadd(fang$path)
  fang.roc = huge.roc(fang.path, Theta_offdiag, verbose = FALSE)
  fang.roc
  fangval = geteval(Theta_real = Theta, Theta_est = fang$opt.icov, 
                    p = p, data = data)
  fang.etime = Sys.time()
  
  
  ###result of CDtrace bioinformatics
  cd.stime = Sys.time()
  cd <- CompDtrace(clr.data = clr.data, lambda = NULL, rho = 0.8,
                   nlambda = 20, lambda.min.ratio = 0.001)
  cd.path = pathadd(cd$path)
  cd.roc = huge.roc(cd.path, Theta_offdiag, verbose = FALSE)
  cd.roc
  cdval = geteval(Theta_real = Theta, Theta_est = cd$Theta[[which.min(cd$bic)]], 
                    p = p, data = data)
  cd.etime = Sys.time()
  
  
  ##CDtrace approximate
  cdtrace.stime_app = Sys.time()
  cdtrace.res_app = cdtrace_path(data, exact = FALSE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                             tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.path_app = lapply(cdtrace.res_app$icovpath, function(x){
    return(1*(x!=0))
  })
  cdtrace.path_app = pathadd(cdtrace.path_app)
  cdtrace.roc_app = huge.roc(cdtrace.path_app, Theta_offdiag, verbose = FALSE)
  cdtrace.roc_app
  cdtraceval_app = geteval(Theta_real = Theta, Theta_est = cdtrace.res_app$icov, 
                       p = p, data = data)
  cdtrace.etime_app = Sys.time()
  
  ##CDtrace
  cdtrace.stime = Sys.time()
  cdtrace.res = cdtrace_path(data, exact = TRUE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                          tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.path = lapply(cdtrace.res$icovpath, function(x){
    return(1*(x!=0))
  })
  cdtrace.path = pathadd(cdtrace.path)
  cdtrace.roc = huge.roc(cdtrace.path, Theta_offdiag, verbose = FALSE)
  cdtrace.roc
  cdtraceval = geteval(Theta_real = Theta, Theta_est = cdtrace.res$icov,
                    p = p, data = data)
  cdtrace.etime = Sys.time()
  
  
  maxF1_score = data.frame(mb = max(mb.roc$F1), gl = max(gl.roc$F1), 
                           fang = max(fang.roc$F1), cd = max(cd.roc$F1),
                           cdtrace_app = max(cdtrace.roc_app$F1),
                           cdtrace = max(cdtrace.roc$F1))
  AUCtab = data.frame(mb = mb.roc$AUC, gl = gl.roc$AUC, 
                      fang = fang.roc$AUC, cd = cd.roc$AUC, 
                      cdtrace_app = cdtrace.roc_app$AUC, 
                      cdtrace = cdtrace.roc$AUC)
  valtab = rbind(mbval, glval, fangval, cdval, cdtraceval_app, cdtraceval)
  rownames(valtab) = c("mb", "gl", "fang","CDtrace", "cdtrace_app", "cdtrace")
  timetab = data.frame(mb = mb.etime - mb.stime,
                       gl = gl.etime - gl.stime,
                       fang = fang.etime - fang.stime,
                       cd = cd.etime - cd.stime,
                       cdtrace_app = cdtrace.etime_app - cdtrace.stime_app,
                       cdtrace = cdtrace.etime - cdtrace.stime)
  
  save(seed, p, n, valtab, maxF1_score, AUCtab, timetab,
       fang.roc, cdtrace.roc, mb.roc, gl.roc, cdtrace.roc_app,
       file = paste0(model, "-p", p, "-n", n, "-seed", seed, ".rdata"))
  return(NULL)
}
