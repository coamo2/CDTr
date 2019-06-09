######assumption exactly/approxiamtely right#########
starttime = Sys.time()

path = '/data/heshun/CDtrace/CDtracesim'
#path = "C:\\Users\\asus\\Desktop\\CDtrace\\CDtracesim"
setwd(path)
source("gCoda-master/R/gcoda.R")
source("bioinfo_CDtrace/mycd.R")
require(huge)
source("simfunc.R")

require(doSNOW)
require(parallel)
cl = makeSOCKcluster(rep("localhost", times = 6))
registerDoSNOW(cl)

p = 50
n = 150
asbias_all = NULL
 model = "band-exact"
# model = "band-approx1"
# model = "band-approx2"
model = "cluster-exact"
model = "cluster-approx1"
model = "cluster-approx2"
for(model in c("band-exact", "band-approx1", "band-approx2", 
               "cluster-exact", "cluster-approx1", "cluster-approx2")){
  if(model == "cluster-exact"){
    g = 5; upper = 0.1; lower = 0.1; model_gen = "cluster"
  }else if(model == "cluster-approx1"){
    g = 5; upper = 0.15; lower = 0.05; model_gen = "cluster"
  }else if(model == "cluster-approx2"){
    g = 5; upper = 0.2; lower = 0; model_gen = "cluster"
  }else if(model == "band-exact"){
    g = 2; upper = 0.1; lower = 0.1; model_gen = "band"
  }else if(model == "band-approx1"){
    g = 2; upper = 0.15; lower = 0.05; model_gen = "band"
  }else if(model == "band-approx2"){
    g = 2; upper = 0.2; lower = 0; model_gen = "band"
  }
  ##generate network
  if(model_gen == "cluster"){
    set.seed(1)
    Theta = matrix(0,p,p)
    subcl = matrix(1, p/g, p/g)
    subcl[1:(p/g/2),1:(p/g/2)] = -1
    subcl[(p/g/2+1):(p/g),(p/g/2+1):(p/g)] = -1
    diag(subcl) = 0
    for(i in 1:g){
      Theta[((i-1)*p/g+1):(i*p/g),((i-1)*p/g+1):(i*p/g)] = subcl
    }
    Theta = Theta * matrix(runif(p^2, lower, upper), p, p)
    Theta[lower.tri(Theta)] = 0
    Theta = Theta + t(Theta)
  }else if(model_gen == "band"){
    set.seed(1)
    Theta = matrix(0,p,p)
    onecol = c(0,1,-1, rep(0, times = p-1-2*g), -1,1)
    for(i in 1:(p-1)){
      Theta[i,(i+1):p] = onecol[2:(p-i+1)]
    }
    Theta = Theta * matrix(runif(p^2, lower, upper), p, p)
    Theta[lower.tri(Theta)] = 0
    Theta = Theta + t(Theta)
  }
  diag(Theta) = abs(min(eigen(Theta)$values))+0.3
  
  Sigma = solve(Theta)
  mu = runif(n = p, min  = -0.5, max = 0.5)
  
  GG = diag(p) - 1/p*matrix(1,p,p)
  asbias = c(norm(GG %*% Theta - Theta %*% GG, 'F'), norm(GG %*% Sigma - Sigma %*% GG, 'F'))
  asbias
  asbias_all = rbind(asbias_all,c(model,asbias))
  
  setwd(path)
  modelpath = paste0(path, "/", model)
  dir.create(modelpath, showWarnings = FALSE)
  setwd(modelpath)

  save(Theta, Sigma, mu, asbias, file = paste0(model, "-p", p, "-net.rdata"))

  for(n in c(50, 100, 150, 200)){
    cat("graph:", model, "  n:", n, "\n")
    void = foreach(seed = 1:200)%dopar%{
      path1 = paste0(path, "/SpiecEasi-master/R")
      setwd(path1)
      for(f in list.files(path1)){
        source(f)
      }
      path2 = paste0(path, "/CDtrace")
      setwd(path2)
      for(f in list.files(path2)){
        source(f)
      }
      setwd(path)
      source("gCoda-master/R/gcoda.R")
      source("bioinfo_CDtrace//mycd.R")
      require(huge)
      source("simfunc.R")

      #seed=1
      setwd(modelpath)
      errind = try(onesim(Sigma = Sigma, Theta = Theta, mu = mu,
                          p = p, n = n, seed = seed, model = model),
                   silent = TRUE)

      if(inherits(errind, "try-error")){
        return(seed)
      }else{
        return(0)
      }
    }

    setwd(modelpath)
    save(void, file = paste0(model, "-p", p, "-n", n, "-errseed.rdata"))
  }
}

stopCluster(cl)
endtime = Sys.time()
setwd(path)
save(asbias_all, file = paste0(p,"-asbias-right.rdata"))






######assumption not right####
starttime = Sys.time()
path = '/data/heshun/CDtrace/CDtracesim'
#path = "C:/Users/asus/Desktop/CDtrace/CDtracesim"
setwd(path)
source("gCoda-master/R/gcoda.R")
source("bioinfo_CDtrace/mycd.R")
require(huge)
source("simfunc.R")

require(doSNOW)
require(parallel)
cl = makeSOCKcluster(rep("localhost", times = 10))
registerDoSNOW(cl)

asbias_all = NULL
p = 50
n = 200
#model = "neighbor"
#model = "block"
#model = "random"
#model = "hub"
# model = "scale-free"
# model = "band"
 model = "cluster"
for(model in c("random", "hub", "neighbor", "block", "band", "scale-free")){
  if(model == "random"){
    g = NULL; prob = 0.1; upper = 0.2; lower = 0.1
  }else if(model == "hub"){
    g = 3; prob = NULL; upper = 0.2; lower = 0.1
  }else if(model == "cluster"){
    g = 3; prob = 0.1; upper = 0.2; lower = 0.1
  }else if(model == "band"){
    g = 4; prob = NULL; upper = NULL; lower = NULL
  }else if(model == "scale-free"){
    g = NULL; prob = NULL; upper = 0.2; lower = 0.1
  }
  if(model == "neighbor"){
    set.seed(1)
    Theta = matrix(0,p,p)
    locmat = matrix(0, p, p)
    d  = runif(n = p, 0, 1)
    for(i in 1:p){
      x = abs(d[i] - d)
      loc = order(x)[2:6]
      locmat[i, loc] = 1; locmat[loc, i] = 1
    }
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(locmat[i,j] == 1){
          val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
            runif(1, 0.05, 0.15)
          Theta[i, j] = val; Theta[j, i] = val
        }
      }
    }
  }else if(model == "block"){
    set.seed(1)
    Theta = matrix(0,p,p)
    block = sample(p, size = p, replace = FALSE)
    blocksize = p / 5
    block = lapply(1:5, function(i){
      sort(block[(blocksize*(i-1)+1):(blocksize*i)])
    })
    
    for(i in 1:5){
      for(j in block[[i]]){
        for(k in block[[i]][block[[i]] > j]){
          if(rbinom(n = 1, size = 1, prob = 0.3)){
            val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
              runif(n = 1, 0.1, 0.2)
            Theta[k,j] = val
            Theta[j,k] = val
          }}}}
    
    for(i in 1:4){
      for(j in block[[i]]){
        for(k in unlist(lapply(i:5, function(l)block[[l]]))){
          if(rbinom(n = 1, size = 1, prob = 0.1)){
            val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
              runif(n = 1, 0.1, 0.2)
            Theta[k,j] = val
            Theta[j,k] = val
          }}}}
  }else{
    ##generate network
    set.seed(1)
    net_generator = huge.generator(n = 2*p, d = p, graph = model, v = NULL, u = NULL, 
                                   g = g, prob = prob, vis = FALSE, verbose = FALSE)
    #plot(net_generator, align = F)
    net = as.matrix(net_generator$theta)
    
    if(model == "band"){
      link_strength = matrix(0,p,p)
      for(i in 1:(p-1))link_strength[i,i+1] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.2, 0.25)
      for(i in 1:(p-2))link_strength[i,i+2] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.15, 0.2)
      for(i in 1:(p-3))link_strength[i,i+3] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.1, 0.15)
      for(i in 1:(p-4))link_strength[i,i+4] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.05, 0.1)
    }else{
      link_strength = matrix((2 * rbinom(n = p^2, size = 1, prob = 0.5) - 1) * 
                               runif(n = p^2, lower, upper), p, p)
    }
    link_strength[lower.tri(link_strength, diag = TRUE)] =  0
    link_strength = link_strength + t(link_strength)
    
    Theta = net * link_strength
  }
  
  diag(Theta) = abs(min(eigen(Theta)$values))+0.3
  
  Sigma = solve(Theta)
  mu = runif(n = p, min  = -0.5, max = 0.5)
  
  GG = diag(p) - 1/p*matrix(1,p,p)
  asbias = c(norm(GG %*% Theta - Theta %*% GG, 'F'), norm(GG %*% Sigma - Sigma %*% GG, 'F'))
  asbias_all = rbind(asbias_all,c(model,asbias))
  
  setwd(path)
  modelpath = paste0(path, "/", model)
  dir.create(modelpath, showWarnings = FALSE)
  setwd(modelpath)

  save(Theta, Sigma, mu, asbias, file = paste0(model, "-p", p, "-net.rdata"))

  for(n in c(50, 100, 150, 200)){
    cat("graph:", model, "  n:", n, "\n")
    void = foreach(seed = 1:200)%dopar%{
      path1 = paste0(path, "/SpiecEasi-master/R")
      setwd(path1)
      for(f in list.files(path1)){
        source(f)
      }
      path2 = paste0(path, "/CDtrace")
      setwd(path2)
      for(f in list.files(path2)){
        source(f)
      }
      setwd(path)
      source("gCoda-master/R/gcoda.R")
      source("bioinfo_CDtrace/mycd.R")
      require(huge)
      source("simfunc.R")

      #seed=1
      setwd(modelpath)
      errind = try(onesim(Sigma = Sigma, Theta = Theta, mu = mu,
                          p = p, n = n, seed = seed, model = model),
                   silent = TRUE)

      if(inherits(errind, "try-error")){
        return(seed)
      }else{
        return(0)
      }
    }

    setwd(modelpath)
    save(void, file = paste0(model, "-p", p, "-n", n, "-errseed.rdata"))
  }
}

stopCluster(cl)
endtime = Sys.time()
setwd(path)
save(asbias_all, file = paste0(p,"-asbias-notright.rdata"))
