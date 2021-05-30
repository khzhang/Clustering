### Get bandwidth
h.dD = function(D){
  N = nrow(D)
  d = dim(D)[2]
  D_r = c()
  sd.sum = 0
  r = matrix(nrow = d,ncol = 2)
  IQR = c()
  
  for (i in 1:d){
    D_r[i] = max(D[,i])-min(D[,i])
  }
  D_rm = min(D_r)
  
  for (i in 1:d){
    sd.sum = sd(D[,i])/D_r[i]+sd.sum
    r[i,] = quantile(D[,i], c(0.25, 0.75))
    IQR[i] = (r[i,2]-r[i,1])/D_r[i]*D_rm
  }
  
  sd.r = sd.sum /d * D_rm
  IQR.r = mean(IQR)
  con = min(IQR.r/1.34, sd.r)
  return((1/N)^(1/(8+d))*con)
}

### Get gradient of f
grad = function(X, gr, h){
  n = dim(X)[1]
  d = dim(X)[2]
  m = dim(gr)[1]
  u = list()
  result = matrix(NA, nrow=m, ncol=d)
  
  gr = matrix(gr, ncol = d, byrow = F)
  for(i in 1:m){
    S = 1
    
    for(j in 1:d){
      u[[j]] = (X[,j]-gr[i,j])/h
      S =  dnorm(u[[j]])*S
    }
    
    for(j in 1:d){
      result[i,j] = 1/(n*h^(d+2))*sum(u[[j]]* h * S)
    }
    
  }
  return(as.matrix(result))
}

### Get Hessian of f
Hess = function(X, gr, h){
  
  n = dim(X)[1]
  result = list()
  d = dim(X)[2]
  gr = matrix(gr, ncol = d, byrow = F)
  m = dim(gr)[1]
  u = list()
  Hes = matrix(0, nrow = d, ncol=d)
  
  for(i in 1:m){
    
    S = 1
    for(j in 1:d){
      u[[j]] = (X[,j]-gr[i,j])/h
      S =  dnorm(u[[j]])*S
    }
    
    for(j in 1:d){
      Hes[j,j] = 1/(n*h^d)*sum(((u[[j]]/h)^2 - 1/(h^2))* S) 
      for(k in 1:d){
        if(k != j){
          Hes[k,j] = 1/(n*h^d)* sum((u[[j]]/h * u[[k]]/h) * S) 
        }
      }
    }
    result[[i]] = Hes
  }
  return(result)
}

### Algorithmn 1
GD_flow2 = function(X, gr, h, k.max = 12000, eps = 1e-12){
  gr = matrix(gr, ncol = 2,byrow = F)
  d = ncol(gr)
  x_test_old = gr
  x_test = matrix(0,dim(x_test_old)[1],dim(x_test_old)[2])
  
  test_Hess = Hess(X, gr = matrix(x_test_old,ncol = d, byrow = F), h=h)
  stepsize = min(unlist(lapply(1:dim(x_test_old)[1], function(i){1/max(eigen(test_Hess[[i]] %*% test_Hess[[i]])$values)})))
  w = 1:nrow(x_test_old)
  k = 0
  while(1){
    if (k == k.max || length(w) == 0){
      break
    }
    
    if (k > 0){
      test_Hess = Hess(X, gr = matrix(x_test,ncol = d, byrow = F), h=h)
      stepsize = abs(min(unlist(lapply(1:dim(x_test)[1], function(i){1/max(eigen(test_Hess[[i]] %*% test_Hess[[i]])$values)}))))
    }
    
    test_grad_w = grad(X, gr=matrix(x_test_old[w,],ncol = d, byrow = F), h=h) 
    test_Hess_w = Hess(X, gr = matrix(x_test_old[w,],ncol = d, byrow = F), h=h) 
    
    flowChange2 = lapply(1:length(w),function(i){
      test_grad_w [i,] %*% test_Hess_w[[i]]
    })
    flowChange2 = stepsize * matrix(unlist(flowChange2),ncol = d, byrow = T)
    
    x_test = x_test_old[w,] - flowChange2
    tmp = x_test_old
    x_test_old[w,] = x_test
    w_err = rowSums(x_test_old - tmp)^2
    
    w = which(w_err>eps)
    k = k+1
  }
  return(x_test_old)
}
