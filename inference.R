library(TDA)

CenterFind = function(X, x_test_i, n_grid = 51){
  h = h.dD(X)
  gr0 = expand.grid(seq(from=min(X[,1]),to=max(X[,1]), length.out=n_grid),
                    seq(from=min(X[,2]),to=max(X[,2]), length.out=n_grid))
  gr0 = as.matrix(gr0)
  # gr0_kde = kde(X, Grid=gr0, h=h)
  
  X_den0 = bkde2D(X, h, c(n_grid, n_grid))
  Hession = Hess(X, gr0, h)
  gr_Hess_det = unlist( lapply(1:dim(gr0)[1],function(i){det(Hession[[i]])})  )
  
  
  D_den = X_den0$fhat
  CL_den = contourLines(D_den, levels=(0:60)/60 *max(summary(c(D_den))),x=X_den0$x1,y=X_den0$x2)
  
  CL_denInd = x_test_i[,1] < max(CL_den[[15]]$x) & x_test_i[,1] > min(CL_den[[15]]$x)
  upperPart = which.min(CL_den[[15]]$x):which.max(CL_den[[15]]$x)
  downPart = (1:length(CL_den[[15]]$x))[-upperPart]
  
  tmpup = nn2( CL_den[[15]]$x[upperPart],x_test_i[CL_denInd,1],k=1)
  tmpdown = nn2(CL_den[[15]]$x[downPart],x_test_i[CL_denInd,1],k=1)
  
  downOutInd = which(x_test_i[CL_denInd,2] - CL_den[[15]]$y[downPart[tmpdown$nn.idx]] <0)
  upOutInd = which(x_test_i[CL_denInd,2] - CL_den[[15]]$y[upperPart[tmpup$nn.idx]] >0)
  
  # the centerpart
  x_test_inarea = x_test_i[CL_denInd,][-c(downOutInd,upOutInd),]
  Usefulind = c(1:dim(x_test_i)[1])[CL_denInd][-c(downOutInd,upOutInd)]
  
  return(list(x_test_inarea, Usefulind))
}


CenterFind2 = function(X, x_test_i){
  quantileV = apply(X,2,function(t){quantile(t,probs = c(0.05,0.95))})
  
  # the centerpart
  x_test_inarea = x_test_i[(x_test_i[,1] >= quantileV[1,1]) & 
                             (x_test_i[,1] <= quantileV[2,1]) &
                             (x_test_i[,2]  >= quantileV[1,2])&
                             (x_test_i[,2] <= quantileV[2,2]) 
                             ,]
  Usefulind = c(1:dim(x_test_i)[1])[(x_test_i[,1] >= quantileV[1,1]) & 
                                      (x_test_i[,1] <= quantileV[2,1]) &
                                      (x_test_i[,2]  >= quantileV[1,2])&
                                      (x_test_i[,2] <= quantileV[2,2]) ]
  
  return(list(x_test_inarea, Usefulind))
}


ProportionFinder = function(x_test_i, x_test_inarea, Usefulind,h){
  hclust_x_test_inarea = hclust(dist(x_test_inarea), method="single")
  gr1_clu = cutree(hclust_x_test_inarea,h = 1.2*h)
  #useful cluster id
  notusefulId =which(table(gr1_clu) <=1)
  
  ######calculate distance to centers for each group
  N = dim(x_test_i)[1]
  Xlab = rep(10, N)
  Xlab[Usefulind] = gr1_clu
  Xlab[Xlab %in% notusefulId] = 10
  X_lab1 = Xlab[1:(N/2)]
  X_lab2 = Xlab[(N/2+1):N]
  
  grp1Count = data.frame(table(X_lab1))
  grp2Count = data.frame(table(X_lab2))
  grp1Name = rownames(grp1Count)
  grp2Name = rownames(grp2Count)
  grpName = c(grp1Name,grp1Name)
  contigenTable = merge(grp1Count, grp2Count,by.x = "X_lab1",by.y = "X_lab2",all = T)
  contigenTable[is.na(contigenTable)] = 0
  
  return(contigenTable)
}

zFun = function(p,p0,ni){
  zscore = (p - p0)/sqrt(p0*(1 - p0)/ni)
  pval = dnorm(zscore)
  return(cbind(zscore,pval))
}

TwosampeTest = function(X,x_test_i){
  N = dim(X)[1]
  h = 1.2 * h.dD(X)
  ress = CenterFind2(X, x_test_i)
  restmp = ProportionFinder(x_test_i, ress[[1]], ress[[2]], h)
  Prop = restmp[restmp$X_lab1 != 10,]
  print(dim(Prop)[1])
  N1 = N/2
  N2 = N/2
  p0 = N1/(N1 + N2)
  ni = Prop$Freq.x + Prop$Freq.y
  pProp = Prop$Freq.x /(Prop$Freq.x + Prop$Freq.y)
    
  pvalue = c()
  for (i in 1:dim(Prop)[1]){
    pvalue[i] = zFun(pProp[i],p0,ni[i])[2]
  }
  # if true, reject H0
  holm = sum(p.adjust(pvalue , method="holm") <= 0.05) >0
  Fdr = sum(p.adjust(pvalue , method="BH") <= 0.05) >0
  return(list(holm,Fdr))
}

ciFun = function(p, ni, zc = qnorm(0.975)){
  lower = p - zc * sqrt(p*(1 - p)/ni)
  upper = p + zc * sqrt(p*(1 - p)/ni)
  if (upper > 1){
    upper = 1
  }
  if(lower <0){
    lower = 0
  }
  c(lower,upper)
}


