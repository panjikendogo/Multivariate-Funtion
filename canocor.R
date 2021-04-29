canocor = function(x,y, standardize = FALSE){
  p=ncol(y)
  q=ncol(x)
  k = min(p,q)
  data=as.matrix(data.frame(y,x))
  
  if(standardize == TRUE){
    data = scale(data)
  }
  S = cov(data)
  s11 = S[1:p, 1:p]
  s12 = S[1:p, -(1:p)]
  s21 = S[-(1:p), 1:p]
  s22 = S[-(1:p),-(1:p)]
  msrs11=eigen(s11)$vectors%*%diag(1/sqrt(eigen(s11)$value))%*%t(eigen(s11)$vectors)
  msrs22=eigen(s22)$vectors%*%diag(1/sqrt(eigen(s22)$value))%*%t(eigen(s22)$vectors)
  #msr artinya min square root matriks
  rhomax=msrs11%*%s12%*%solve(s22)%*%s21%*%msrs11
  d=eigen(rhomax)$vectors
  c=msrs22%*%s21%*%msrs11%*%d
  
  # corr=diag((t(c)%*%msrs22%*%s21%*%msrs11%*%d)/(sqrt(t(c)%*%c)%*%sqrt(t(d)%*%d)))
  rhosq=solve(s11)%*%s12%*%solve(s22)%*%s21
  ev_rhosq=eigen(rhosq)$values
  corr = c()
  for (i in 1:k){
    corr[i]=sqrt(ev_rhosq[i])
  }
  
  hasil_u = matrix(0,p,k)
  for(i in (1:k)){
    hasil_u[,i] = msrs11%*%d[,i]
  }
  
  b1=matrix(0,q,k)
  for (i in (1:k)){
    b1[,i] = (solve(s22))%*%s21%*%hasil_u[,i]
  }
  
  b2 = matrix(0,k,1)
  for (i in (1:k)){
    b = b1[,i]
    b2[i,] = t(b)%*%s22%*%b
  }
  
  b11 = matrix(0,q,k)
  for (j in (1:k)){
    b11[,j] = (1/sqrt(b2[j,]))*b1[,j]
  }
  
  n=nrow(data)
  method = "Wilk's Lambda Test"
  data.name.x = deparse(substitute(x))
  data.name.y = deparse(substitute(y))
  L = (det(S))/(det(s11)*det(s22))
  chisq = -(n-(0.5*(p+q+3)))*log(L)
  df = p*q
  p.value = pchisq(chisq, df, lower.tail=FALSE)
  names(chisq) = "Chi-squared"
  names(df) = "df"
  o = structure(list(method=method, data.name=c(data.name.x , data.name.y),
                     statistic=chisq, parameter=df, p.value=p.value), class="htest")
  
  return(list(r=corr, xcoef=b11, ycoef=hasil_u, HypothesisTest=o ))
}
