#ini script untuk function pca
PCA=function(x, standardize=FALSE){
  data=as.matrix(x)
  S = cov(data)
  if(standardize == TRUE){
    S = cor(data)
  }
  eigen_val = eigen(S)$values
  eigen_vec = eigen(S)$vector
  
  n = length(eigen_val)
  prop = c()
  for (i in (1:n)){
    prop[i] = (eigen_val[i])/sum(eigen_val)
  }
  
  q = length(prop)
  propcum = c()
  for (i in 1:q){
    propcum[i]=sum(prop[1:i])
  }
  
  p = nrow(eigen_vec)
  corr = matrix(0,p,p)
  for (i in (1:p)){
    for (j in (1:p)){
      corr[i,j] = (eigen_vec[i,j]*sqrt(eigen_val[i]))/S[j,j]
    }
  }
  
  plot(eigen_val, main="Scree Plot", type="o")
  hasil = list("Matriks Varkov/Korelasi"=S, "Eigen Value"=eigen_val, 
               "Eigen Vector"=eigen_vec, "Proporsi Komponen"=prop, 
               "Proporsi Kumulatif Komponen"=propcum, "Matriks Korelasi Y dan X"=corr)
  print(hasil)
}

baku=function(X){
  x=as.matrix(X)
  n=nrow(x)
  p=ncol(x)
  a=matrix(0,n,p)
  for (j in (1:p)) {
    for (i in (1:n)) {
      a[i,j]=(x[i,j]-mean(x[,j]))/sd(x[,j])
    }
  }
  print(a)
}

data=read.csv("DATA PCA.csv", sep=";")
PCA(data, standardize=T)
library(MVN)