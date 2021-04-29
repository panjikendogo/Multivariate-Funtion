#MULTIDIMENSIONAL SCALLING
mdsn=function(dataset,k){
  #Input dataset berupa data yang kolom pertama merupakan data kualitatif yang
  #menunjukkan nama dari setiap observasi dan kolom selanjutnya adalah variabel yang sudah
  #diukur
  data=dataset[,-1]
  obs=dataset[,1]
  rownames(data)=obs
  n=nrow(data)
  D=dist(data)
  matD=as.matrix(D)
  A=(-0.5)*matD^2
  I=diag(1,n)
  J=matrix(1, nrow=n, ncol=n)
  C=I-((1/n)*J)
  B=C%*%A%*%C
  ne=eigen(B)$values[1:k]
  ve=eigen(B)$vectors[,1:k]
  F=ve%*%(sqrt(diag(ne)))
  Dhat=dist(F)
  S2=(sum((D-Dhat)^2)/(sum(D^2)))
  S=sqrt(S2)*100
  dim1=ve[,1]
  dim2=ve[,2]
  plot(dim1,dim2,main="Nonmetric MDC", type="p", pch=2)
  abline(v=0)
  abline(h=0)
  text(dim1,dim2,labels=row.names(data), cex=0.7)

  while (S>10) {
    D=Dhat
    matD=as.matrix(D)
    A=(-0.5)*matD^2
    I=diag(1,n)
    J=matrix(1, nrow=n, ncol=n)
    C=I-((1/n)*J)
    B=C%*%A%*%C
    ne=eigen(B)$values[1:k]
    ve=eigen(B)$vectors[,1:k]
    F=ve%*%(sqrt(diag(ne)))
    Dhat=dist(F)
    S2=(sum((D-Dhat)^2)/(sum(D^2)))
    S=sqrt(S2)*100
    dim1=ve[,1]
    dim2=ve[,2]
    plot(dim1,dim2,main="Euclidian Distance Model", type="p", pch=2)
    abline(v=0)
    abline(h=0)
    text(dim1,dim2,labels=row.names(data), cex=0.7)
  }
  return(list(Koordinat=F,StressHat=S))
}
