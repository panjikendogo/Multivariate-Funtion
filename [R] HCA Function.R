data=read.csv("hcajabar1.csv", header = T, sep = ";")

#UJI ASUMSI
attach(data)
##UJI HOMOGENTIAS VARIANS
Variabel=as.factor(rep(c("X1","X2","X3","X4","X5"),each=27))
predictor = c(X1,X2,X3,X4,X5)
data2 = data.frame(Variabel,predictor)
bartlett.test(predictor~Variabel,data2)
##UJI MULTIKOLINEARITAS
diag(solve(cor(data[,-1])))
##Semua variabel menghasilkan nilai VIF yang kurang dari 10 (Non Multikol)

#ANALISIS KLASTER
library(cluster)

HCA = function(dataset, method=c("single","complete","average","ward"), standardize = FALSE){
  data=dataset[,-1]
  obs=dataset[,1]
  rownames(data)=obs
  
  if (standardize == TRUE){
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
      a
    }
    data = baku(data)
    rownames(data)=obs
  }
  
  if(method=="single"){
    single=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=min(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    single(data)
    jarak = dist(data,method="euclidean")
    dendogram=as.dendrogram(hclust(jarak,method="single"))
    plot(dendogram,type="rectangle",main="Dendrogram: Single Linkage", ylab="Distance")
    
    hc = agnes(data, method = "single")
    if (standardize == TRUE){
      return(list(StandardizedData=data, AgglomerativeCoefficient=hc$ac))
    }
    return(list(AgglomerativeCoefficient=hc$ac))
  }
  
  else if(method=="complete") {
    complete=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=max(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    complete(data)
    jarak = dist(data,method="euclidean")
    dendogram=as.dendrogram(hclust(jarak,method="complete"))
    plot(dendogram,type="rectangle",main="Dendrogram: Complete Linkage",ylab="Distance")
    
    hc = agnes(data, method = "complete")
    if (standardize == TRUE){
      return(list(StandardizedData=data, AgglomerativeCoefficient=hc$ac))
    }
    return(list(AgglomerativeCoefficient=hc$ac))
  }
  
  else if(method=="average"){
    average=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=mean(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    average(data)
    jarak = dist(data,method="euclidean")
    dendogram=as.dendrogram(hclust(jarak,method="average"))
    plot(dendogram,type="rectangle",main="Dendrogram: Average", ylab="Distance")
    
    hc = agnes(data, method = "average")
    if (standardize == TRUE){
      return(list(StandardizedData=data, AgglomerativeCoefficient=hc$ac))
    }
    return(list(AgglomerativeCoefficient=hc$ac))
  }
  
  else if(method=="ward"){
    ward=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(data,jarak){
        data=as.matrix(data)
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        
        d1=(data[b1,]+data[k1,])/2
        data1=data[-c(b1,k1),]
        
        data=rbind(d1,data1)
        jarak=dist(data,method="euclidean")
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(data,jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]]),dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    ward(data)
    jarak = dist(data,method="euclidean")
    dendogram=as.dendrogram(hclust(jarak,method="ward.D"))
    plot(dendogram,type="rectangle",main="Dendrogram: Ward's Method",ylab="Distance")
    
    hc = agnes(data, method = "ward")
    if (standardize == TRUE){
      return(list(StandardizedData=data, AgglomerativeCoefficient=hc$ac))
    }
    return(list(AgglomerativeCoefficient=hc$ac))
  }
}
HCA(data, method="single")
HCA(data, method="single", standardize = TRUE)
HCA(data,method="complete")
HCA(data,method="complete", standardize = TRUE)
HCA(data,method="average")
HCA(data,method="average", standardize = TRUE)
HCA(data,method="ward")
HCA(data,method="ward", standardize = TRUE)
