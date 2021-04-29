data=read.csv("diskriminan1.csv",header=T)
data
y=data[,1]
x=data[,2:7]

##Variabel y berisi nilai 1 dan 2 yang menunjukkan kategori
andis2=function(x,y){
  
  data=data.frame(x,y)
  n=nrow(data)
  p=ncol(data)-1
  
  k1=data[which(y==1),][1:p]
  k2=data[which(y==2),][1:p]
  
  n1=nrow(k1)
  n2=nrow(k2)
  mu1=colMeans(k1)
  mu2=colMeans(k2)
  S1=cov(k1)
  S2=cov(k2)
  Spl=((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
  
  a=solve(Spl)%*%(mu1-mu2)
  
  z1=as.matrix(k1)%*%a
  z2=as.matrix(k2)%*%a
  
  z1.bar=mean(z1)
  z2.bar=mean(z2)
  m=(n1*z2.bar+n2*z1.bar)/(n1+n2)
  
  kat1=matrix(0,n1,1)
  for (i in (1:n1)){
    kat1[i]=if(z1[i] >= m) 1 else 2
  }
  kat2=matrix(0,n2,1)
  for (i in (1:n2)){
    kat2[i]=if(z2[i] >= m) 1 else 2
  }
  
  b=length(kat1[which(kat1==1)])
  c=length(kat1[which(kat1==2)])
  d=length(kat2[which(kat2==1)])
  e=length(kat2[which(kat2==2)])
  
  Kes=matrix(c(b,c,d,e),2,2,byrow = F)
  dimnames(Kes)=list(Estimasi=c("K1","K2"),Observasi=c("K1","K2"))
  
  hitratio=(b+e)/n
  
  hasil=list("Koefisien Kombinasi Linier Fisher"=a,"Skor Diskriminan Kelompok 1"=z1,
             "Skor Diskriminan Kelompok 2"=z2,"Cutting Scor"=m,
             "Matriks Estimasi dan Observasi"=Kes,"Hit Ratio"=hitratio)
  print(hasil)
}

andis2(x,y)

