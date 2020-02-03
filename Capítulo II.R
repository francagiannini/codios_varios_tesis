# Multispati-PCA

setwd("C:/Users/Franca/Dropbox/Franca/Doctorado/Analisis de datos/SIG-PID")
pred=read.table("predsinlim.txt", header=T, sep= "\t")
install.packages("spdep")
library(ade4)
library(spdep)

pca=dudi.pca(pred[,3:13],center=T, scannf=F, nf=10)
coor_2=coordinates(pred[,1:2])
gri_2=dnearneigh(coor_2,0,50000)
lw_2=nb2listw(gri_2, style="W")
ms2=multispati(pca,lw_2,scannf=F, nfposi = 10)

pred_cpssinlim=cbind(pred,ms2$li)
a=summary(ms2)
eigenvalues=ms2$eig
b=summary(pca)
ms2$ls
ms2$c1

##nbclust

install.packages("NbClust")
library(NbClust)
setwd("C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict")
getwd()

pict=componentespict
class(pict)
str(pict)
p=pict[,26:33]
p
as.matrix(p)
res=NbClust(data =p, diss = NULL, distance = "euclidean", 
            min.nc =2 , max.nc = 15,method = "kmeans", index = "all", alphaBeale = 0.1)
All.index=res$All.index
Best.nc=res$Best.nc
All.CriticalValues=res$All.CriticalValues
Best.partition=res$Best.partition

#exportar
save(All.index, file = 
       "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/All.index")
save (Best.nc, file=
        "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/Best.nc")
save(All.CriticalValues,file=
       "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/All.CriticalValues")
save (Best.partition, file=
        "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/Best.partition")

set.seed(1234)
# Solución con dos conglomerados
C2<- kmeans(p, 2);C2
C3<- kmeans(p, 3);C3 
C4<- kmeans(p, 4);C4
C5<- kmeans(p,5);C5 
C6<- kmeans(p, 6);C6 
C7<- kmeans(p, 7);C7 
C8<- kmeans(p, 8);C8 
C9<- kmeans(p, 9);C9 
C10<-kmeans(p, 10);C10 
C11<-kmeans(p, 11);C11 
C12<-kmeans(p, 12);C12 
C13<-kmeans(p, 13);C13 
C14<-kmeans(p, 14);C14 
C15<-kmeans(p, 15);C15 

# conjunto de datos la clasificación para los distintos grupos

pictclsificacionNbclust<-data.frame(pict,C2$cluster,C3$cluster,C4$cluster,C5$cluster,C6$cluster,C7$cluster,C8$cluster,C9$cluster,C10$cluster,C11$cluster,C12$cluster,C11$cluster,C12$cluster,C13$cluster,C14$cluster,C15$cluster);pictclsificacionNbclust
save(pictclsificacionNbclust, file= 
       "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/pictclsificacionNbclust.txt")

write.table(pictclsificacionNbclust, "C:/Users/Franca/Desktop/Franca/Doctorado/Analisis de datos/SIG-pict/pictclsificacionNb.txt",sep="\t")