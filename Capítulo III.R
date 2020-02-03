#Kd atrazina

library(INLA)
library(devtools)
library(ggplot2)
library(faraway)
library(gridExtra)
library(brinla)
library(nlme)
library(geoR)
library(sp)
library(mapview)
library(raster)
library(survival)
library(splines)
library(parallel)
library(lattice)
library(gbm)
library(dismo)
library(caret)
library(BiocParallel)
library(doMC)
registerDoMC(10)
##########variable selection######
atr <- read.table(file="~/Documents/atr.txt", header=T, sep = "\t")

gbmGrid <-  expand.grid(interaction.depth = c(10:20), 
                        n.trees = (1:30)*100, 
                        shrinkage = c(0.001,0.01, 0.1),
                        n.minobsinnode = c(15,10,5))
control<- trainControl(method="repeatedcv", number=5,repeats =5)


gbm_optim_Kda<- train(x=atr[,c(4:24)],y=atr[,26], 
                      method = "gbm",
                      trControl = control,
                      verbose = FALSE,
                      #importance=T,
                      metric="RMSE",
                      tuneGrid = gbmGrid)

atrbrt <- as.data.frame(summary(gbm_optim_Kda))
colnames(atrbrt) <- c("Variable", "BRT")
plot(gbm_optim_Kda)

gbm_optim_Kda$results[which.min(gbm_optim_Kda$results$RMSE),]

library(ggregplot) 


INLAselect2 <- INLAModelSel(Response ="LN_Kda", Explanatory = c("TvsPP", "PPanual", "SOC", "CEC", "WHC", "Clay", "TN", 
                                                                "Elevation", "Zn", "Nappm", "Cappm", "Mn", "Sand", "pH", "Kppm",
                                                                "Mgppm", "Tm", "EC", "Cu","Silt","Pppm"), 
                            Family = "gaussian", Data = atr, Delta = 1)

atrdic <- as.data.frame(read.table(file="~/Documents/var_sel_atr.txt", header=T, sep = "\t"))

Variables <- as.data.frame(c("Ca", "CEC", "Clay", "Cu",  "EC",  "Elevation", "K",  "Mg", "Mn", "Na", "pH", "pp", 
                             "P", "Sand", "Silt", "SOC", "Tm", "TN", "TvsPP", "WHC", "Zn"))
colnames(Variables) <- "Variables"

var_sel_atr <-cbind(Variables,merge(atrbrt, atrdic))

colnames(var_sel_atr) <- c("Variables","Variable", "BRT","DIC")


ggplotBRT <- ggplot(var_sel_atr,aes(reorder(Variables, BRT))) +
  geom_bar(aes(y = BRT),stat="identity", group = 1, fill="#8B0000", width = 0.7) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ggplotDIC <-   ggplot(var_sel_atr,aes(reorder(Variables, BRT))) +
  geom_bar(aes(y=DIC), stat="identity", group = 1,fill="#191970", width = 0.7) + 
  scale_x_discrete("Variables")+
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))


library(ggpubr)

ggarrange(ggplotBRT, ggplotDIC, ncol = 1)




dfggplot <- stack(var_sel_atr)
ggplot(var_sel_atr,aes(reorder(Variable, BRT))) +
  geom_bar(aes(y = BRT),stat="identity", group = 1, fill="#8B0000") 

ggplot(var_sel_atr,aes(reorder(Variable, BRT))) +
  geom_bar(aes(y=DIC), stat="identity", group = 1,fill="#191970")



#"Silt", "Pppm"
#HostModelSel <- INLAModelSel(resp, covar, "ID", "iid", "nbinomial", TestHosts)
#http://localhost:456/graphics/plot_zoom_png?width=940&height=896


##############model fitting#####
atr <- read.table(file="~/Documents/vida_media_usos.txt", header=T, sep = "\t", na.strings = ".")
kda <- as.data.frame(cbind(atr,"SOC2"=(atr[,"SOC"])^2)) #, LN_Kda2=atr[,26]))        
#formula specification 
formula_spde <-LN_Kda ~ TvsPP +pp+ SOC +SOC2 + Clay + f(node, model = spde, diagonal = 1e-6)#

loc.obs <- cbind(kda$X.UTM20, kda$Y.UTM20)

mesh1 <- inla.mesh.2d(loc.obs, cutoff = 50,
                      max.edge = 200000)

#spde 0.5
spde <- inla.spde2.matern(mesh=mesh1, alpha=1)


inla_pred_spde_atr <- inla(formula_spde, family = 'gaussian', data=kda
                           ,control.predictor = list(link = 1, compute = TRUE)
                           ,control.inla = list(h=0.0001)
                           ,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE))

write.table(cbind(inla_pred_spde_atr$summary.fitted.values,kda$ID_2), file = "~/Documents/pred60LN_kda.txt")
### range y sigma2 by Bangliardo Cammeletti et al pag 207

spdeparam <-  inla_pred_spde_atr$summary.hyperpar
range <- sqrt(8/exp(spdeparam$mean[2]))#theta1
sigma2 <- 1/(4*pi*((exp(spdeparam$mean[2]))^2)*((exp(spdeparam$mean[3]))^2))
rsv <- range/sigma2*100

##########predictive model######
inlakdpred <-function(i){
  
  library(devtools)
  library(INLA)
  library(faraway)
  library(gridExtra)
  library(brinla)
  library(nlme)
  library(geoR)
  library(geoRglm)
  library(sp)
  
  atr <- read.table(file="~/Documents/atr.txt", header=T, sep = "\t")
  #kda <- cbind(atr[,-27])
  kda <- as.data.frame(cbind(atr,"SOC2"=(atr[,"SOC"])^2, LN_Kda2=atr[,26]))
  colnames(kda)
  
  #structure residual table
  tablepred <- matrix(NA, nrow = nrow(kda), ncol = 1)
  
  colnames(tablepred)<- c("inlaspde_at_1")#,"inlaspde_at_1.005")#,
  #"inlaspde_at_0.97", "inlaspde_at_0.98")
  
  train     <- kda[-i, ]
  test     <- kda[i, ]
  
  iNA <- as.data.frame(cbind(test[,1:28],"LN_Kda2"=NA))
  
  testi <- as.data.frame(rbind(train,iNA))
  
  loc.obs <- cbind(testi$Xt, testi$Yt)
  
  mesh1 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 200000)
  node <- mesh1$idx$loc
  
  #formula specification 
  formula_spde <-LN_Kda2 ~ TvsPP +PPanual+ SOC +SOC2 + Clay + f(node, model = spde, diagonal = 1e-6)#
  
  #spde 1
  spde <- inla.spde2.matern(mesh=mesh1, alpha=1)
  try({
    
    
    inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.inla = list(h=0.0001)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde_at_1"] <- (inla_pred_spde$summary.fitted.values$mean[154])
  })
  #############################  # #spde 1#################
  # spde <- inla.spde2.matern(mesh=mesh1, alpha=1.005)
  # 
  # try({
  #   
  #   inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
  #                          ,control.predictor = list(link = 1, compute = TRUE)
  #                          #,control.inla = list(h=0.0001)
  #                          #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
  #   )
  #   
  #   tablepred[i, "inlaspde_at_1.005"] <- (inla_pred_spde$summary.fitted.values$mean[154])
  # })
  # 
  # #spde 1.5
  # spde <- inla.spde2.matern(mesh=mesh1, alpha=0.98)
  # try({
  # 
  #   inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
  #                          ,control.predictor = list(link = 1, compute = TRUE)
  #                          #,control.inla = list(h=0.0001)
  #                          #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
  #   )
  #   
  #   tablepred[i, "inlaspde_at_0.98"] <- (inla_pred_spde$summary.fitted.values$mean[154])
  # })
  # 
  # #spde 2
  # spde <- inla.spde2.matern(mesh=mesh1, alpha=1.15)
  # 
  # try({
  # 
  #   inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
  #                          ,control.predictor = list(link = 1, compute = TRUE)
  #                          #,control.inla = list(h=0.0001)
  #                          #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
  #   )
  #   
  #   tablepred[i, "inlaspde_at_1.15"] <- (inla_pred_spde$summary.fitted.values$mean[154])
  # })
  
  #####
  try({
    
    tabla=tablepred[!apply(tablepred,1, function(X){all(is.na(X))}),]
    return(tabla)
    return(tablepred)
  })
  
}

results_inla_at1 <-do.call("rbind",bplapply(1:154, inlakdpred, BPPARAM=SnowParam(workers=10, progressbar=TRUE, type="SOCK")))


RMSPE=apply(results_inla_at1, 2, function (x) {sqrt(mean((x-kda$LN_Kda)^2))})

write.table(cbind(atr,results_inla_at1), file="~/Documents/pred152_atr.txt", sep="\t")

############mapping#########

##defining mesh and smothness parameter
param <- expand.grid(alpha=c(0.5, 1, 1.5, 2), maxedge=c(100000,10000,2000))
x <- param[1,]

mapinlaoptim <-apply(param, 1, function (x){ 
  
  loc.obs <- cbind(kda$Xt, kda$Yt)
  mesh <- inla.mesh.2d(loc.obs, cutoff = 200,
                       max.edge = x["maxedge"] )
  node <- mesh$idx$loc
  spde<- inla.spde2.matern(mesh=mesh1, x["alpha"])
  formula <-LN_Kda~1+ f(node, model = spde, diagonal = 1e-6)
  inla_spde <- inla(formula, family = 'gaussian', data = kda
                    ,control.predictor = list(link = 1, compute = TRUE)
                    #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
  )
  
  list(x,"hyperpar"=inla_spde$summary.hyperpar)
  
}
)
####mesh
#spdeparam <- read.table(file = "clipboard", sep = "\t", header=TRUE)

a <- as.data.frame(do.call(rbind, mapinlaoptim))
b <- do.call(rbind, a[,2])
names<- do.call(rbind, a[,1])
write.table(b, file="~/Documents/c.txt", sep = "\t", row.names = T)


spdeparam=read.table(file="~/Documents/spde_param_mesh_atr.txt", header = T)

### range y sigma2 by Bangliardo Cammeletti et al pag 207
range <- sqrt(8/exp(spdeparam$Theta1))
sigma2 <- 1/(4*pi*((exp(spdeparam$Theta1))^2)*((exp(spdeparam$Theta2))^2))
rsv <- range/sigma2

spdeparam=cbind(spdeparam,range,sigma2,rsv)

##selected parameters
meshparam <- spdeparam[which.max(spdeparam$rsv),]

## Kd Glifosato
gly <- read.table(file="~/Documents/gly.txt", header=T, sep = "\t")

gbmGrid <-  expand.grid(interaction.depth = c(10:20), 
                        n.trees = (1:30)*100, 
                        shrinkage = c(0.001,0.01, 0.1),
                        n.minobsinnode = c(10,5,3))

control<- trainControl(method="repeatedcv", number=5,repeats =5)


gbm_optim_gly<- train(x=gly[,c(4:26)],y=gly[,28], 
                      method = "gbm",
                      trControl = control,
                      verbose = FALSE,
                      #importance=T,
                      metric="RMSE",
                      tuneGrid = gbmGrid)

glybrt <- as.data.frame(summary(gbm_optim_gly))
colnames(glybrt) <- c("Variable", "BRT")

plot(gbm_optim_gly)

gbm_optim_gly$results[which.min(gbm_optim_gly$results$RMSE),]

library(ggregplot) 


INLAselectgly <- INLAModelSel(Response ="LN_Kdg", Explanatory = c("X.Al","Clay","Sand","PPanual","pH","X..Fe",
                                                                  "WHC","Cappm","SOC","EC","Mn","Mgppm","CEC","TvsPP",
                                                                  "Nappm","Elevation","Silt","Pppm","Kppm",
                                                                  "Tm","Zn","TN","Cu"), Family = "gaussian", Data = gly)


glydic <- as.data.frame(read.table(file="~/Documents/var_sel_gly.txt", header=T, sep = "\t"))

var_sel_gly <- merge(glybrt,glydic)

colnames(var_sel_gly) <- c("Variable", "BRT","DIC")

Variables <- as.data.frame(c( "Ca", "CEC","Clay", "Cu", "EC",  "Elevation",  "K", "Mg", "Mn",  "Na",    
                              "pH", "pp",  "P",    "Sand","Silt", "SOC", "Tm", "TN",  "TvsPP", "WHC", 
                              "Fe(Ox)", "Al(Ox)", "Zn"))
colnames(Variables) <- "Variables"

var_sel_gly <- cbind(Variables, var_sel_gly)



ggplotBRT <- ggplot(var_sel_gly,aes(reorder(Variables, BRT))) +
  geom_bar(aes(y = BRT),stat="identity", group = 1, fill="#8B0000", width = 0.7)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ggplotDIC <-   ggplot(var_sel_gly,aes(reorder(Variables, BRT))) +
  geom_bar(aes(y=DIC), stat="identity", group = 1,fill="#191970", width = 0.7) + 
  scale_x_discrete("Variables")+
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

library(ggpubr)

ggarrange(ggplotBRT, ggplotDIC, ncol = 1)

#https://community.rstudio.com/t/assign-2-geom-lines-to-different-y-axis/14797

ggplot(var_sel_gly,aes(reorder(Variable, BRT))) + 
  geom_line(aes(y = BRT, group=), group = 1,color="#8B0000") + 
  geom_text(aes(y = BRT, 
                label = paste(round(BRT, digits = 1),'%')),
            vjust = 1.4, color = "#8B0000", size = 3)+
  #scale_discrete_manual(aesthetics, ..., values)
  geom_line(aes(y=DIC*0.2-25, group=), group = 1,color="#191970") + 
  geom_text(aes(y = DIC*0.2-25, label = paste(round(DIC, digits = 1))),
            vjust = 1.4, color = "#191970", size = 3)+ 
  scale_y_continuous(sec.axis = sec_axis(~., name = "DIC"))+
  scale_x_discrete("Variables")+
  ggtitle("Glyphosate variable selection") +
  scale_colour_manual(name="Selection criteria",values=cols)


cols <- c("BRT"="#8B0000","DIC"="#191970")
ggplot(data=var_sel_gly,aes(x=reorder(Variable, BRT))) + 
  geom_line(aes(y=BRT,group=1, colour="BRT"),size=1.0) +   #red
  # geom_point(aes(y=BRT, colour="BRT"),size=3) + 
  geom_text(aes(y = BRT, 
                label = paste(round(BRT, digits = 1),'%')),
            vjust = 1.4, color = "#8B0000", size = 3)+#red
  geom_line(aes(y=DIC*0.2-25,group=1,colour="DIC"),size=1.0) +   #blue 
  # geom_point(aes(y=DIC*0.7-30,colour="DIC"),size=3) +  #blue
  geom_text(aes(y = DIC*0.2-25, label = paste(round(DIC, digits = 1))),
            vjust = 1.4, color = "#191970", size = 3)+
  scale_colour_manual(name="Selection criteria",values=cols) +
  ggtitle("Glyphosate variable selection") 
# theme_bw() +
# theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
# theme(axis.title.y = element_text(size = 15, vjust=0.3)) +

#######fitting
kdg <-read.table(file = "~/Documents/gly.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")

kdg <- as.data.frame(cbind(kdg, "X.Al2"=kdg[,"X.Al"]^2,"LN_Kdg"=kdg[,"LN_Kdg"]))

###building mesh
loc.obs <- cbind(testi$Xt, testi$Yt)

mesh2 <- inla.mesh.2d(loc.obs, cutoff = 200,
                      max.edge = 200000)

node <- mesh2$idx$loc

#spde model
spde1085 <- inla.spde2.matern(mesh=mesh2, alpha=1.085)

#define formula
formula_spde1085 <-LN_Kdg ~Clay+X.Al+X.Al2+pH+ f(node, model = spde1085, diagonal = 1e-6)
#inla model
inla_pred_spde1085 <- inla(formula_spde1085, family = 'gaussian', data = kdg
                           ,control.predictor = list(link = 1, compute = TRUE)
                           ,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE))

######
inlakdpredgly <-function(i){
  
  library(devtools)
  library(INLA)
  library(faraway)
  library(gridExtra)
  library(brinla)
  library(nlme)
  library(geoR)
  library(geoRglm)
  library(sp)
  #Loading data
  
  #setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")
  kdg <-read.table(file = "~/Documents/gly.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")
  
  kdg <- as.data.frame(cbind(kdg, "X.Al2"=kdg[,"X.Al"]^2,"LN_Kdg"=kdg[,"LN_Kdg"]))
  
  head(kdg)
  
  #structure residual table
  tablepred <- matrix(NA, nrow = nrow(kdg), ncol = 1)
  
  colnames(tablepred)<- c("inlaspde")
  #test-train set
  train     <- kdg[-i, ]
  train1    <- kdg[-i, ]
  test     <- kdg[i, ]
  test1     <- kdg[i, ]
  #defining missing value
  iNA <- as.data.frame(cbind(test[,c(1:29)],"LN_Kdg"=NA))
  
  testi <- as.data.frame(rbind(train,iNA))
  
  ###building mesh
  loc.obs <- cbind(testi$Xt, testi$Yt)
  
  mesh2 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 200000)
  
  node <- mesh2$idx$loc
  
  #spde model
  spde1085 <- inla.spde2.matern(mesh=mesh2, alpha=1.085)
  try({
    #define formula
    formula_spde1085 <-LN_Kdg ~Clay+X.Al+X.Al2+pH+ f(node, model = spde1085, diagonal = 1e-6)
    #inla model
    inla_pred_spde1085 <- inla(formula_spde1085, family = 'gaussian', data = testi#data=kdg
                               ,control.predictor = list(link = 1, compute = TRUE)
                               #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde"] <- (inla_pred_spde1085$summary.fitted.values$mean[89])
  })
  
  try({
    #result table
    tabla=tablepred[!apply(tablepred,1, function(X){all(is.na(X))}),]
    return(tabla)
    return(tablepred)
  })
  
}

#predictive inla values
results_inla_gly<-do.call("rbind",bplapply(1:89, inlakdpredgly, BPPARAM=SnowParam(workers=2, progressbar=TRUE, type="SOCK")))
#Prediction errors RMSPE
apply(results_inla_gly, 2 ,function (x) {(sqrt(mean((x-kdg$LN_Kdg)^2)))})


write.table(cbind(gly,results_inla_gly), file="~/Documents/pred89_gly.txt", sep="\t")


