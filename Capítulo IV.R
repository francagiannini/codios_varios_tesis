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
tmedia <- as.data.frame(read.table(file="~/Documents/vida_media.txt", header=T, sep = "\t"))

gbmGrid <-  expand.grid(interaction.depth = c(10:20), 
                        n.trees = (1:40)*100, 
                        shrinkage = c(0.001,0.01),
                        n.minobsinnode = c(10,5,3))
control<- trainControl(method="repeatedcv", number=5,repeats =5)


gbm_optim_tmedia<- train(x=tmedia[,c(4:26)],y=tmedia[,37], 
                         method = "gbm",
                         trControl = control,
                         verbose = FALSE,
                         #importance=T,
                         metric="RMSE",
                         tuneGrid = gbmGrid)

tmedia_brt <- as.data.frame(summary(gbm_optim_tmedia))
colnames(tmedia_brt) <- c("Variable", "BRT")
summary(gbm_optim_tmedia)
plot(gbm_optim_tmedia)

tmedia_brt$results[which.min(tmedia_brt$results$RMSE),]

library(ggregplot) 


INLAselect_tmedia <- INLAModelSel(Response ="LN_vida_media", Explanatory = c("Fe.Ox.", "Kppm","EC","Elevation", "Sand","Cappm","pH", "Nappm",  
                                                                             "Mn", "TvsPP" , "Clay", "Silt" , "Tm","Al.Ox." ,"Zn" , "TN",      
                                                                             "CEC", "Mgppm", "SOC", "WHC" , "pp", "Cu", "Pppm" ), 
                                  Family = "gaussian", Data = tmedia, Delta = 1)

INLAselect_tmedia$DIC

var_sel_tmedia <- as.data.frame(read.table(file="~/Documents/var_sel_vida_media.txt", header=T, sep = "\t"))
colnames(var_sel_tmedia) <- c("Variable", "DIC", "BRT")

ggplot(var_sel_tmedia,aes(reorder(Variable, BRT))) +
  geom_line(aes(y = BRT, group=), group = 1,color="#8B0000") +
  geom_text(aes(y = BRT, 
                label = paste(round(BRT, digits = 1),'%')),
            vjust = 1.4, color = "#8B0000", size = 3)+
  #scale_discrete_manual(aesthetics, ..., values)
  geom_line(aes(y=DIC*0.7-30, group=), group = 1,color="#191970") + 
  geom_text(aes(y = DIC*0.7-30, label = paste(round(DIC, digits = 1))),
            vjust = 1.4, color = "#191970", size = 3)+ 
  scale_y_continuous(sec.axis = sec_axis(~.+40, name = "DIC"))+
  scale_x_discrete("Variables")+
  ggtitle("Atrazine variable selection") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
#"Silt", "Pppm"
#HostModelSel <- INLAModelSel(resp, covar, "ID", "iid", "nbinomial", TestHosts)
#http://localhost:456/graphics/plot_zoom_png?width=940&height=896


##############model fitting####
#formula specification 
formula_spde <-LN_vida_media ~ Kppm+TvsPP+Clay+Clas_graminea+EC + f(node, model = spde, diagonal = 1e-6)#

#spde 0.5
spde <- inla.spde2.matern(mesh=mesh1, alpha=2)


inla_pred_spde_tmed <- inla(formula_spde, family = 'gaussian', data=tmedia
                            ,control.predictor = list(link = 1, compute = TRUE)
                            #,control.inla = list(h=0.0001)
                            ,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE))

### range y sigma2 by Bangliardo Cammeletti et al pag 207

spdeparam <-  inla_pred_spde_atr$summary.hyperpar
range <- sqrt(8/exp(spdeparam$mean[2]))#theta1
sigma2 <- 1/(4*pi*((exp(spdeparam$mean[2]))^2)*((exp(spdeparam$mean[3]))^2))
rsv <- range/sigma2*100

##########predictive model######
inla_tmedia_pred <-function(i){
  
  library(devtools)
  library(INLA)
  library(faraway)
  library(gridExtra)
  library(brinla)
  library(nlme)
  library(geoR)
  library(geoRglm)
  library(sp)
  
  tmedia <- read.table(file="~/Documents/vida_media_usos.txt", header=T, sep = "\t")
  
  #structure predicted table
  tablepred <- matrix(NA, nrow = nrow(tmedia), ncol = 4)
  
  colnames(tablepred)<- c("inlaspde_at_0.5","inlaspde_at_1","inlaspde_at_1.5", "inlaspde_at_2")
  
  train     <- tmedia[-i, ]
  test     <- tmedia[i, ]
  
  iNA <- as.data.frame(cbind(test[,1:39],"LN_vida_media"=NA))
  
  testi <- as.data.frame(rbind(train,iNA))
  
  loc.obs <- cbind(testi$X.UTM20, testi$Y.UTM20)
  
  mesh1 <- inla.mesh.2d(loc.obs, cutoff = 200, max.edge = 200000)
  
  node <- mesh1$idx$loc
  
  #formula specification 
  formula_spde <-LN_vida_media ~ Kppm+TvsPP+Clay+Clas_graminea+EC + f(node, model = spde, diagonal = 1e-6)
  
  #spde 1
  spde <- inla.spde2.matern(mesh=mesh1, alpha=0.5)
  try({
    
    
    inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = tmedia#testi
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.inla = list(h=0.0001)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde_at_0.5"] <- (inla_pred_spde$summary.fitted.values$mean[60])
  })
  #############################  # #spde 1#################
  spde <- inla.spde2.matern(mesh=mesh1, alpha=1)
  
  try({
    
    inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.inla = list(h=0.0001)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde_at_1"] <- (inla_pred_spde$summary.fitted.values$mean[60])
  })
  
  #spde 1.5
  spde <- inla.spde2.matern(mesh=mesh1, alpha=1.5)
  try({
    
    inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.inla = list(h=0.0001)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde_at_1.5"] <- (inla_pred_spde$summary.fitted.values$mean[60])
  })
  
  #spde 2
  spde <- inla.spde2.matern(mesh=mesh1, alpha=2)
  
  try({
    
    inla_pred_spde <- inla(formula_spde, family = 'gaussian', data = testi
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.inla = list(h=0.0001)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspde_at_2"] <- (inla_pred_spde$summary.fitted.values$mean[60])
  })
  
  #####
  try({
    
    tabla=tablepred[!apply(tablepred,1, function(X){all(is.na(X))}),]
    return(tabla)
    return(tablepred)
  })
  
}

results_tmedia <-do.call("rbind",bplapply(1:60, inla_tmedia_pred, BPPARAM=SnowParam(workers=10, progressbar=TRUE, type="SOCK")))



RMSPE=apply(results_tmedia, 2, function (x) {sqrt(mean((x-tmedia$LN_vida_media)^2))})


RMSPE/mean(tmedia$LN_vida_media)*100

#plot(results_tmedia-tmedia$LN_vida_media, tmedia$Clay)


write.table(cbind(tmedia,results_tmedia, SSE), file="~/Documents/tmedia_atr.txt", sep="\t")

SSE <- as.numeric((results_tmedia[,4]-tmedia$LN_vida_media)/tmedia$LN_vida_media*100)

SSEcat <- cut(abs(SSE), breaks=c(-Inf, 25, Inf), labels=c("low","high"))

summary(SSEcat)

SEEplot <- cbind(SEE,SEEcat)
colnames(SEEplot) <- c("SEE", "Xt","Yt", "SEEcat")

############mapping#########
######defining mesh and smothness parameter#####
param <- expand.grid(alpha=c(0.5, 1, 1.5, 2), maxedge=c(100000,10000,2000))
x <- param[1,]

mapinlaoptim <-apply(param, 1, function (x){ 
  
  loc.obs <- cbind(tmedia$X.UTM20, tmedia$Y.UTM20)
  mesh <- inla.mesh.2d(loc.obs, cutoff = 200,
                       max.edge = x["maxedge"] )
  node <- mesh$idx$loc
  spde<- inla.spde2.matern(mesh=mesh1, x["alpha"])
  formula <-LN_vida_media~1+ f(node, model = spde, diagonal = 1e-6)
  inla_spde <- inla(formula, family = 'gaussian', data = tmedia
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
#write.table(b, file="~/Documents/cc.txt", sep = "\t", row.names = T)

spdeparam=read.table(file="~/Documents/cc.txt", header = T)


### range y sigma2 by Bangliardo Cammeletti et al pag 207
range <- sqrt(8/exp(spdeparam$Theta1))
sigma2 <- 1/(4*pi*((exp(spdeparam$Theta1))^2)*((exp(spdeparam$Theta2))^2))
rsv <- range/sigma2

spdeparam=cbind(spdeparam,range,sigma2,rsv)

##selected parameters
meshparam <- spdeparam[which.max(spdeparam$rsv),]


#########mapping prediction ####

tmedia2 <- as.data.frame(read.table(file="~/Documents/vida_media_usos2.txt", header=T, sep = "\t", na.strings = "."))

library(devtools)
library(INLA)
library(faraway)
library(gridExtra)
library(brinla)
library(nlme)
library(geoR)
library(geoRglm)
library(sp)
library(fields)

grida <- cbind(as.data.frame(read.table(file = "~/Documents/pred_usos2.txt", header=TRUE, sep = "\t")))
grid1 <- grida[,c("X","Y")]
colnames(grid1) <- c("Xt","Yt")
grid <- grid1

loc.obs <- cbind(tmedia2$X.UTM20, tmedia2$Y.UTM20)
mesh1 <- inla.mesh.2d(loc.obs, cutoff = 50,
                      max.edge = 90000)

plot(mesh1)

A1 <- inla.spde.make.A(mesh1, loc=loc.obs)

spde <- inla.spde2.matern(mesh=mesh1, alpha=2)

stk <- inla.stack(data=list(resp=tmedia2$LN_vida_media),A=list(A1,1), 
                  effects=list(i=1:spde$n.spde,
                               data.frame(Intercept=1,
                                          Kppm=tmedia2$Kppm,
                                          TvsPP=tmedia2$TvsPP, 
                                          EC=tmedia2$EC, 
                                          pH=tmedia2$pH,
                                          Clay=tmedia2$Clay,
                                          adapt=tmedia2$adapt)),
                  tag='est')

#LN_vida_media ~ Kppm + pH + TvsPP + Clay +EC + Clas_graminea + f(node, model = spde, diagonal = 1e-6)
res <- inla(resp ~ 0 + Intercept+Kppm+TvsPP+Clay+f(i, model=spde),#+ Kppm +Clay + pH +EC + Clas_graminea + TvsPP  +
            data=inla.stack.data(stk),
            control.predictor=list(A=inla.stack.A(stk)))

project <- inla.mesh.projector(mesh1, loc = cbind(grid$Xt, grid$Yt)) 

spa.mean <- inla.mesh.project(project, res$summary.ran$i$mean)

quilt.plot(grid$Xt, grid$Yt, spa.mean, nx = 150, ny = 150, main="spatial (s)",xlab="Longitude", ylab="Latitude", asp=1)


##building stacks
stkgrid <- inla.stack(data=list(resp=NA), A=list(project$proj$A, 1),
                      effects=list(i=1:spde$n.spde,
                                   data.frame(Intercept=1,
                                              Kppm=grida$Kppm,
                                              TvsPP=grida$TvsPP, 
                                              EC=grida$EC, 
                                              pH=grida$pH,
                                              Clay=grida$pH,
                                              adapt=grida$adapt)),
                      
                      tag='prd.grd')

stk.all <- inla.stack(stk, stkgrid)

##res2 with parameters which maxmize rsv

resa <- inla(resp ~ 0 + Intercept+Kppm+TvsPP+Clay+f(i , model=spde),
             data=inla.stack.data(stk.all),
             control.predictor=list(A=inla.stack.A(stk.all),
                                    compute=TRUE),
             control.results=list(return.marginals.random=FALSE,
                                  return.marginals.predictor=FALSE))

igr <- inla.stack.index(stk.all, 'prd.grd')$data

res_tmedia_pred_table <-cbind("Xt"=project$loc[,1],"Yt"=project$loc[,2] , as.data.frame(resa$summary.fitt[igr,1:5]), "ID"=grida$Columna1)

quilt.plot(res_tmedia_pred_table[,1], res_tmedia_pred_table[,2], exp(res_tmedia_pred_table[,3]), nx = 150, ny = 150, main="tmedia",xlab="Longitude", ylab="Latitude", asp=1)

quilt.plot(res_tmedia_pred_table[,1], res_tmedia_pred_table[,2], exp(res_tmedia_pred_table[,4]), nx = 150, ny = 150, main="tmedia SD",xlab="Longitude", ylab="Latitude", asp=1)

quilt.plot(res_tmedia_pred_table[,1][exp(res_tmedia_pred_table[,4])<10]
           , res_tmedia_pred_table[,2][exp(res_tmedia_pred_table[,4])<10]
           , exp(res_tmedia_pred_table[,4])[exp(res_tmedia_pred_table[,4])<10], nx = 150, ny = 150, main="tmedia SD",xlab="Longitude", ylab="Latitude", asp=1)

quilt.plot(res_tmedia_pred_table[,1], res_tmedia_pred_table[,2], (exp(res_tmedia_pred_table[,7])-exp(res_tmedia_pred_table[,5])), nx = 150, ny = 150, main="IC range tmedia",xlab="Longitude", ylab="Latitude", asp=1)

write.table(res_tmedia_pred_table, file="~/Documents/res_tmedia_pred_table.txt", sep="\t")
