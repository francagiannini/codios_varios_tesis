#Comparación de predicciones obtenidas desde distintos marcos teóricos ML-REML, PLRS, GBR, RFR

library(BiocParallel) 

i=1
#k=1

prediccion <-function(i){
  library(gstat)
  library(sp)
  library(automap)
  library(spdep)
  library(randomForest)
  library(gbm)
  library(dismo)
  library(spm)
  library(MASS)
  library(pls)
  library(nlme)
  library(caret) 
  setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")
  
  kdg <-read.table(file = "kdg_selected_pff_cen.txt", header=TRUE, sep= "\t",  na.strings = ".")
  #kdg <- kdg[,c("Yt","Xt","X.Al","ARC","X.Fe","CO","Mn","CC", "pH",  "ARENA", "LN_KdGlifo")]
  #dim(kdg)
  head(kdg)
  
  kdg1 <-kdg
  class(kdg1)
  
  #leave one out cross validation
  #generate structure
  tablaresiduos <- matrix(NA, nrow = nrow(kdg), ncol = 9)
  colnames(tablaresiduos)<- c("REML","LM" ,"RK","GBM", "GBMR", "PLS" ,"PLSR","RF","RFR")
  
  #set train and validation set
  train     <- kdg[-i, ]
  train1    <- kdg1[-i, ] 
  test     <- kdg[i, ]
  test1     <- kdg1[i, ]
  
  coordinates(train)<-~Xt+Yt
  coordinates(test)<-~Xt+Yt
  
  ##########    REML    ##############
  try({
    
    REML <-gls(LN_KdGlifo ~1+ CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_Mn+CENc_pH+CENc_POT_pH+Xt+Yt
               ,correlation=corExp(form=~as.numeric(as.character(Xt))+as.numeric(as.character(Yt))
                                   ,metric="euclidean",nugget=FALSE) ,method="REML",na.action=na.omit,data=train)
    
    pred_REML = predict(REML, newdata = test)
    tablaresiduos[i, "REML"] <- (test$LN_KdGlifo - pred_REML)
  })
  ####### Regression Kriging  ######
  try({
    RK<- lm(LN_KdGlifo ~1+ CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_Mn+CENc_pH+CENc_POT_pH+Xt+Yt, data=train)
    VF_VUT_KED <- autofitVariogram(RK$residuals~1,train, c("Sph", "Exp", "Gau", "Mat", "Ste"),cressie=T)
    
    RK_pred <- predict.lm(RK, test)
    tablaresiduos[i, "LM"] <- (test$LN_KdGlifo - RK_pred)
    
    RK_res <- krige(residuals(RK)~1, train, newdata=test, model=VF_VUT_KED$var_model)
    
    RK_pred_reg <- RK_pred + RK_res$var1.pred
    tablaresiduos[i, "RK"] <- (test$LN_KdGlifo - RK_pred_reg)
  })
  
  ############# PLS - Regression ######################
  try({
    PLS <- plsr(formula=LN_KdGlifo ~1+ CENc_.Al+CENc_ARENA+CENc_Mn+CENc_pH+Xt+Yt,  data = train1, validation = "LO",jackknife = TRUE)
    predPLS <- predict(PLS, ncomp = 6, newdata = test1)
    residuoPLS = test$LN_KdGlifo - predPLS
    
    predPLS_train <- predict(PLS, ncomp = 6, newdata = train1)
    resPLS <- train$LN_KdGlifo - predPLS_train
    names(resPLS) <- "resPLS"
    VF_PLS_res <- autofitVariogram(resPLS~1,train, c("Sph", "Exp", "Gau", "Mat", "Ste"),cressie=T)
    RKPLS_pred <- krige(resPLS~1, train, newdata=test, model=VF_PLS_res$var_model,nmax=25)
    
    predPLS_test <- predict(PLS, ncomp = 6, newdata = test1)
    tablaresiduos[i, "PLS"] <- (test1$LN_KdGlifo - predPLS_test)
    
    RKPLS_pred<- predPLS_test + RKPLS_pred$var1.pred
    
    tablaresiduos[i, "PLSR"] <- (test1$LN_KdGlifo - RKPLS_pred)
  })
  ############# GBM- Regression ######################
  try({
    GBM <- gbm(formula=LN_KdGlifo ~1+ CENc_.Al+CENc_ARENA+CENc_Mn+CENc_pH+Xt+Yt,
               distribution = "gaussian",data=train1,
               n.trees = 550,interaction.depth = 2, n.minobsinnode = 4, shrinkage = 0.01,
               bag.fraction = 0.5, train.fraction = 1, cv.folds = 0,
               keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL, n.cores = 2)
    
    GBMpred_train <- predict(GBM, newdata=train1, n.trees=550, type="response")
    
    resGBM <-GBMpred_train -train$LN_KdGlifo
    train$GBMres <-resGBM
    vgm_GBM <- autofitVariogram(GBMres~1,train, c("Sph", "Exp", "Gau", "Mat", "Ste"),cressie=T)
    
    pred_res_ok_gbm <- krige(GBMres~1, locations=train, newdata=test, model=vgm_GBM$var_model)
    
    predGBM_test <- predict(GBM, newdata=test1, n.trees=550, type="response")
    tablaresiduos[i, "GBM"] <- test$LN_KdGlifo - predGBM_test
    
    RKGBM_pred <- predGBM_test + pred_res_ok_gbm$var1.pred
    
    tablaresiduos[i, "GBMR"] <- test$LN_KdGlifo - RKGBM_pred
    
  })
  ############# RF- Regression ######################
  try({
    rftrain <- train(LN_KdGlifo ~1+ CENc_.Al+CENc_ARENA+CENc_Mn+CENc_pH+Xt+Yt, method="rf", data=train1)
    
    RF <- randomForest(formula=LN_KdGlifo ~1+ CENc_.Al+CENc_ARENA+CENc_Mn+CENc_pH+Xt+Yt, mtry=2 ,data=train1)
    
    RFpred_train <- predict(RF, newdata=train1)
    
    resRF <- RFpred_train -train$LN_KdGlifo
    train$RFres <-resRF
    vgm_RF <- autofitVariogram(RFres~1,train, c("Sph", "Exp", "Gau", "Mat", "Ste"),cressie=T)
    
    pred_res_ok_rf <- krige(RFres~1, locations=train, newdata=test, model=vgm_RF$var_model)
    
    predRF_test <- predict(RF, newdata=test1, n.trees=550, type="response")
    tablaresiduos[i, "RF"] <- test$LN_KdGlifo - predRF_test
    
    RFGBM_pred <- predRF_test + pred_res_ok_gbm$var1.pred
    
    tablaresiduos[i, "RFR"] <- test$LN_KdGlifo - RKGBM_pred
    
  })
  ######results######
  try({
    
    tabla=na.omit(data.frame(tablaresiduos))
    tabla=tablaresiduos[!apply(tablaresiduos,1,function(X){all(is.na(X))}),]
    return(tabla)
    return(tablaresiduos)
  })
  
}

####function call 
residuos_modelos_con_coor<-do.call("rbind",bplapply(1:89, prediccion, BPPARAM=SnowParam(workers=2, progressbar=TRUE, type="SOCK")))
apply(residuos_modelos_con_coor, 2 ,function (x) {(sqrt(mean(x^2)))/ mean(kdg$LN_KdGlifo)*100})


##Regresión bayesiana INLA-SPDE

library(BiocParallel)
i=1

inlakdpred <-function(i){
  
  library(devtools)
  library(INLA)
  library(faraway)
  library(gridExtra)
  library(brinla)
  library(sp)
  #Loading data
  
  setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")
  kdg <-read.table(file = "kdg_selected_pff_cen_inla.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")
  
  #structure residual table
  tablepred <- matrix(NA, nrow = nrow(kdg), ncol = 3)
  
  colnames(tablepred)<- c("inlaspdem1","inlaspdem2", "inlaspdem3")
  
  train     <- kdg[-i, ]
  train1    <- kdg[-i, ]
  test     <- kdg[i, ]
  test1     <- kdg[i, ]
  
  iNA <- as.data.frame(cbind(test[,c(1:12)],"LN_KdGlifo"=NA))
  
  testi <- as.data.frame(rbind(train,iNA))
  
  loc.obs <- cbind(testi$Xt, testi$Yt)
  
  mesh1 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 100000)
  mesh2 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 10000)
  mesh3 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 1000)
  
  node <- mesh1$idx$loc
  
  #spde m1
  spde1085 <- inla.spde2.matern(mesh=mesh1, alpha=1.085)
  try({
    formula_spde1085 <-LN_KdGlifo ~1+ CENc_.Al+CENc_pH+CENc_ARENA+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH #+ f(node, model = spde1085, diagonal = 1e-6)
    
    inla_pred_spde1085 <- inla(formula_spde1085, family = 'gaussian', data = kdg
                               ,control.predictor = list(link = 1, compute = TRUE)
                               #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspdem1"] <- (inla_pred_spde1085$summary.fitted.values$mean[89])
  })
  #spde m2
  spde1085 <- inla.spde2.matern(mesh=mesh2, alpha=1.085)
  try({
    formula_spde1085sc <-LN_KdGlifo ~1+ CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH+f(node, model = spde1085, diagonal = 1e-6)
    
    inla_pred_spde1085sc <- inla(formula_spde1085sc, family = 'gaussian', data = testi
                                 ,control.predictor = list(link = 1, compute = TRUE)
                                 #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspdem2"] <- (inla_pred_spde1085sc$summary.fitted.values$mean[89])
  })
  #spde m3 
  spde1085 <- inla.spde2.matern(mesh=mesh3, alpha=1.085)
  try({
    formulacc <-LN_KdGlifo ~1+ CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH+f(node, model = spde1085, diagonal = 1e-6)
    
    inla_cc <- inla(formulacc, family = 'gaussian', data = testi
                    ,control.predictor = list(link = 1, compute = TRUE)
                    #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )
    
    tablepred[i, "inlaspdem3"] <- (inla_cc$summary.fitted.values$mean[89])
  })
  try({
    
    tabla=tablepred[!apply(tablepred,1, function(X){all(is.na(X))}),]
    return(tabla)
    return(tablepred)
  })
  
}

results_inla<-do.call("rbind",bplapply(11:50, inlakdpred, BPPARAM=SnowParam(workers=2, progressbar=TRUE, type="SOCK")))

results_inla_fran1a10=results_inla_fran
apply(results_inla_fran, 2 ,function (x) {(sqrt(mean((x-kdg$LN_KdGlifo)^2)))/ mean(kdg$LN_KdGlifo)*100})