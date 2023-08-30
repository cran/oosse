## ----install, eval = FALSE----------------------------------------------------
#  install.packages("oosse")

## ----load---------------------------------------------------------------------
library(oosse)

## ----loadBrassica-------------------------------------------------------------
data(Brassica)

## ----linModelFit--------------------------------------------------------------
fitFunLM = function(y, x){lm.fit(y = y, x = cbind(1, x))}

## ----linModelPredict----------------------------------------------------------
predFunLM = function(mod, x) {cbind(1,x) %*% mod$coef}

## ----multithread--------------------------------------------------------------
nCores = 2 # For CRAN build max 2
library(BiocParallel)
if(.Platform$OS.type == "unix"){
    #On unix-based systems, use MulticoreParam
    register(MulticoreParam(nCores))
} else {
    #On windows, use makeCluster
    library(doParallel)
    Clus = makeCluster(nCores)
    registerDoParallel(Clus)
    register(DoparParam(), default = TRUE)
}

## ----LMpred-------------------------------------------------------------------
R2lm = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:5], 
               fitFun = fitFunLM, predFun = predFunLM)

## ----lmests-------------------------------------------------------------------
#R2
R2lm$R2
#MSE
R2lm$MSE
#MST
R2lm$MST

## ----confintlm----------------------------------------------------------------
# R2
buildConfInt(R2lm)
#MSE, 90% confidence interval
buildConfInt(R2lm, what = "MSE", conf = 0.9)
#MST, based on chi-square distribution
buildConfInt(R2lm, what = "MST")

## ----lmBoot-------------------------------------------------------------------
R2lm632jn = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:5], 
                    fitFun = fitFunLM, predFun = predFunLM, methodMSE = "bootstrap",
                    methodCor = "jackknife")

## ----reglinmod----------------------------------------------------------------
fitFunReg = function(y, x, ...) {cv.glmnet(y = y, x = x, ...)}
predFunReg = function(mod, x, ...){predict(mod, newx = x, ...)}

## ----reglinmodR2--------------------------------------------------------------
nFolds = 5; cvReps = 1e2; nBoots = 4e1;numFeat = 25
if(require(glmnet)){
    if(onWindows <- (.Platform$OS.type == "windows")){
        clusterEvalQ(Clus, require(glmnet))
    }
    R2pen = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(numFeat)], #Subset genes for speed
                    nFolds = nFolds, cvReps = cvReps, nBootstrapsCor = nBoots,
                        fitFun = fitFunReg, predFun = predFunReg, alpha = 1)#Lasso model
    R2pen$R2
}

## ----predRf-------------------------------------------------------------------
if(require(randomForest)){
   if(onWindows){
        clusterEvalQ(Clus, require(randomForest))
    }
    fitFunrf = function(y, x, ...){randomForest(y = y, x, ...)}
    predFunrf = function(mod, x, ...){predict(mod, x, ...)}
    R2rf = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(numFeat)],
                     nFolds = nFolds, cvReps = cvReps, nBootstrapsCor = nBoots,
                        fitFun = fitFunrf, predFun = predFunrf)
    R2rf$R2
}
if(onWindows){
    stopCluster(Clus)
}

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

