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
library(BiocParallel)
nCores = 2 # For CRAN build
register(MulticoreParam(nCores))

## ----LMpred-------------------------------------------------------------------
R2lm = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:10], 
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
R2lm632jn = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:10], 
                    fitFun = fitFunLM, predFun = predFunLM, methodMSE = "bootstrap",
                    methodCor = "jackknife")

## ----reglinmod----------------------------------------------------------------
fitFunReg = function(y, x, ...) {cv.glmnet(y = y, x = x, ...)}
predFunReg = function(mod, x, ...){predict(mod, newx = x, ...)}

## ----reglinmodR2--------------------------------------------------------------
if(require(glmnet)){
    R2pen = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(5e1)], #Subset genes for speed
                    nFolds = 4, cvReps = 1e2, nBootstrapsCor = 30,
                        fitFun = fitFunReg, predFun = predFunReg, alpha = 1)#Lasso model
    R2pen$R2
}

## ----predBart-----------------------------------------------------------------
if(require(randomForest)){
    fitFunrf = function(y, x, ...){randomForest(y = y, x, ...)}
    predFunrf = function(mod, x, ...){predict(mod, x, ...)}
    R2rf = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(5e1)],
                     nFolds = 4, cvReps = 1e2, nBootstrapsCor = 30,
                        fitFun = fitFunrf, predFun = predFunrf)
    R2rf$R2
}

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

