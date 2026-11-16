## ----install, eval = FALSE----------------------------------------------------
# install.packages("oosse")

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
param <- if (.Platform$OS.type == "unix") {
    # On unix-based systems (linux and macOS), use MulticoreParam
    MulticoreParam(nCores, stop.on.error = FALSE)
} else {
    # On windows, build vignette without multithreading
    SerialParam(stop.on.error = FALSE)
}
register(param)

## ----LMpred-------------------------------------------------------------------
nFolds = 5; cvReps = 1e2; nBoots = 3e1;cvRepsCor = 5L
R2lm = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:5], 
               fitFun = fitFunLM, predFun = predFunLM, nFolds = nFolds, 
               cvReps = cvReps, nBootstrapsCor = nBoots, cvRepsCor = cvRepsCor)

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

## ----lmBoot, eval = FALSE-----------------------------------------------------
# R2lm632jn = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, 1:5],
#                     fitFun = fitFunLM, predFun = predFunLM, methodMSE = "bootstrap",
#                     methodCor = "jackknife", nFolds = nFolds, cvReps = cvReps,
#                     nBootstrapsCor = nBoots, cvRepsCor = cvRepsCor)

## ----modelMatrix--------------------------------------------------------------
#We construct a fake data frame also containing genotypes
fakeDf = data.frame(Brassica$Expr[, 1:5], "genotype" = sample(c("Genotype1", "Genotype2", "Genotype3"), replace = TRUE, nrow(Brassica$Expr)))
#Build the design matrix. The model.matrix variables automatically constructs dummy variables
designMatrix = model.matrix(~ .  , data = fakeDf)[, -1] #Include no intercept as fitting function already does this
#Now run oosse
R2modMat = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = designMatrix, 
                   nFolds = nFolds, cvReps = cvReps, nBootstrapsCor = nBoots, 
               cvRepsCor = cvRepsCor, fitFun = fitFunLM, predFun = predFunLM)

## ----reglinmod----------------------------------------------------------------
fitFunReg = function(y, x, ...) {cv.glmnet(y = y, x = x, ...)}
predFunReg = function(mod, x, ...){predict(mod, newx = x, ...)}

## ----reglinmodR2--------------------------------------------------------------
numFeat = 20;
if(require(glmnet)){
    R2pen = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(numFeat)], #Subset genes for speed
                    nFolds = nFolds, cvReps = cvReps, nBootstrapsCor = nBoots, cvRepsCor = cvRepsCor,
                        fitFun = fitFunReg, predFun = predFunReg, alpha = 1)#Lasso model
    R2pen$R2
}

## ----predRf-------------------------------------------------------------------
if(require(randomForest)){
    fitFunrf = function(y, x, ...){randomForest(y = y, x, ...)}
    predFunrf = function(mod, x, ...){predict(mod, x, ...)}
    R2rf = R2oosse(y = Brassica$Pheno$Leaf_8_width, x = Brassica$Expr[, seq_len(numFeat)],
                     nFolds = nFolds, cvReps = cvReps, nBootstrapsCor = nBoots,
                   cvRepsCor = cvRepsCor, fitFun = fitFunrf, predFun = predFunrf)
    R2rf$R2
}

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

