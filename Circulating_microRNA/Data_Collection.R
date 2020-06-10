
setwd('~/bigdata/EarlyCancerDetection/miRNA/')

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)


getPhenoFun <- function(phenoData) {
  
  info <- c('geo_accession','title','source_name_ch1','organism_ch1', 'platform_id',
            'contact_institute','contact_name','contact_email')
  
  idx1 <- match(info, colnames(phenoData))
  idx2 <- grep(':ch1', colnames(phenoData))
  idx <- c(idx1, idx2)
  
  phenoData <- phenoData[,idx]
  colnames(phenoData) <- gsub(':ch1', '', colnames(phenoData))
  colnames(phenoData) <- gsub(' ', '.', colnames(phenoData))
  colnames(phenoData)[1:8] <- c('Accession','Title','Source','Organism','Platform','Contact.Institute','Contact.Name','Contact.Email')
  
  return(phenoData)
  
}


#########################################################################
############ GSE112264

gse <- 'GSE112264'

seriesMatrix <- getGEO(gse, AnnotGPL = FALSE, getGPL = TRUE, GSEMatrix = TRUE, destdir = 'data/fromGEO/') # AnnotGPL = TRUE

length(seriesMatrix)

seriesMatrix <- seriesMatrix[[1]]


### Phenotype
phenoData <- pData(seriesMatrix)
phenoData <- getPhenoFun(phenoData)


### Annotation
annoData <- seriesMatrix@featureData@data
dim(annoData)

platform <- seriesMatrix@annotation
platform

saveRDS(phenoData, file=paste0('data/rData/', gse, '_', platform, '_Sample_Information.RDS'))


### Expression
# from series matrix
exprData <- exprs(seriesMatrix)
exprData[1:5,1:5]

rownames(phenoData) == colnames(exprData)

# raw
filePaths = getGEOSuppFiles(gse, baseDir = 'data/fromGEO', makeDirectory = FALSE, filter_regex = 'RAW')
untar(paste0('data/fromGEO/', gse, '_RAW.tar'), exdir = paste0('data/fromGEO/', gse, '_RAW'))

celFiles = list.celfiles(paste0('data/fromGEO/', gse, '_RAW'), full.names=T, listGzipped=T)
celFiles

rawData = read.celfiles(celFiles, pkgname = 'pd.huex10st.hs.gencodeg')

probesetData = oligo::rma(rawData)

exprData = exprs(probesetData)
colnames(exprData)
rownames(exprData)
colnames(exprData) <- unlist(lapply(colnames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
rownames(exprData) <- unlist(lapply(rownames(exprData), function(x) strsplit(x, '.', fixed=T)[[1]][1]))
colnames(exprData) <- gsub('_GBX', '', colnames(exprData))

saveRDS(exprData, file=paste0('data/rData/', gse, '_Expression_CDF24_GENCODE_RMA.RDS'))


#########################################################################


library(caret)

# We set a seed for the sake of reproducibility
set.seed(777)

# First, we'll pick off the training data: 
inTrain<-createDataPartition(y=phenoData$disease.state, p = 0.60, list=FALSE)
inTrain

geno <- data.frame(t(exprData[1:20,]),Type=phenoData$disease.state)

training<-geno[inTrain,]

# Then what's left is testing data.
testing<-geno[-inTrain,]
View(geno)

set.seed(777)
rf.fit <- train(Type ~ ., data=training, method="rf") # svmLinear
confusionMatrix(predict(rf.fit, training), reference=training$Type) # positive='Prostate Cancer'
confusionMatrix(predict(rf.fit, testing), reference=testing$Type)

pred <- predict(rf.fit, testing)
pred

testing$Type == pred





