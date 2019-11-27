setwd('~/Projects/Infrastructure_20190109/PCa/')

library(GDCRNATools)
library(rtracklayer)
library(tibble)


#===============================================================================================
project <- 'TCGA-PRAD'
rnadir <- paste('data', project, 'RNAseq', sep='/')
mirdir <- paste('data', project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)


####### Download mature miRNA data #######
gdcRNADownload(project.id     = project, 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)



####### Download clinical data #######
clinicaldir <- paste('data', project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = project, 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)



####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = project,
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)



####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

saveRDS(rnaCounts, file='data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')

saveRDS(mirCounts, file='data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')



####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

View(clinicalDa)


saveRDS(clinicalDa, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_11262019.RDS')


sum(metaMatrix.MIR$sample_type=='PrimaryTumor')


pheno <- read.table('data/TCGA-PRAD/phenotype.TCGA-PRAD.final.txt', header=T, stringsAsFactors = F,
                    sep='\t', row.names = 1)

pheno

saveRDS(pheno, file='data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')

sum(!is.na(pheno$days_to_first_biochemical_recurrence))



### Annotation

getGENCODEAnnotation <- function(species='human', release='31', type='gene') {
  
  # species: human, mouse
  # release: human 31, mouse M20
  # type: gene, transcript
  
  baseurl <- 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
  
  gtf.file <- paste0(baseurl, 'Gencode_', species, '/release_', release, '/gencode.v', release, '.annotation.gtf.gz')
  gtf <- readGFF(gtf.file, version=2L)
  
  if (type!='all') {
    gtf <- gtf[gtf$type==type,]
    ensembl <- sapply(gtf$gene_id, function(x) strsplit(x, '.', fixed=T)[[1]][1])
    gtf <- add_column(gtf, ensembl, .before = 'gene_id')
  }
  
  return(gtf)
}


gtf <- getGENCODEAnnotation(species='human', release='32', type='gene')
gtf

rownames(gtf) <- ifelse(duplicated(gtf$ensembl), gtf$gene_id, gtf$ensembl)

gtf <- add_column(.data = gtf, .after = 'end', length=gtf$end-gtf$start+1)

saveRDS(gtf, file='data/GENCODE_Annotation_Human_V32.RDS')










#===================================================================================================

##############################################
rnaCounts <- readRDS('data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')
metaMatrix.RNA <- readRDS('data/TCGA-PRAD/Metadata_RNAseq_TCGA_PRAD.RDS')
clinicalDa <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')
clinicalDa

genecodeAnnotation <- readRDS('data/GENCODE_Annotation_Human_V32.RDS')

dge <-  DGEList(counts = rnaCounts)

### TMM normalization
dge = calcNormFactors(dge, method = 'TMM')

### Filter out low-expression genes (cpm>1 in at least 50% of the samples)
keep <- rowSums(edgeR::cpm(dge) > 1) >= 0.5*ncol(rnaCounts)
sum(keep)
dge <- dge[keep,,keep.lib.sizes = TRUE]

### Voom normalization
v <- voom(dge, design=NULL, plot = FALSE)

exprAfterVoom <- v$E ### for visualization
exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
exprLogCPM

### Prepare comparison matrix
group <- factor(metaMatrix.RNA$sample_type)
group

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(contrasts='PrimaryTumor - SolidTissueNormal',
                                 levels=design)
contrast.matrix

### Differential gene expression analysis (limma)

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


### Report DEGs
dgeTable <- topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by='p')
dgeTable


dgeTable$Symbol <- genecodeAnnotation[rownames(dgeTable),'gene_name']
dgeTable$GeneType <- genecodeAnnotation[rownames(dgeTable),'gene_type']
View(dgeTable)

unique(genecodeAnnotation$gene_type)
table(dgeTable$GeneType[abs(dgeTable$logFC)>=1 & dgeTable$adj.P.Val<0.01])

saveRDS(dgeTable, file='DE_Genes_RNAseq_logFC1_FDR0.01.RDS')


dgeTable[which(duplicated(dgeTable$Symbol)),]


#############

daysToDeath <- as.numeric(clinicalDa$days_to_first_biochemical_recurrence)
daysToDeath

nonComplt <- is.na(daysToDeath)
nonComplt

vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
vitalStatus

daysToDeath[nonComplt] <- as.numeric(clinicalDa$days_to_last_followup[nonComplt])


preopPSA <- clinicalDa$preop_psa

gleason1 <- clinicalDa$primary_pattern
gleason2 <- clinicalDa$secondary_pattern

gleason <- paste(gleason1, gleason2, sep='+')

gleason <- gleason1+gleason2



stage <- clinicalDa$pathologic_stage
stage
stage[stage %in% c('StageI','StageIA','StageIB')] <- 'StageI'
stage[stage %in% c('StageII','StageIIA','StageIIB')] <- 'StageII'
stage[stage %in% c('StageIII','StageIIIA','StageIIIB')] <- 'StageIII'
stage[stage %in% c('StageIV','StageIVA','StageIVB')] <- 'StageIV'
stage <- as.numeric(as.factor(stage))
stage

Tstage <- clinicalDa$pathologic_T
Tstage

age <- as.numeric(as.character(clinicalDa$age))
age


traits <- data.frame(preopPSA, age, )
traits



coxphTRTs <- c()
for (i in 1:ncol(traits)) {
  TRT <- traits[,i]
  coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ TRT)
  #coeffs <- coef(summary(coxtest))[1,]
  #coxphTRTs <- rbind(coxphTRTs, coeffs)
  
  summcph <- summary(coxtest)
  coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,-1], 
              summcph$coefficients[,4:5])
  coxphTRTs <- rbind(coxphTRTs, coeffs)
}


rownames(coxphTRTs) <- colnames(traits)
coxphTRTs

sortedTRTs <- coxphTRTs[order(coxphTRTs[,7]),]
sortedTRTs

write.table(sortedTRTs, file = 'NewResults/univariate.coxph.clinical.TCGA.txt', sep='\t', quote=F)



coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ riskScore+age+gender+smoking+stage)
coxtest
coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ riskScore+stage)
coxtest

summcph <- summary(coxtest)

coeffs <- data.frame(summcph$coefficients[,1:2], summcph$conf.int[,-1], 
                     summcph$coefficients[,4:5])
coeffs

write.table(coeffs, file = 'NewResults/multivariate.coxph.clinical.TCGA.txt', sep='\t', quote=F)

coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ riskScore3)
summcph <- summary(coxtest)

c(summcph$coefficients[,1:2], summcph$conf.int[,-1], 
  summcph$coefficients[,4:5])






coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ factor(Tstage))
summary(coxtest)



library(survival)
library(survminer)
                   
#risk+tumorRecurrence+
#                   tnm.t+radio+tnm.m+residualTumor+tobHist)


riskGroup <- Tstage

riskGroup <- ifelse(preopPSA < 10, 0, 1)


sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ riskGroup)
pval1 <- format(pchisq(sdf$chisq, length(sdf$n)-1, lower.tail = FALSE),digits=5)
pval1

HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))

HR
upper95
lower95


survData <- data.frame(daysToDeath, vitalStatus, riskGroup)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
ggplot2::annotate("text", 
                  x = 2200, y = 1, # x and y coordinates of the text
                  label = "p = 4.196e-06", size = 4)

p




#############################################

######## CoxPH for DEGs


ovlp <- intersect(colnames(exprLogCPM), rownames(clinicalDa))
ovlp

expr <- exprLogCPM[,ovlp]
clinical <- clinicalDa[ovlp,]

DEGs <- rownames(dgeTable)[abs(dgeTable$logFC)>=1 & dgeTable$adj.P.Val<0.01]

daysToDeath <- as.numeric(clinical$days_to_first_biochemical_recurrence)
nonComplt <- is.na(daysToDeath)
vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
daysToDeath[nonComplt] <- as.numeric(clinical$days_to_last_followup[nonComplt])

coeffs <- c()
for (gene in DEGs) {
  coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ unlist(expr[gene,]))
  summcph <- summary(coxtest)
  
  coeffs <- rbind(coeffs, c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
                            summcph$coefficients[,4:5]))
  
}

rownames(coeffs) <- DEGs
colnames(coeffs) <- c('coef','HR','lower95','upper95', 'z', 'p')
coeffs <- data.frame(coeffs, stringsAsFactors = F)
coeffs$FDR <- p.adjust(coeffs$p, method='BH')

coeffs$Symbol <- dgeTable[DEGs,'Symbol']
coeffs$GeneType <- dgeTable[DEGs,'GeneType']

saveRDS(coeffs, file='data/Univariate_CoxPH_RNAseq_DEGs.RDS')

View(coeffs[coeffs$FDR<0.05,])


g <- 'ENSG00000225937'

riskGroup <- ifelse(expr[g, ] < median(expr[g,]), 0, 1)

survData <- data.frame(daysToDeath, vitalStatus, riskGroup)

fit <- survfit(Surv(daysToDeath, vitalStatus) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p




coeffsLNC <- coeffs[which(coeffs$GeneType=='lncRNA' & coeffs$FDR<0.01),]
coeffsLNC

exprLNC <- expr[rownames(coeffsLNC),]


coeffsALL <- coeffs[which(coeffs$FDR<0.01),]
coeffsALL

exprALL <- expr[rownames(coeffsALL),]

###############################################

library(glmnet)

?cv.glmnet()


####### only for lasso
filter <- which(is.na(daysToDeath) | daysToDeath==0)
filter

samples <- ovlp[-filter]
samples


daysToDeathLasso <- daysToDeath[-filter]
vitalStatusLasso <- vitalStatus[-filter]
exprLasso <- exprLNC[,-filter]
exprLasso <- exprALL[,-filter]
exprLasso

x <- as.matrix(t(exprLasso))
x[1:7,1:7]

y <- as.matrix(data.frame(time=daysToDeathLasso,status=vitalStatusLasso))
y


nfold<-10

foldidID<-lapply(1:11,function(i){
  n<-nrow(y)
  sample(rep(1:nfold,ceiling(n/nfold))[1:n])
})



cv.fit<-cv.glmnet(x=x, y=y, family='cox', foldid=foldidID[[3]], alpha=1) # Lasso: alpha=1; Ridge: alpha=0
plot(cv.fit)


coef.min<-coef(cv.fit,s="lambda.min")
coef.min
cvfit$lambda.min

active.min <- which(as.numeric(coef.min) !=0)
active.min

lassoGene <- coef.min@Dimnames[[1]][active.min]
lassoGene

output <- data.frame(gene=coef.min@Dimnames[[1]],coeffs=matrix(coef.min), stringsAsFactors = F, 
                     symbol=genecodeAnnotation[coef.min@Dimnames[[1]], 'gene_name'],gene_type=genecodeAnnotation[coef.min@Dimnames[[1]], 'gene_type'],
                     row.names=coef.min@Dimnames[[1]])

output <- output[lassoGene,]
output

exprLassoModel <- exprLasso[lassoGene,]

riskScore = apply(exprLassoModel, 2, function(x) sum(x*output$coeffs))
riskScore



coxtest <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskScore)
summary(coxtest)


riskGroup <- ifelse(riskScore < median(riskScore), 0, 1)

survData <- data.frame(daysToDeathLasso, vitalStatusLasso, riskGroup)

fit <- survfit(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p



###



multiCox <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ 
                    exprLassoModel[1,]+exprLassoModel[2,]+exprLassoModel[3,]+exprLassoModel[4,]+
                    exprLassoModel[5,]+exprLassoModel[6,]+exprLassoModel[7,]+exprLassoModel[8,]+
                    exprLassoModel[9,]+exprLassoModel[10,]+exprLassoModel[11,]+exprLassoModel[12,]+
                    exprLassoModel[13,]+exprLassoModel[14,]+exprLassoModel[15,]+exprLassoModel[16,])


summcph <- summary(multiCox)

exprLassoModel



multiCox <- data.frame(cbind(summcph$coefficients[,1:2], summcph$conf.int[,3:4], summcph$coefficients[,4:5]),
                       stringsAsFactors = F)
rownames(multiCox) <- rownames(exprLassoModel)
multiCox

colnames(multiCox) <- c('coef','HR','lower95','upper95', 'z', 'p')



riskScore = apply(exprLassoModel, 2, function(x) sum(x*multiCox$coef))
riskScore


coxtest <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskScore)
summary(coxtest)


riskGroup <- ifelse(riskScore < median(riskScore), 0, 1)

survData <- data.frame(daysToDeathLasso, vitalStatusLasso, riskGroup)

fit <- survfit(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p


#======================================================================================================

###### GSVA

signature <- list('All'=lassoGene)
signature

gsvaScore2 <- gsva(expr = as.matrix(exprALL[,-filter]),
                  gset.idx.list = signature,
                  method="gsva")

length(gsvaScore2)
length(daysToDeathLasso)

gsvaScore <- as.numeric(gsvaScore)

coxtest <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ as.numeric(gsvaScore))
summary(coxtest)


riskGroup <- ifelse(as.numeric(gsvaScore) < median(gsvaScore), 0, 1)

survData <- data.frame(daysToDeathLasso, vitalStatusLasso, riskGroup)

fit <- survfit(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p

cor.test(gsvaScore, gsvaScore2)

plot(gsvaScore, gsvaScore2)

#======================================================================================================

cr <- read_xlsx('data/Classifiers.xlsx', sheet='Oncotype')
cr

cr <- data.frame(cr, row.names=1, stringsAsFactors = F)
cr

rownames(cr) %in% genecodeAnnotation$gene_name

idx <- match(rownames(cr), genecodeAnnotation$gene_name)
rownames(genecodeAnnotation)[idx]

cr$symbol <- rownames(cr)
rownames(cr) <- rownames(genecodeAnnotation)[idx]


exprCRModel <- exprLogCPM[rownames(cr),samples]
exprCRModel


riskScore = apply(exprCRModel, 2, function(x) sum(x*cr$coeffs))
riskScore


coxtest <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskScore)
summary(coxtest)


riskGroup <- ifelse(riskScore < median(riskScore), 0, 1)

survData <- data.frame(daysToDeathLasso, vitalStatusLasso, riskGroup)

fit <- survfit(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p


###

multiCox <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ 
                    exprCRModel[1,]+exprCRModel[2,]+exprCRModel[3,]+exprCRModel[4,]+
                    exprCRModel[5,]+exprCRModel[6,]+exprCRModel[7,]+exprCRModel[8,]+
                    exprCRModel[9,]+exprCRModel[10,]+exprCRModel[11,]+exprCRModel[12,]+
                    exprCRModel[13,]+exprCRModel[14,]+exprCRModel[15,]+exprCRModel[16,]+
                    exprCRModel[17,]+exprCRModel[18,]+exprCRModel[19,]+exprCRModel[20,]+
                    exprCRModel[21,]+exprCRModel[22,]+exprCRModel[23,]+exprCRModel[24,])


summcph <- summary(multiCox)

multiCox <- data.frame(cbind(summcph$coefficients[,1:2], summcph$conf.int[,3:4], summcph$coefficients[,4:5]),
                       stringsAsFactors = F)
rownames(multiCox) <- rownames(exprCRModel)
multiCox

colnames(multiCox) <- c('coef','HR','lower95','upper95', 'z', 'p')



riskScore = apply(exprCRModel, 2, function(x) sum(x*multiCox$coef))
riskScore


coxtest <- coxph(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskScore)
summary(coxtest)


riskGroup <- ifelse(riskScore < median(riskScore), 0, 1)

survData <- data.frame(daysToDeathLasso, vitalStatusLasso, riskGroup)

fit <- survfit(Surv(daysToDeathLasso, vitalStatusLasso) ~ riskGroup, data=survData)
fit

#########
p <- ggsurvplot(fit, pval = TRUE, pval.coord = c(2200, 1), 
                font.main = c(14, 'bold', 'blue'), conf.int = FALSE, legend = c(0.15, 0.2), 
                #legend.labs = c('Low risk', 'High risk'),  legend.title='',
                xlab = 'OS (days)', ylab = 'Survival probability',
                font.x = c(12), font.y = c(12),xlim=c(0,5000)) #+
p







