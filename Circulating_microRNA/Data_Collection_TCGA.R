setwd('~/bigdata/EarlyCancerDetection/CCMA/')


library(GDCRNATools)
library(edgeR)

#clinicalData <- read.delim('data/fromTCGA/clinical_PANCAN_patient_with_followup.tsv', header = T, stringsAsFactors = F, sep = '\t')
#View(clinicalData)

#clinicalData[clinicalData=='[Not Applicable]'] <- NA
#clinicalData[clinicalData=='[Not Available]'] <- NA
#clinicalData[clinicalData=='[Unknown]'] <- NA

#clinicalData[clinicalData=='[Not Evaluated]'] <- NA
#clinicalData[clinicalData=='[Not Reported]'] <- NA
#clinicalData[clinicalData=='[Completed]'] <- NA
#clinicalData[clinicalData=='[Discrepancy]'] <- NA
#dim(clinicalData)

cdrData <- read_excel('data/fromTCGA/TCGA-CDR-SupplementalTableS1.xlsx', sheet = 'TCGA-CDR')
cdrData <- data.frame(cdrData, stringsAsFactors = F)

cdrData[cdrData=='[Not Applicable]'] <- NA
cdrData[cdrData=='[Not Available]'] <- NA
cdrData[cdrData=='[Unknown]'] <- NA

cdrData[cdrData=='[Not Evaluated]'] <- NA
cdrData[cdrData=='[Not Reported]'] <- NA
cdrData[cdrData=='[Completed]'] <- NA
cdrData[cdrData=='[Discrepancy]'] <- NA

cdrData <- cdrData[,-1]

View(cdrData)


projects <- read.table('data/fromTCGA/projects-table.2020-06-20.tsv', header = T, sep = '\t', stringsAsFactors = F)
projects$Project

####### Download clinical data #######
clinicaldir <- paste('data/fromTCGA/Clinical')

for (project in projects$Project) {
  
  gdcClinicalDownload(project.id     = project, 
                      write.manifest = FALSE,
                      method         = 'gdc-client',
                      directory      = clinicaldir)
}

####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
clinicalDa[1:6,5:10]

View(clinicalDa)


saveRDS(clinicalDa, file='data/fromTCGA/Clinical_All_06202020.RDS')

View(clinicalDa)


mir21 <- readRDS('data/miRBase/mir21.RDS')
idx <- match(rownames(mirCounts), mir21$Name)
idx

mir.id <- mir21$ID[idx]



#############################

for (project in projects$Project) {
  if (project=='TCGA-PRAD') {
    next
    
  }
  
  
  mirdir <- paste('data/fromTCGA', project, sep='/')
  mirdir
  
  ####### Download mature miRNA data #######
  gdcRNADownload(project.id     = project, 
                 data.type      = 'miRNAs', 
                 write.manifest = FALSE,
                 method         = 'gdc-client',
                 directory      = mirdir)
  
  
  ####### Parse miRNAs metadata #######
  metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                     data.type  = 'miRNAs', 
                                     write.meta = FALSE)
  
  ####### Filter duplicated samples in miRNAs metadata #######
  metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
  
  if (project!='TCGA-LAML') {
    ####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
    metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)
    metaMatrix.MIR
  }
  

  
  idx1 <- match(metaMatrix.MIR$patient, rownames(clinicalDa))
  idx1
  
  idx2 <- match(metaMatrix.MIR$patient, cdrData$bcr_patient_barcode)
  idx2
  
  phenoData <- data.frame(metaMatrix.MIR[,c(3:5,7)], clinicalDa[idx1,c(1:8,17:19,24:28)], cdrData[idx2, c(25:32)], project_id=project)
  
  saveRDS(phenoData, file=paste0('data/fromTCGA/rData/Clinical_', project, '.RDS'))
  
  
  ####### Merge miRNAs data #######
  mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                           path      = mirdir, # the folder in which the data stored
                           organized = FALSE, # if the data are in separate folders
                           data.type = 'miRNAs')
  
  rownames(mirCounts) <- mir.id
  mirCounts[1:5,1:5]
  
  
  saveRDS(mirCounts, file=paste0('data/fromTCGA/rData/miRNA_Counts_', project, '.RDS'))
  
  
  ##################
  
  #rnaCounts <- readRDS('data/TCGA-PRAD/RNAseq_Counts_TCGA_PRAD.RDS')
  #metaMatrix.RNA <- readRDS('data/TCGA-PRAD/Metadata_RNAseq_TCGA_PRAD.RDS')
  #mirCounts <- readRDS('data/TCGA-PRAD/miRNA_Counts_TCGA_PRAD.RDS')
  #metaMatrix.MIR <- readRDS('data/TCGA-PRAD/Metadata_miRNAs_TCGA_PRAD.RDS')
  #clinicalDa <- readRDS('data/TCGA-PRAD/Clinical_TCGA_PRAD_With_PreopPSA_and_BCR.RDS')
  
  dge <-  DGEList(counts = mirCounts)
  
  ### TMM normalization
  dge = calcNormFactors(dge, method = 'TMM')
  
  exprLogCPM <- edgeR::cpm(dge,log = TRUE) ### for visualization
  dim(exprLogCPM)
  exprLogCPM[1:5,1:5]
  
  saveRDS(exprLogCPM, file=paste0('data/fromTCGA/rData/miRNA_LogCPM_', project, '.RDS'))
  
  
}


