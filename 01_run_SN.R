library(data.table)

source("/mnt/lareaulab/reliscu/code/SampleNetwork/SampleNetwork_1.08.r")

setwd("/mnt/lareaulab/reliscu/projects/NSF_GRFP/analyses/bulk/GTEx/frontal_cortex")

datExprT <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/GTEx/frontal_cortex/GTEx_frontal_cortex_gene_TPM.csv", data.table=FALSE)
sampleinfo1 <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/GTEx/frontal_cortex/GTEx_frontal_cortex_sampleinfo.csv", data.table=FALSE)

sampleinfo1$Mean_age <- sapply(strsplit(sampleinfo1$AGE, "-"), function(x) mean(as.numeric(x)))
sampleinfo1$grouplabels1 <- "All"

# Order columns by variable they will be colored by:

sampleinfo1 <- sampleinfo1[order(sampleinfo1$Mean_age),]
datExprT <- datExprT[, c(1, match(sampleinfo1[,1], colnames(datExprT)))]

projectname1 <- "GTEx_frontal_cortex_TPM"
skip1 <- 1 # An integer describing the number of feature information columns
indices1 <- list(seq(2, ncol(datExprT)))
samplelabels1 <- 1 # An integer that points to the column number in sampleinfo1 containing the sample labels that will appear in plots.  Note: these 					sample labels must be identical to the sample column headers in datExprT or an error will be triggered.
grouplabels1 <- grep("grouplabels1", colnames(sampleinfo1))
subgroup1 <- grep("AGE", colnames(sampleinfo1)) # grep("SMCENTER", colnames(sampleinfo1))

btrait1 <- c(
  grep("SMCENTER", colnames(sampleinfo1)),
  grep("SMNABTCH", colnames(sampleinfo1)),
  grep("SEX", colnames(sampleinfo1)),
  grep("DTHHRDY", colnames(sampleinfo1)),
  grep("SMRIN", colnames(sampleinfo1)),
  grep("Mean_age", colnames(sampleinfo1))
)

asfactors1 <- c(
  grep("SMCENTER", colnames(sampleinfo1)),
  grep("SMNABTCH", colnames(sampleinfo1)),
  grep("SEX", colnames(sampleinfo1)),
  grep("DTHHRDY", colnames(sampleinfo1))
)

SampleNetwork(
  datExprT=datExprT,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=skip1,
  indices1=indices1,
  sampleinfo1=sampleinfo1,
  subgroup1=subgroup1,
  samplelabels1=samplelabels1,
  grouplabels1=grouplabels1,
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=btrait1,
  trait1=NULL,
  asfactors1=NULL,
  projectname1=projectname1,
  cexlabels1=0.7,
  normalize1=TRUE,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
