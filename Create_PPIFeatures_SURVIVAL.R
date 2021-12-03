#For linkedOmics data
#Create Feature Matrices

Methy <- read.table("Human__TCGA_OV__JHU_USC__Methylation__Meth450__01_28_2016__BI__Gene__Firehose_Methylation_Prepocessor.cct", sep="\t", header=TRUE)
mRNA  <- read.table("Human__TCGA_OV__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", sep="\t", header=TRUE)
#Mut   <- read.table("Human__TCGA_OV__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt", sep="\t", header=TRUE)

# Methylation
Methy_genes <- Methy[,1]
Methy <- Methy[,-1]
Methy <- t(Methy)
colnames(Methy) <- Methy_genes
Methy_patients <- rownames(Methy)

# mRNA
mRNA_genes <- mRNA[,1]
mRNA <- mRNA[,-1]
mRNA <- t(mRNA)
colnames(mRNA) <- mRNA_genes
mRNA_patients  <- rownames(mRNA)

# Mutation (SNV)
#Mut_genes <- Mut[,1]
#Mut <- Mut[,-1]
#Mut <- t(Mut)
#colnames(Mut) <- Mut_genes
#Mut_patients  <- rownames(Mut)


#Harmonize Multi-Omics data
common_genes    <- intersect(Methy_genes, mRNA_genes)
common_patients <- intersect(Methy_patients, mRNA_patients)
#common_genes    <- intersect(intersect(Methy_genes, mRNA_genes), Mut_genes)
#common_patients <- intersect(intersect(Methy_patients, mRNA_patients), Mut_patients)

mRNA2  <- mRNA[common_patients,common_genes]
Methy2 <- Methy[common_patients,common_genes]
#Mut2   <- Mut[common_patients,common_genes]

# READ IN THE PPI
PPI <- read.table("~/RF-on-Graphs-Project/Application/PPI_processed2.txt")

# HARMONIZE WITH PPI (1)
ids1    <- match(PPI[,1], common_genes)
ids2    <- match(PPI[,2], common_genes)

na.ids1 <- which(is.na(ids1))
na.ids2 <- which(is.na(ids2))

na.IDS  <- unique(c(na.ids1,na.ids2))

PPI2    <- PPI[-na.IDS,]


# HARMONIZE WITH THE PPI (2)
PPI_genes   <- unique(c(PPI2[,1],PPI2[,2]))
ids         <- match(common_genes, PPI_genes)
na.ids      <- which(is.na(ids))
mRNA3       <- mRNA2[,-na.ids]
Methy3      <- Methy2[,-na.ids]
#Mut3        <- Mut2[, -na.ids]

# READ IN MEDICAL SURVIVAL DATA
CLIN <- read.table("Human__TCGA_OV__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi.txt", sep="\t", header=TRUE)

status_id <- which(CLIN[,1]=="status")
survival <- CLIN[status_id,-1]

ids      <- match(rownames(mRNA3), names(survival))
survival_final <- survival[ids]
na.ids <- which(is.na(survival_final))

del.patients <- names(survival_final)[na.ids]
del.ids      <- match(del.patients, rownames(mRNA3))

# These are the final objects
mRNA4        <- mRNA3[-del.ids,]
Methy4       <- Methy3[-del.ids,]
#Mut4         <- Mut3[-del.ids,]
SURVIVAL     <- survival_final[-na.ids]
PPI_FINAL    <- PPI2


write.table(mRNA4,     file="OV_mRNA_FEATURES.txt")
write.table(Methy4,    file="OV_Methy_FEATURES.txt")
#write.table(Mut4,      file="BRCA_Mut_FEATURES.txt")
write.table(SURVIVAL,  file="OV_SURVIVAL.txt")
write.table(PPI_FINAL, file="OV_PPI.txt")

