# random TCGA
# Harmonize two cancer types

KIDNEY_PPI <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")
KIDNEY_mRNA <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")
KIDNEY_Methy <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")

# OV_PPI      <- read.table("~/LinkedOmics/OV/OV_PPI.txt")
# OV_mRNA     <- read.table("~/LinkedOmics/OV/OV_mRNA_FEATURES.txt")
# OV_Methy    <- read.table("~/LinkedOmics/OV/OV_Methy_FEATURES.txt")

LUAD_PPI <- read.table("~/LinkedOmics/LUAD/LUAD_PPI.txt")
LUAD_mRNA <- read.table("~/LinkedOmics/LUAD/LUAD_mRNA_FEATURES.txt")
LUAD_Methy <- read.table("~/LinkedOmics/LUAD/LUAD_Methy_FEATURES.txt")

BRCA_PPI <- read.table("~/LinkedOmics/BRCA/BRCA_PPI.txt")
BRCA_mRNA <- read.table("~/LinkedOmics/BRCA/BRCA_mRNA_FEATURES.txt")
BRCA_Methy <- read.table("~/LinkedOmics/BRCA/BRCA_Methy_FEATURES.txt")

# TARGET   <- read.table("~/LinkedOmics/KIRC/KIDNEY_SURVIVAL.txt"

KIDNEY_GENENAMES <- colnames(KIDNEY_mRNA)
# OV_GENENAMES     <- colnames(OV_mRNA)
LUAD_GENENAMES <- colnames(LUAD_mRNA)
BRCA_GENENAMES <- colnames(BRCA_mRNA)


# COMMON_NAMES  <- intersect(intersect(intersect(KIDNEY_GENENAMES, OV_GENENAMES),
# 					LUAD_GENENAMES),BRCA_GENENAMES)

COMMON_NAMES <- intersect(intersect(
    KIDNEY_GENENAMES,
    LUAD_GENENAMES
), BRCA_GENENAMES)

# Filter data accordingly

KIDNEY_mRNA2 <- KIDNEY_mRNA[, COMMON_NAMES]
KIDNEY_Methy2 <- KIDNEY_Methy[, COMMON_NAMES]

# OV_mRNA2  <- OV_mRNA[,COMMON_NAMES]
# OV_Methy2 <- OV_Methy[,COMMON_NAMES]

LUAD_mRNA2 <- LUAD_mRNA[, COMMON_NAMES]
LUAD_Methy2 <- LUAD_Methy[, COMMON_NAMES]

BRCA_mRNA2 <- BRCA_mRNA[, COMMON_NAMES]
BRCA_Methy2 <- BRCA_Methy[, COMMON_NAMES]

# HARMONIZE WITH PPI (1) # here just one cancer PPI is needed
ids1 <- match(KIDNEY_PPI[, 1], COMMON_NAMES)
ids2 <- match(KIDNEY_PPI[, 2], COMMON_NAMES)

na.ids1 <- which(is.na(ids1))
na.ids2 <- which(is.na(ids2))

na.IDS <- unique(c(na.ids1, na.ids2))

PPI2 <- KIDNEY_PPI[-na.IDS, ]

# HARMONIZE WITH THE PPI (2)
PPI_genes <- unique(c(PPI2[, 1], PPI2[, 2]))
ids <- match(COMMON_NAMES, PPI_genes)
na.ids <- which(is.na(ids))

KIDNEY_mRNA3 <- KIDNEY_mRNA2[, -na.ids]
KIDNEY_Methy3 <- KIDNEY_Methy2[, -na.ids]

# OV_mRNA3       <- OV_mRNA2[,-na.ids]
# OV_Methy3      <- OV_Methy2[,-na.ids]

LUAD_mRNA3 <- LUAD_mRNA2[, -na.ids]
LUAD_Methy3 <- LUAD_Methy2[, -na.ids]

BRCA_mRNA3 <- BRCA_mRNA2[, -na.ids]
BRCA_Methy3 <- BRCA_Methy2[, -na.ids]


ALL_mRNA <- rbind(KIDNEY_mRNA3, LUAD_mRNA3, BRCA_mRNA3)
ALL_Methy <- rbind(KIDNEY_Methy3, LUAD_Methy3, BRCA_Methy3)

# TARGET <- numeric(dim(KIDNEY_mRNA3)[1]+dim(OV_mRNA3)[1]+dim(KIDNEY_mRNA3)[1]+dim(OV_mRNA3)[1])
# TARGET[1:dim(KIDNEY_mRNA3)[1]] <- 1
# names(TARGET) <- rownames(ALL_mRNA)

# randomly select 200 patients
ids <- sample(dim(ALL_mRNA)[1], 200)

SAMP_mRNA <- ALL_mRNA[ids, ]
SAMP_Methy <- ALL_Methy[ids, ]


write.table(SAMP_mRNA, file = "RANDOM_mRNA_FEATURES.txt")
write.table(SAMP_Methy, file = "RANDOM_Methy_FEATURES.txt")
# write.table(t(TARGET),  file="KIDNEY_OV_TARGET.txt")
write.table(PPI2, file = "RANDOM_PPI.txt")
