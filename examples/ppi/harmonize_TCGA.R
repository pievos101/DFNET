# Harmonize two cancer types

KIDNEY_PPI <- read.table("~/LinkedOmics/KIRC/KIDNEY_PPI.txt")
KIDNEY_mRNA <- read.table("~/LinkedOmics/KIRC/KIDNEY_mRNA_FEATURES.txt")
KIDNEY_Methy <- read.table("~/LinkedOmics/KIRC/KIDNEY_Methy_FEATURES.txt")

OV_PPI <- read.table("~/LinkedOmics/RANDOM/RANDOM_PPI.txt")
OV_mRNA <- read.table("~/LinkedOmics/RANDOM/RANDOM_mRNA_FEATURES.txt")
OV_Methy <- read.table("~/LinkedOmics/RANDOM/RANDOM_Methy_FEATURES.txt")

# TARGET   <- read.table("~/LinkedOmics/KIRC/KIDNEY_SURVIVAL.txt"

KIDNEY_GENENAMES <- colnames(KIDNEY_mRNA)
OV_GENENAMES <- colnames(OV_mRNA)

COMMON_NAMES <- intersect(KIDNEY_GENENAMES, OV_GENENAMES)

# Filter data accordingly

KIDNEY_mRNA2 <- KIDNEY_mRNA[, COMMON_NAMES]
KIDNEY_Methy2 <- KIDNEY_Methy[, COMMON_NAMES]

OV_mRNA2 <- OV_mRNA[, COMMON_NAMES]
OV_Methy2 <- OV_Methy[, COMMON_NAMES]

# HARMONIZE WITH PPI (1)
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

if (length(na.ids) > 0) {
    KIDNEY_mRNA3 <- KIDNEY_mRNA2[, -na.ids]
    KIDNEY_Methy3 <- KIDNEY_Methy2[, -na.ids]
    OV_mRNA3 <- OV_mRNA2[, -na.ids]
    OV_Methy3 <- OV_Methy2[, -na.ids]
} else {
    KIDNEY_mRNA3 <- KIDNEY_mRNA2
    KIDNEY_Methy3 <- KIDNEY_Methy2
    OV_mRNA3 <- OV_mRNA2
    OV_Methy3 <- OV_Methy2
}

ALL_mRNA <- rbind(KIDNEY_mRNA3, OV_mRNA3)
ALL_Methy <- rbind(KIDNEY_Methy3, OV_Methy3)

TARGET <- numeric(dim(KIDNEY_mRNA3)[1] + dim(OV_mRNA3)[1])
TARGET[1:dim(KIDNEY_mRNA3)[1]] <- 1
names(TARGET) <- rownames(ALL_mRNA)

write.table(ALL_mRNA, file = "KIDNEY_RANDOM_mRNA_FEATURES.txt")
write.table(ALL_Methy, file = "KIDNEY_RANDOM_Methy_FEATURES.txt")
write.table(t(TARGET), file = "KIDNEY_RANDOM_TARGET.txt")
write.table(PPI2, file = "KIDNEY_RANDOM_PPI.txt")
