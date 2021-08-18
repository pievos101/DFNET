# Create Data Set

# Read in PPI
Protein     <- read.table("9606.protein.links.v11.0 (1).txt", header=TRUE)

# Translate Gene names
Genes       <- read.table("human.name_2_string.csv")

# Read in selected genes from the Nature publication
Known_Genes       <- read.table("multiomics_features.tsv", sep="\t")
cancer_gene_names <- Known_Genes[,1]

ids <- match(cancer_gene_names, Genes[,2])
ids <- ids[!is.na(ids)]

cancer_gene_names_filtered <- Genes[ids,2]
cancer_gene_IDs_filtered   <- Genes[ids,3]

ids1 <- match(Protein[,1], cancer_gene_IDs_filtered)
ids2 <- match(Protein[,2], cancer_gene_IDs_filtered)

Protein_NEW     <- Protein
Protein_NEW[,1] <- cancer_gene_names_filtered[ids1]
Protein_NEW[,2] <- cancer_gene_names_filtered[ids2]

# -- REVISE THIS PART ---
# Protein_NEW <- Protein
#for(xx in 1:length(cancer_gene_IDs_filtered)){
#
# cat(xx, " of ", length(cancer_gene_IDs_filtered), "\n")
# 
# ids1 <- which(Protein[,1]==cancer_gene_IDs_filtered[xx])
# Protein_NEW[ids1,1] <- cancer_gene_names_filtered[xx]
# 
# ids2 <- which(Protein[,2]==cancer_gene_IDs_filtered[xx])
# Protein_NEW[ids2,2] <- cancer_gene_names_filtered[xx]

#}
# --- END REVISE THIS PART ---

# CLEAN DATA

na.ids1 <- which(is.na(Protein_NEW[,1]))
na.ids2 <- which(is.na(Protein_NEW[,2]))

na.IDS  <- unique(c(na.ids1, na.ids2)) 

Protein_NEW_cleaned <- Protein_NEW[-na.IDS,]