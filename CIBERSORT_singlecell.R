# Aim to extract appropriate data from the scRNAseq COMBAT dataset 
#in order to create reference panel for deconvolution using CIBERSORTx

library(tidyverse)


### Step 1: Get the COMBAT rds object ###

rds <- readRDS("COMBAT_sepsis_HV_GEX_metadata.rds")
#there are multiple levels of annotation
#there are 43 samples

dim(rds)
#20615 genes and 243385 cells

#there are no reductions or clustering in this object
Reductions(rds)

### Step 2: select only the HV samples ###

rds_HV <- subset(rds, subset = Source == "HV")
dim(rds_HV)
#20615 genes and 92205 samples

#remove the original object to save space
rm(rds)

#how many different HVs are there?
rds_HV[["gen1ID"]] %>% unique()
#there are 10 HV samples as expected

#filter out any cells with no gene expression
test <- subset(rds_HV, subset = nCount_RNA > 0)
dim(test)
rm(test)
# as expected this does not change the dimensions because filtering has already been done

# how are the cells annotated broadly?
names(rds_HV[[]])
rds_HV[["MM_annot_major_cell_type"]] %>% unique()
#so we have "CD4", "CD8" for T cells
# cMono and ncMono for classical and nonclassical monocytes.
rds_HV[["MM_annot_lineage"]] %>% unique()
# MNP for all mononuclear cells


### Step 3: create separate objects for each cell type ###

#select the CD4 T cells, CD8 T cells, CD14 monocytes and MNP
names(rds_HV[[]])
rds_CD4 <- subset(rds_HV, subset = MM_annot_major_cell_type == "CD4")
rds_CD8 <- subset(rds_HV, subset = MM_annot_major_cell_type == "CD8")
rds_cMono <- subset(rds_HV, subset = MM_annot_major_cell_type == "cMono")
rds_MNP <- subset(rds_HV, subset = MM_annot_lineage == "MNP")
rm(rds_HV)


dim(rds_CD4) #20615 genes and 31225 cells
dim(rds_CD8) #20615 genes and 15745 cells
dim(rds_cMono) #20615 genes and 14806 cells
dim(rds_MNP) # 20615 genes and 20909 cells

#what are the names of the subpopulations
CD4_pseudobulk <- rds_CD4[["MM_annot_pseudobulk"]] %>% unique()
CD4_fine <- rds_CD4[["MM_annot_fine"]] %>% unique()
CD8_pseudobulk <- rds_CD8[["MM_annot_pseudobulk"]] %>% unique()
CD8_fine <- rds_CD8[["MM_annot_fine"]] %>% unique()
cMono_pseudobulk <- rds_cMono[["MM_annot_pseudobulk"]] %>% unique()
cMono_fine <- rds_cMono[["MM_annot_fine"]] %>% unique()
MNP_pseudobulk <- rds_MNP[["MM_annot_pseudobulk"]] %>% unique()
MNP_fine <- rds_MNP[["MM_annot_fine"]] %>% unique()

# there are two levels of annotation, the "pseudobulk" and the "fine" annoations
# refer to Figure 1C in the combat paper

### Step 4: extract the count matrices ###

#create a function
generate_matrix <- function(rds, annotation){
    #save the CD4 count matrix
    counts <- rds@assays$RNA@counts %>% as.matrix()
    print(counts[1:10, 1:10])
    #transpose & convert to df
    counts_t <- t(counts) %>% 
        as.data.frame() %>%
        rownames_to_column(var = "barcode")
    print(counts_t[1:10, 1:10])
    #extract annotations from metadata
    metadata <- rds[[]] %>% 
        select(annotation) %>%
        rownames_to_column(var = "barcode")
    #join the two df
    matrix <- counts_t %>%
        left_join(metadata, by = "barcode") %>%
        select(annotation, 2:length(counts_t))
    print(matrix [1:10, 1:10])
    #transpose back
    matrix_t <- t(matrix) %>% 
        as.data.frame()
    #put the annotation as the colname, and put gene names back as first column
    colnames(matrix_t) <- matrix_t[1,]
    # should add a step here to convert NA to "NA"
    matrix_t <- matrix_t[-1,] %>% 
        rownames_to_column(var = "GeneSymbol")
    print(matrix_t[1:10, 1:10])
    return(matrix_t)
}

# Run the function for each cell type and annotation level
# annotation is either MM_annot_pseudobulk or MM_annot_fine

CD4_matrix_pseudobulk <- generate_matrix(rds_CD4, "MM_annot_pseudobulk")
CD4_matrix_fine <- generate_matrix(rds_CD4, "MM_annot_fine")
CD8_matrix_pseudobulk <- generate_matrix(rds_CD8, "MM_annot_pseudobulk")
CD8_matrix_fine <- generate_matrix(rds_CD8, "MM_annot_fine")
cMono_matrix_pseudobulk <- generate_matrix(rds_cMono, "MM_annot_pseudobulk")
cMono_matrix_fine <- generate_matrix(rds_cMono, "MM_annot_fine")
MNP_matrix_pseudobulk <- generate_matrix(rds_MNP, "MM_annot_pseudobulk")
MNP_matrix_fine <- generate_matrix(rds_MNP, "MM_annot_fine")


# save the output to run in CIBERSORTx
write.table(CD4_matrix_pseudobulk, "CD4_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CD4_matrix_fine, "CD4_matrix_fine.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CD8_matrix_pseudobulk, "CD8_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CD8_matrix_fine, "CD8_matrix_fine.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(cMono_matrix_pseudobulk, "cMono_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(cMono_matrix_fine, "cMono_matrix_fine.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(MNP_matrix_pseudobulk, "MNP_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(MNP_matrix_fine, "MNP_matrix_fine.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# randomly select 10% of the cells for CD4 pseudobulk
CD4_pseudobulk_subset <- CD4_matrix_pseudobulk[,sample(length(CD4_matrix_pseudobulk), length(CD4_matrix_pseudobulk)*0.1)]
# bind back on the GeneSymbol column
CD4_pseudobulk_subset <- cbind (GeneSymbol = CD4_matrix_pseudobulk$GeneSymbol, CD4_pseudobulk_subset)

#save this to test out CIBERSORTx on a smaller matrix
write.table(CD4_pseudobulk_subset, "CD4_pseudobulk_subset.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)


### fixing errors in CD8 pseudobulk

counts <- rds_CD8@assays$RNA@counts %>% as.matrix()
print(counts[1:10, 1:10])
#transpose & convert to df
counts_t <- t(counts) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "barcode")
print(counts_t[1:10, 1:10])
#extract annotations from metadata
metadata <- rds_CD8[[]] %>% 
    select(MM_annot_pseudobulk) %>%
    rownames_to_column(var = "barcode")
#join the two df
matrix <- counts_t %>%
    left_join(metadata, by = "barcode") %>%
    select(MM_annot_pseudobulk, 2:length(counts_t))
print(matrix [1:10, 1:10])
#transpose back
matrix_t <- t(matrix) %>% 
    as.data.frame()
matrix_t[1:10, 1:10]
dim(matrix_t)
#put the annotation as the colname, and put gene names back as first column
colnames(matrix_t) <- matrix_t[1,]
matrix_t2 <- matrix_t[-1,] 
matrix_t2[1:10, 1:10]
test <- matrix_t %>% rownames_to_column(var = "GeneSymbol")
# this throws an error because there is a column called "NA"
which(is.na(colnames(matrix_t2))) # shows the empty column name is at position 600
length(matrix_t2) # there are 15745 columns
#rename column 600 to "NA"
colnames(matrix_t2)[600] <- "NA" 
#convert rownames to column
matrix_t3 <- matrix_t2 %>% rownames_to_column(var = "GeneSymbol")
dim(matrix_t3)
print(matrix_t3[1:10, 1:10])

write.table(matrix_t3, "CD8_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
