# Aim to extract appropriate data from the scRNAseq COMBAT dataset 
#in order to create reference panel for deconvolution using CIBERSORTx

library(tidyverse)
setwd("/project/tendonhca/ccohen/CIBERSORTx/ankspond")


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
rds_CD14 <- subset(rds_HV, subset = MM_annot_major_cell_type == "cMono" | MM_annot_major_cell_type == "ncMono")
rm(rds_HV)


dim(rds_CD4) #20615 genes and 31225 cells
dim(rds_CD8) #20615 genes and 15745 cells
dim(rds_cMono) #20615 genes and 14806 cells
dim(rds_MNP) # 20615 genes and 20909 cells
dim(rds_CD14) # 20615 genes and 18879 cells

#what are the names of the subpopulations
CD4_pseudobulk <- rds_CD4[["MM_annot_pseudobulk"]] %>% unique()
CD4_fine <- rds_CD4[["MM_annot_fine"]] %>% unique()
CD8_pseudobulk <- rds_CD8[["MM_annot_pseudobulk"]] %>% unique()
CD8_fine <- rds_CD8[["MM_annot_fine"]] %>% unique()
cMono_pseudobulk <- rds_cMono[["MM_annot_pseudobulk"]] %>% unique()
cMono_fine <- rds_cMono[["MM_annot_fine"]] %>% unique()
MNP_pseudobulk <- rds_MNP[["MM_annot_pseudobulk"]] %>% unique()
MNP_fine <- rds_MNP[["MM_annot_fine"]] %>% unique()
CD14_pseudobulk <- rds_CD14[["MM_annot_pseudobulk"]] %>% unique()
CD14_fine <- rds_CD14[["MM_annot_fine"]] %>% unique()

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
CD14_matrix_pseudobulk <- generate_matrix(rds_CD14, "MM_annot_pseudobulk")
CD14_matrix_fine <- generate_matrix(rds_CD14, "MM_annot_fine")

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
write.table(CD14_matrix_pseudobulk, "CD14_matrix_pseudobulk.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(CD14_matrix_fine, "CD14_matrix_fine.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# randomly select 10% of the cells for CD4 pseudobulk
CD4_pseudobulk_subset <- CD4_matrix_pseudobulk[,sample(length(CD4_matrix_pseudobulk), length(CD4_matrix_pseudobulk)*0.1)]
# bind back on the GeneSymbol column
CD4_pseudobulk_subset <- cbind (GeneSymbol = CD4_matrix_pseudobulk$GeneSymbol, CD4_pseudobulk_subset)

#save this to test out CIBERSORTx on a smaller matrix
write.table(CD4_pseudobulk_subset, "CD4_pseudobulk_subset.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)


### Step 5: Read in these outputs to continue formatting the data ###

CD4_matrix_pseudobulk <- read.table("CD4_matrix_pseudobulk.txt", row.names = 1, sep = "\t", header = TRUE)
CD4_matrix_fine <- read.table("CD4_matrix_fine.txt", row.names = 1, sep = "\t", header = TRUE)
CD8_matrix_pseudobulk <- read.table("CD8_matrix_pseudobulk.txt", row.names = 1, sep = "\t", header = TRUE)
CD8_matrix_fine <- read.table("CD8_matrix_fine.txt", row.names = 1, sep = "\t", header = TRUE)
cMono_matrix_pseudobulk <- read.table("cMono_matrix_pseudobulk.txt", row.names = 1, sep = "\t", header = TRUE)
cMono_matrix_fine <- read.table("cMono_matrix_fine.txt", row.names = 1, sep = "\t", header = TRUE)
MNP_matrix_pseudobulk <- read.table("MNP_matrix_pseudobulk.txt", row.names = 1, sep = "\t", header = TRUE)
MNP_matrix_fine <- read.table("MNP_matrix_fine.txt", row.names = 1, sep = "\t", header = TRUE)
CD14_matrix_pseudobulk <- read.table("CD14_matrix_pseudobulk.txt", row.names = 1, sep = "\t", header = TRUE)
CD14_matrix_fine <- read.table("CD14_matrix_fine.txt", row.names = 1, sep = "\t", header = TRUE)

CD8_matrix_pseudobulk_test <- CD8_matrix_pseudobulk[1:100, 1:800]

#create a function to clean up the matrix
cleanup_matrix <- function (matrix){
    # convert "." to "_" in colnames
    colnames(matrix) <- colnames(matrix) %>% str_replace_all("\\.(\\D)", "_\\1")
    print(head(colnames(matrix)))
    
    # remove any columns with "NA"
    print("dimensions before NA removal")
    print(dim(matrix))
    matrix_2 <- matrix %>% select(!starts_with("NA"))
    print("dimensions after NA removal")
    print(dim(matrix_2))
    
    # remove any rows that sum to zero
    matrix_2 <- matrix_2 %>% slice (-which(rowSums(matrix_2) == 0))
    print("dimensions following removal of rows that sum to zero")
    print(dim(matrix_2))
    
    # convert rownames to col 1
    matrix_2 <- matrix_2 %>% rownames_to_column(var = "GeneSymbol")
    print(matrix_2[1:10, 1:10])
    
    return(matrix_2)
}

test <- cleanup_matrix(CD8_matrix_pseudobulk_test)
CD8_matrix_pseudobulk_clean <- cleanup_matrix(CD8_matrix_pseudobulk)
CD4_matrix_pseudobulk_clean <- cleanup_matrix(CD4_matrix_pseudobulk)
cMono_matrix_pseudobulk_clean <- cleanup_matrix(cMono_matrix_pseudobulk)
MNP_matrix_pseudobulk_clean <- cleanup_matrix(MNP_matrix_pseudobulk)
CD14_matrix_pseudobulk_clean <- cleanup_matrix(CD14_matrix_pseudobulk)
CD8_matrix_fine_clean <- cleanup_matrix(CD8_matrix_fine)
CD4_matrix_fine_clean <- cleanup_matrix(CD4_matrix_fine)
cMono_matrix_fine_clean <- cleanup_matrix(cMono_matrix_fine)
MNP_matrix_fine_clean <- cleanup_matrix(MNP_matrix_fine)

dir.create("clean_matrices", showWarnings = FALSE)
write.table (CD8_matrix_pseudobulk_clean, "./clean_matrices/CD8_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (CD4_matrix_pseudobulk_clean, "./clean_matrices/CD4_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (cMono_matrix_pseudobulk_clean, "./clean_matrices/cMono_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (MNP_matrix_pseudobulk_clean, "./clean_matrices/MNP_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (CD14_matrix_pseudobulk_clean, "./clean_matrices/CD14_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (cMono_matrix_fine_clean, "./clean_matrices/MNP_matrix_pseudobulk.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#function to check how many columns of each cell type are present

celltypes_count <- function(df){
    celltypes <- colnames (df) %>% 
        str_split_fixed("\\.", n=2) %>%
        as.data.frame() %>%
        select(V1)%>%
        unique() 
    celltypes <- celltypes[-1,]
    
    print(celltypes)
    print(paste0("there are ", length(celltypes), " cell types"))
    
    column_names <- colnames (df) %>% 
        str_split_fixed("\\.", n=2) %>%
        as.data.frame() %>%
        select(V1) 
    column_names <- column_names[-1,]
    
    #print(head(column_names))

    for (i in seq(1:length(celltypes))){
        celltype <- celltypes[i]
        print(celltype)
        regex <- paste0("^", celltype, "$")
        #print(regex)
        length <- grep (regex, column_names) %>% length()
        print(length)
        i <- i+1
    }    
}
    
celltypes_CD4_pseudobulk <- celltypes_count(CD4_matrix_pseudobulk_clean)
#smallest subset is CD4_TEFF_prolif 228
celltypes_CD8_pseudobulk <- celltypes_count(CD8_matrix_pseudobulk_clean)
celltypes_cMono_pseudobulk <- celltypes_count(cMono_matrix_pseudobulk_clean)
#cMono_cyc 117
celltypes_MNP_pseudobulk <- celltypes_count(MNP_matrix_pseudobulk_clean)
#117
celltypes_CD14_pseudobulk <- celltypes_count(CD14_matrix_pseudobulk_clean)
#117
celltypes_CD4_fine <- celltypes_count(CD4_matrix_fine_clean)
#smallest subset is 2
celltypes_CD8_fine <- celltypes_count(CD8_matrix_fine_clean)
celltypes_cMono_fine <- celltypes_count(cMono_matrix_fine_clean)
#2
celltypes_MNP_fine <- celltypes_count(MNP_matrix_fine_clean)
#2

#basically I can't use the fine subsets as there are not enough replicates of the different cell types

#sample the matrix to get 20 columns of each cell type 
#create function to clean matrix and downsample

matrix_sample <- function (matrix){
    # convert "." to "_" in colnames
    colnames(matrix) <- colnames(matrix) %>% str_replace_all("\\.(\\D)", "_\\1")
    print(head(colnames(matrix)))
    
    # remove any columns with "NA"
    print("dimensions before NA removal")
    print(dim(matrix))
    matrix_2 <- matrix %>% select(!starts_with("NA"))
    print("dimensions after NA removal")
    print(dim(matrix_2))
    
    # remove any rows that sum to zero
    matrix_2 <- matrix_2 %>% slice (-which(rowSums(matrix_2) == 0))
    print("dimensions following removal of rows that sum to zero")
    print(dim(matrix_2))

    #transpose
    matrix_2t <- t(matrix_2) %>% as.data.frame()
    #add column with celltype
    names <- rownames(matrix_2t) %>% 
        str_split("\\.", simplify = TRUE) %>%
        as.data.frame() %>%
        select(1) %>%
        rename("celltype" = "V1")
    matrix_2t <- matrix_2t %>% 
        mutate(celltype = names, .before = 1) %>% # add column with cell type
        group_by(celltype) %>% # group by cell type
        slice_sample(n=20, replace = FALSE)#select 20 random samples from each group
    
    #transpose back
    matrix_2tt <- t(matrix_2t) %>%
        as.data.frame()
    
    # convert row "celltype" to colnames
    colnames(matrix_2tt) <- matrix_2tt [1,]
    matrix_2tt <- matrix_2tt[-1,]
    #print(matrix_2tt[1:10, 1:10])
    
    #convert rownames to first column
    matrix_final <- matrix_2tt %>% rownames_to_column(var = "GeneSymbol")
    print(matrix_final[1:5,1:5])
    print("dimensions following downsampling")
    print(dim(matrix_final))
    
    return(matrix_final)
}

test_sample <- matrix_sample(CD8_matrix_pseudobulk_test)
CD8_matrix_pseudobulk_sample <- matrix_sample(CD8_matrix_pseudobulk)
CD4_matrix_pseudobulk_sample <- matrix_sample(CD4_matrix_pseudobulk)
cMono_matrix_pseudobulk_sample <- matrix_sample(cMono_matrix_pseudobulk)
MNP_matrix_pseudobulk_sample <- matrix_sample(MNP_matrix_pseudobulk)
CD14_matrix_pseudobulk_sample <- matrix_sample(CD14_matrix_pseudobulk)
cMono_matrix_fine_sample <- matrix_sample(cMono_matrix_fine)

#Remove celltype cMono_MKI67_RRM2_cyc from cMono_matrix_fine as there are only 2 columns
cMono_matrix_fine_sample <- cMono_matrix_fine_sample %>% select(!starts_with("cMono_MKI67"))


dir.create("downsampled_matrices", showWarnings = FALSE)

write.table (CD8_matrix_pseudobulk_sample, "./downsampled_matrices/CD8_matrix_pseudobulk_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (CD4_matrix_pseudobulk_sample, "./downsampled_matrices/CD4_matrix_pseudobulk_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (cMono_matrix_pseudobulk_sample, "./downsampled_matrices/cMono_matrix_pseudobulk_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (MNP_matrix_pseudobulk_sample, "./downsampled_matrices/MNP_matrix_pseudobulk_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (CD14_matrix_pseudobulk_sample, "./downsampled_matrices/CD14_matrix_pseudobulk_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table (cMono_matrix_fine_sample, "./downsampled_matrices/cMono_matrix_fine_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
