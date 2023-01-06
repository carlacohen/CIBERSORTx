## Prepare counts file for bulk RNAseq data file as input for CIBERSORTx

library (tidyverse)
library(biomaRt)


## counts matrix

counts_CD4 <- read.table("./data/bulk_counts_matrices/RNA_CD4_raw_counts.txt.gz", header = TRUE, row.names = 1)
counts_CD8 <- read.table("./data/bulk_counts_matrices/RNA_CD8_raw_counts.txt.gz", header = TRUE, row.names = 1)
counts_CD14 <- read.table("./data/bulk_counts_matrices/RNA_CD14_raw_counts.txt.gz", header = TRUE, row.names = 1)

#helpfully these are just the samples from the paper

#convert ensemblIDs to gene names

ensembl <- try(useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                 mart = ensembl)

IDs_all <- str_split(rownames(counts_CD4), "\\.", simplify = TRUE)
counts_CD4$ensembl_gene_id <- IDs_all[,1]
test <- counts_CD4 %>% inner_join(mapping) %>% 
  dplyr::select(hgnc_symbol, 1:(length(counts_CD4)-1))%>% 
  distinct(hgnc_symbol, .keep_all = TRUE) 
colnames(test)[1] <- "Gene"
#remove rows with empty row names
test_3 <- test %>% slice (-which(test$Gene == ""))
write.table(test, "test.txt", quote = FALSE, sep = "\t", row.names = FALSE)

bulk_data <- function (counts, celltype){
  IDs_all <- str_split(rownames(counts), "\\.", simplify = TRUE)
  counts$ensembl_gene_id <- IDs_all[,1]
  join <- counts %>% inner_join(mapping) %>% 
    dplyr::select(hgnc_symbol, 1:(length(counts)-1))%>% 
    distinct(hgnc_symbol, .keep_all = TRUE) 
  colnames(join)[1] <- "Gene"
  #remove rows with empty row names
  out <- join %>% slice (-which(join$Gene == ""))
  return(out)
}

CD4_out <- bulk_data(counts_CD4, "CD4")
CD8_out <- bulk_data(counts_CD8, "CD8")
CD14_out <- bulk_data(counts_CD14, "CD14")

#change column names to first row
CD4_out2 <- rbind(colnames(CD4_out), CD4_out)
CD8_out2 <- rbind(colnames(CD8_out), CD8_out)
CD14_out2 <- rbind(colnames(CD14_out), CD14_out)

write.table(CD4_out2, "./data/bulk_counts_matrices/RNA_CD4_counts_hgnc_symbol.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(CD8_out2, "./data/bulk_counts_matrices/RNA_CD8_counts_hgnc_symbol.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(CD14_out2, "./data/bulk_counts_matrices/RNA_CD14_counts_hgnc_symbol.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
