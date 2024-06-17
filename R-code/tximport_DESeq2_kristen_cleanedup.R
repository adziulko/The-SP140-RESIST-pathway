library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(plyr)
library(readr)
library(dplyr)
library(tibble)

#clear global environment
rm(list = ls())


#below is code used for iteraction DeSeq2 - accounting for two genotypes (B&S) and two treatments (control&DMXAA)
#I followed info from video to design code: 
#https://www.google.com/search?q=deseq2+on+2+genotypes+and+2+treatments&oq=deseq2+on+2+genotypes+and+2+treatments&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIHCAEQABiiBDIHCAIQABiiBDIHCAMQABiiBDIHCAQQABiiBNIBDjUxMTYyNzQ4MGowajE1qAIAsAIA&sourceid=chrome&ie=UTF-8#fpstate=ive&vld=cid:94ac113d,vid:X6p3E-QTcUc
################################################################################
setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/rna_seq/")

### make RNAmetadata table and make matrix of interactions
RNAmetadata <- read.table("RNAmetadata_DeSeq2Interaction_BvsS_BDvsSD.txt", header = TRUE)

### for sp140 and ifnar kos treated/untreated
RNAmetadata3 <- read.table("RNAmetadata_DeSeq2Interaction_UTSpIfvsIf_DSpIfvsIf.txt", header = TRUE)
model.matrix(~batch+genotype+treatment+genotype:treatment, data = RNAmetadata3)


### Create a files vector containing the full path to each of your quant.sf input files. The purpose of the files vector is to point to the salmon quantification (quant.sf) files. The file.path() function constructs the path to a file, with each comma (,) designating a slash (/). 
files_B_BD_S_SD <- file.path("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/rna_seq", RNAmetadata$sample)
#
files_UTIf_DIf_UTSpIf_DSpIf <- file.path("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/rna_seq", RNAmetadata3$sample)


### Give each quant.sf input file its own "sample 1", "sample 2", etc. 
names(files_B_BD_S_SD) <- paste0("sample", 1:12)
all(file.exists(files_B_BD_S_SD))
files_B_BD_S_SD
#
names(files_UTIf_DIf_UTSpIf_DSpIf) <- paste0("sample", 1:12)
all(file.exists(files_UTIf_DIf_UTSpIf_DSpIf))
files_UTIf_DIf_UTSpIf_DSpIf


### If you already have a tx2gene annotation file prepared, load it in using the readr package
tx2gene <- read_csv("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/rna_seq/tx2gene.gencode.vM18.mm10.csv")
head(tx2gene)




### Load and run tximport using your files and tx2gene objects. Depending on your sample, you might need to add the "dropInfReps = TRUE" argument for tximport to work.
txi_B_BD_S_SD <- tximport(files_B_BD_S_SD, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE, ignoreAfterBar = TRUE)
names(txi_B_BD_S_SD)
head(txi_B_BD_S_SD$counts)
#
txi_UTIf_DIf_UTSpIf_DSpIf <- tximport(files_UTIf_DIf_UTSpIf_DSpIf, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE, ignoreAfterBar = TRUE)
names(txi_UTIf_DIf_UTSpIf_DSpIf)
head(txi_UTIf_DIf_UTSpIf_DSpIf$counts)


### Start DESeq2 with txi$counts
### import ensemble table (will need later for merging with deseq output)
### if don't have this file - make it. go to ensemble.org --> biomoart (then choose species of interest. Then select 'gene stable ID version' and 'gene name'. results. save as csv)
ensemble = read.table("mm10_ensemble_GeneVersion_GeneName.txt", sep=",", header = TRUE)
colnames(ensemble) = c('Gene', 'gene_name')
ensemble$Gene <- gsub("\\..*","",as.character(ensemble$Gene))

### Convert the list into matrix format.
RNAcountdat_B_BD_S_SD <- as.matrix(txi_B_BD_S_SD$counts)
storage.mode(RNAcountdat_B_BD_S_SD) = "integer"
colnames(RNAcountdat_B_BD_S_SD) = c('B1', 'B2', 'B3', 'B1D', 'B2D', 'B3D',
                                    'S1', 'S2', 'S3', 'S1D', 'S2D', 'S3D')
#
RNAcountdat_UTIf_DIf_UTSpIf_DSpIf <- as.matrix(txi_UTIf_DIf_UTSpIf_DSpIf$counts)
storage.mode(RNAcountdat_UTIf_DIf_UTSpIf_DSpIf) = "integer"
colnames(RNAcountdat_UTIf_DIf_UTSpIf_DSpIf) = c('UTIf1', 'UTIf2', 'UTIf3', 'DIf1', 'DIf2', 'DIf3',
                                                'UTSpIf1', 'UTSpIf2', 'UTSpIf3', 'DSpIf1', 'DSpIf2', 'DSpIf3')


### run DEseq2
ddsFull_1.2 = DESeqDataSetFromMatrix(countData = RNAcountdat_B_BD_S_SD,
                                     colData = RNAmetadata,
                                     design = ~batch+genotype+treatment+genotype:treatment)
#
ddsFull_3.2 = DESeqDataSetFromMatrix(countData = RNAcountdat_UTIf_DIf_UTSpIf_DSpIf,
                                     colData = RNAmetadata3,
                                     design = ~batch+genotype+treatment+genotype:treatment)



ddsFull_1.2 = DESeq(ddsFull_1.2)
ddsFull_1.2_NAs = DESeq(ddsFull_1.2, minReplicatesForReplace=Inf) #to prevent filtering out genes with low read counts
#
ddsFull_3.2 = DESeq(ddsFull_3.2)
ddsFull_3.2_NAs = DESeq(ddsFull_3.2, minReplicatesForReplace=Inf) #to prevent filtering out genes with low read counts


### list the coefficients
resultsNames(ddsFull_1.2)
#
resultsNames(ddsFull_3.2)


### Raw/normalized counts table for submitting to GEO ----
### Raw counts
raw_counts_genes_B_BD_S_SD = counts(ddsFull_1.2_NAs)
head(raw_counts_genes_B_BD_S_SD)
write.table(raw_counts_genes_B_BD_S_SD, file="raw_counts_genes_B_BD_S_SD.txt", quote = FALSE, row.names = TRUE, sep = "\t", col.names=NA) 

### Normalized counts
normalized_counts_genes_B_BD_S_SD <- counts(ddsFull_1.2_NAs, normalized=TRUE)
tail(normalized_counts_genes_B_BD_S_SD)
write.table(normalized_counts_genes_B_BD_S_SD, file="normalized_counts_genes_B_BD_S_SD.txt", quote = FALSE, row.names = TRUE, sep = "\t", col.names=NA)


### Raw counts
raw_counts_genes_If_IfD_SpIf_SpIfD = counts(ddsFull_3.2_NAs)
head(raw_counts_genes_If_IfD_SpIf_SpIfD)
write.table(raw_counts_genes_If_IfD_SpIf_SpIfD, file="raw_counts_genes_If_IfD_SpIf_SpIfD.txt", quote = FALSE, row.names = TRUE, sep = "\t", col.names=NA) 

### Normalized counts
normalized_counts_genes_If_IfD_SpIf_SpIfD <- counts(ddsFull_3.2_NAs, normalized=TRUE)
tail(normalized_counts_genes_If_IfD_SpIf_SpIfD)
write.table(normalized_counts_genes_If_IfD_SpIf_SpIfD, file="normalized_counts_genes_If_IfD_SpIf_SpIfD.txt", quote = FALSE, row.names = TRUE, sep = "\t", col.names=NA)
#------


### make results tables
contrast1.2_BDvsB = results(ddsFull_1.2, contrast=c('treatment', 'DMXAA', 'control')) #BDvsB
contrast1.2_SDvsBD = results(ddsFull_1.2, contrast=list(c("genotype_II_vs_I", "genotypeII.treatmentDMXAA"))) #SDvsBD
contrast1.2_interaction = results(ddsFull_1.2, name='genotypeII.treatmentDMXAA')     #interaction term
contrast1.2_SvsB = results(ddsFull_1.2, contrast=c('genotype', 'II', 'I')) #SvsB
contrast1.2_SvsB_2 = results(ddsFull_1.2, contrast=c('treatmentcontrol', 'II', 'I')) #SDvsBD
#
contrast3.2_D_SpIfvsIf = results(ddsFull_3.2, contrast=list(c("genotype_II_vs_I", "genotypeII.treatmentDMXAA"))) #D_SpIfvsIf
contrast3.2_SpIf_interaction = results(ddsFull_3.2, name='genotypeII.treatmentDMXAA') #interaction term (UTIf_DIf_UTSpIf_DSpIf)
contrast3.2_SpIfvsIf = results(ddsFull_3.2, contrast=c('genotype', 'II', 'I')) #SpIfvsIf


### results tables for preventing genes with low read counts resulting in 'padj=NA' 
contrast1.2_SDvsBD_NAs = results(ddsFull_1.2_NAs, cooksCutoff = FALSE, independentFiltering=FALSE, contrast=list(c("genotype_II_vs_I", "genotypeII.treatmentDMXAA"))) #SDvsBD
#
contrast3.2_D_SpIfvsIf_NAs = results(ddsFull_3.2_NAs, cooksCutoff = FALSE, independentFiltering=FALSE, contrast=list(c("genotype_II_vs_I", "genotypeII.treatmentDMXAA"))) #D_SpIfvsIf



### Then, filter out all "N/A" entries using na.omit() and sort your res table by ascending p-value. 
### List the number of entries that have an adjusted p-value less than 0.05.
contrast1.2_BDvsB <- na.omit(contrast1.2_BDvsB)
contrast1.2_BDvsB <- contrast1.2_BDvsB[order(contrast1.2_BDvsB$padj), ]
table(contrast1.2_BDvsB$padj<0.05)
#
contrast1.2_SDvsBD <- na.omit(contrast1.2_SDvsBD)
contrast1.2_SDvsBD <- contrast1.2_SDvsBD[order(contrast1.2_SDvsBD$padj), ]
table(contrast1.2_SDvsBD$padj<0.05)
#
contrast1.2_interaction <- na.omit(contrast1.2_interaction)
contrast1.2_interaction <- contrast1.2_interaction[order(contrast1.2_interaction$padj), ]
table(contrast1.2_interaction$padj<0.05)
#
contrast1.2_SvsB <- na.omit(contrast1.2_SvsB)
contrast1.2_SvsB <- contrast1.2_SvsB[order(contrast1.2_SvsB$padj), ]
table(contrast1.2_SvsB$padj<0.05) 
#
contrast1.2_SDvsBD_NAs <- na.omit(contrast1.2_SDvsBD_NAs)
contrast1.2_SDvsBD_NAs <- contrast1.2_SDvsBD_NAs[order(contrast1.2_SDvsBD_NAs$padj), ]
table(contrast1.2_SDvsBD_NAs$padj<0.05)


contrast3.2_D_SpIfvsIf <- na.omit(contrast3.2_D_SpIfvsIf)
contrast3.2_D_SpIfvsIf <- contrast3.2_D_SpIfvsIf[order(contrast3.2_D_SpIfvsIf$padj), ]
table(contrast3.2_D_SpIfvsIf$padj<0.05)
#
contrast3.2_SpIf_interaction <- na.omit(contrast3.2_SpIf_interaction)
contrast3.2_SpIf_interaction <- contrast3.2_SpIf_interaction[order(contrast3.2_SpIf_interaction$padj), ]
table(contrast3.2_SpIf_interaction$padj<0.05)
#
contrast3.2_SpIfvsIf <- na.omit(contrast3.2_SpIfvsIf)
contrast3.2_SpIfvsIf <- contrast3.2_SpIfvsIf[order(contrast3.2_SpIfvsIf$padj), ]
table(contrast3.2_SpIfvsIf$padj<0.05)
#
contrast3.2_D_SpIfvsIf_NAs <- na.omit(contrast3.2_D_SpIfvsIf_NAs)
contrast3.2_D_SpIfvsIf_NAs <- contrast3.2_D_SpIfvsIf_NAs[order(contrast3.2_D_SpIfvsIf_NAs$padj), ]
table(contrast3.2_D_SpIfvsIf_NAs$padj<0.05)



### Merge your results table with the normalized count data table. 
### To clarify, this adds additional columns to your results table to reflect the normalized 
### counts for your samples - it does not change the values in your results table in any way.
c1.2_resdata_BDvsB <- merge(as.data.frame(contrast1.2_BDvsB), 
                            as.data.frame(counts(ddsFull_1.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1.2_resdata_BDvsB)[1] <- "Gene"
#
c1.2_resdata_SDvsBD <- merge(as.data.frame(contrast1.2_SDvsBD), 
                             as.data.frame(counts(ddsFull_1.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1.2_resdata_SDvsBD)[1] <- "Gene"
#
c1.2_resdata_interaction <- merge(as.data.frame(contrast1.2_interaction), 
                                  as.data.frame(counts(ddsFull_1.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1.2_resdata_interaction)[1] <- "Gene"
#
c1.2_resdata_SvsB <- merge(as.data.frame(contrast1.2_SvsB), 
                           as.data.frame(counts(ddsFull_1.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1.2_resdata_SvsB)[1] <- "Gene"
#
c1.2_resdata_SDvsBD_NAs <- merge(as.data.frame(contrast1.2_SDvsBD_NAs), 
                                 as.data.frame(counts(ddsFull_1.2_NAs, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1.2_resdata_SDvsBD_NAs)[1] <- "Gene"


c3.2_resdata_D_SpIfvsIf <- merge(as.data.frame(contrast3.2_D_SpIfvsIf), 
                                 as.data.frame(counts(ddsFull_3.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c3.2_resdata_D_SpIfvsIf)[1] <- "Gene"
#
c3.2_resdata_SpIf_interaction <- merge(as.data.frame(contrast3.2_SpIf_interaction), 
                                       as.data.frame(counts(ddsFull_3.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c3.2_resdata_SpIf_interaction)[1] <- "Gene"
#
c3.2_resdata_SpIfvsIf <- merge(as.data.frame(contrast3.2_SpIfvsIf), 
                               as.data.frame(counts(ddsFull_3.2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c3.2_resdata_SpIfvsIf)[1] <- "Gene"
#
c3.2_resdata_D_SpIfvsIf_NAs <- merge(as.data.frame(contrast3.2_D_SpIfvsIf_NAs), 
                                     as.data.frame(counts(ddsFull_3.2_NAs, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c3.2_resdata_D_SpIfvsIf_NAs)[1] <- "Gene"


### calculate number of hits left over after filtering padj (false discovery rate, so padj < .05 is a false discovery rate of 5%)
sum(c1.2_resdata_BDvsB$padj < 0.05, na.rm=TRUE )
sum(c1.2_resdata_SDvsBD$padj < 0.05, na.rm=TRUE )
sum(c1.2_resdata_interaction$padj < 0.05, na.rm=TRUE )
sum(c1.2_resdata_SvsB$padj < 0.05, na.rm=TRUE )
#
sum(c3.2_resdata_D_SpIfvsIf$padj < 0.05, na.rm=TRUE )
sum(c3.2_resdata_SpIf_interaction$padj < 0.05, na.rm=TRUE )
sum(c3.2_resdata_SpIfvsIf$padj < 0.05, na.rm=TRUE )


### Filter out entires that have an adjusted p-value greater than 0.05.
c1.2_resdata_BDvsB_padj <- subset(c1.2_resdata_BDvsB, padj <= 0.05)
c1.2_resdata_SDvsBD_padj <- subset(c1.2_resdata_SDvsBD, padj <= 0.05)
c1.2_resdata_interaction_padj <- subset(c1.2_resdata_interaction, padj <= 0.05)
c1.2_resdata_SvsB_padj <- subset(c1.2_resdata_SvsB, padj <= 0.05)
#
c3.2_resdata_D_SpIfvsIf_padj <- subset(c3.2_resdata_D_SpIfvsIf, padj <= 0.05)
c3.2_resdata_SpIf_interaction_padj <- subset(c3.2_resdata_SpIf_interaction, padj <= 0.05)
c3.2_resdata_SpIfvsIf_padj <- subset(c3.2_resdata_SpIfvsIf, padj <= 0.05)


### unfiltered results table
c1.2_resdata_SDvsBD_nf = c1.2_resdata_SDvsBD
c1.2_resdata_SDvsBD_nf$Gene <- gsub("\\..*","",as.character(c1.2_resdata_SDvsBD_nf$Gene))
c1.2_resdata_SDvsBD_nf_ensemble = merge(ensemble, c1.2_resdata_SDvsBD_nf, by = "Gene")
#write.table(c1.2_resdata_SDvsBD_nf_ensemble, file="all_sig_up_down_GENES_DMXAA_SPvsB6_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")
#
c2.2_resdata_SBvsBB_nf = c2.2_resdata_SBvsBB
c2.2_resdata_SBvsBB_nf$Gene <- gsub("\\..*","",as.character(c2.2_resdata_SBvsBB_nf$Gene))
c2.2_resdata_SBvsBB_nf_ensemble = merge(ensemble, c2.2_resdata_SBvsBB_nf, by = "Gene")
#write.table(c2.2_resdata_SBvsBB_nf_ensemble, file="all_sig_up_down_GENES_IFNB_SPvsB6_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")
#
c1.2_resdata_SvsB_nf = c1.2_resdata_SvsB
c1.2_resdata_SvsB_nf$Gene <- gsub("\\..*","",as.character(c1.2_resdata_SvsB_nf$Gene))
c1.2_resdata_SvsB_nf_ensemble = merge(ensemble, c1.2_resdata_SvsB_nf, by = "Gene")
#write.table(c1.2_resdata_SvsB_nf_ensemble, file="all_sig_up_down_GENES_untreated_SPvsB6_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")
#
c3.2_resdata_D_SpIfvsIf_nf = c3.2_resdata_D_SpIfvsIf
c3.2_resdata_D_SpIfvsIf_nf$Gene <- gsub("\\..*","",as.character(c3.2_resdata_D_SpIfvsIf_nf$Gene))
c3.2_resdata_D_SpIfvsIf_nf_ensemble = merge(ensemble, c3.2_resdata_D_SpIfvsIf_nf, by = "Gene")
#write.table(c3.2_resdata_D_SpIfvsIf_nf_ensemble, file="all_sig_up_down_GENES_DMXAA_SP140IFNARvsSP140_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")
#
c3.2_resdata_D_SpIfvsIf_NAs_nf = c3.2_resdata_D_SpIfvsIf_NAs
c3.2_resdata_D_SpIfvsIf_NAs_nf$Gene <- gsub("\\..*","",as.character(c3.2_resdata_D_SpIfvsIf_NAs_nf$Gene))
c3.2_resdata_D_SpIfvsIf_NAs_nf_ensemble = merge(ensemble, c3.2_resdata_D_SpIfvsIf_NAs_nf, by = "Gene")
#write.table(c3.2_resdata_D_SpIfvsIf_NAs_nf_ensemble, file="all_sig_up_down_GENES_DMXAA_SP140IFNARvsSP140_NoIndependentFilter_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")
#
c1.2_resdata_SDvsBD_NAs_nf = c1.2_resdata_SDvsBD_NAs
c1.2_resdata_SDvsBD_NAs_nf$Gene <- gsub("\\..*","",as.character(c1.2_resdata_SDvsBD_NAs_nf$Gene))
c1.2_resdata_SDvsBD_NAs_nf_ensemble = merge(ensemble, c1.2_resdata_SDvsBD_NAs_nf, by = "Gene")
#write.table(c1.2_resdata_SDvsBD_NAs_nf_ensemble, file="all_sig_up_down_GENES_DMXAA_SP140vsB6_NoIndependentFilter_unfiltered_batched_tximport.txt", row.names = FALSE, sep ="\t")


### Separate FC>0 and FC<0 entries into their own objects and sort by FC.
up_c1.2_resdata_BDvsB_padj <- subset(c1.2_resdata_BDvsB_padj, log2FoldChange >= 0)
up_c1.2_resdata_BDvsB_padj <- up_c1.2_resdata_BDvsB_padj[order(-up_c1.2_resdata_BDvsB_padj$log2FoldChange), ] #4926 genes 
up_c1.2_resdata_BDvsB_padj$Gene <- gsub("\\..*","",as.character(up_c1.2_resdata_BDvsB_padj$Gene))
#
up_c1.2_resdata_SDvsBD_padj <- subset(c1.2_resdata_SDvsBD_padj, log2FoldChange >= 0)
up_c1.2_resdata_SDvsBD_padj <- up_c1.2_resdata_SDvsBD_padj[order(-up_c1.2_resdata_SDvsBD_padj$log2FoldChange), ] #2716 genes 
up_c1.2_resdata_SDvsBD_padj$Gene <- gsub("\\..*","",as.character(up_c1.2_resdata_SDvsBD_padj$Gene))
#
up_c1.2_resdata_interaction_padj <- subset(c1.2_resdata_interaction_padj, log2FoldChange >= 0)
up_c1.2_resdata_interaction_padj <- up_c1.2_resdata_interaction_padj[order(-up_c1.2_resdata_interaction_padj$log2FoldChange), ] #1430 genes 
up_c1.2_resdata_interaction_padj$Gene <- gsub("\\..*","",as.character(up_c1.2_resdata_interaction_padj$Gene))
#
up_c1.2_resdata_SvsB_padj <- subset(c1.2_resdata_SvsB_padj, log2FoldChange >= 0)
up_c1.2_resdata_SvsB_padj <- up_c1.2_resdata_SvsB_padj[order(-up_c1.2_resdata_SvsB_padj$log2FoldChange), ] #55 genes 
up_c1.2_resdata_SvsB_padj$Gene <- gsub("\\..*","",as.character(up_c1.2_resdata_SvsB_padj$Gene))


up_c3.2_resdata_D_SpIfvsIf_padj <- subset(c3.2_resdata_D_SpIfvsIf_padj, log2FoldChange >= 0)
up_c3.2_resdata_D_SpIfvsIf_padj <- up_c3.2_resdata_D_SpIfvsIf_padj[order(-up_c3.2_resdata_D_SpIfvsIf_padj$log2FoldChange), ] #82 genes 
up_c3.2_resdata_D_SpIfvsIf_padj$Gene <- gsub("\\..*","",as.character(up_c3.2_resdata_D_SpIfvsIf_padj$Gene))
#
up_c3.2_resdata_SpIf_interaction_padj <- subset(c3.2_resdata_SpIf_interaction_padj, log2FoldChange >= 0)
up_c3.2_resdata_SpIf_interaction_padj <- up_c3.2_resdata_SpIf_interaction_padj[order(-up_c3.2_resdata_SpIf_interaction_padj$log2FoldChange), ] #11 genes 
up_c3.2_resdata_SpIf_interaction_padj$Gene <- gsub("\\..*","",as.character(up_c3.2_resdata_SpIf_interaction_padj$Gene))
#
up_c3.2_resdata_SpIfvsIf_padj <- subset(c3.2_resdata_SpIfvsIf_padj, log2FoldChange >= 0)
up_c3.2_resdata_SpIfvsIf_padj <- up_c3.2_resdata_SpIfvsIf_padj[order(-up_c3.2_resdata_SpIfvsIf_padj$log2FoldChange), ] #52 genes 
up_c3.2_resdata_SpIfvsIf_padj$Gene <- gsub("\\..*","",as.character(up_c3.2_resdata_SpIfvsIf_padj$Gene))




down_c1.2_resdata_BDvsB_padj <- subset(c1.2_resdata_BDvsB_padj, log2FoldChange <= 0)
down_c1.2_resdata_BDvsB_padj <- down_c1.2_resdata_BDvsB_padj[order(down_c1.2_resdata_BDvsB_padj$log2FoldChange), ] #4540 genes 
down_c1.2_resdata_BDvsB_padj$Gene <- gsub("\\..*","",as.character(down_c1.2_resdata_BDvsB_padj$Gene))
#
down_c1.2_resdata_SDvsBD_padj <- subset(c1.2_resdata_SDvsBD_padj, log2FoldChange <= 0)
down_c1.2_resdata_SDvsBD_padj <- down_c1.2_resdata_SDvsBD_padj[order(down_c1.2_resdata_SDvsBD_padj$log2FoldChange), ] #1374 genes 
down_c1.2_resdata_SDvsBD_padj$Gene <- gsub("\\..*","",as.character(down_c1.2_resdata_SDvsBD_padj$Gene))
#
down_c1.2_resdata_interaction_padj <- subset(c1.2_resdata_interaction_padj, log2FoldChange <= 0)
down_c1.2_resdata_interaction_padj <- down_c1.2_resdata_interaction_padj[order(down_c1.2_resdata_interaction_padj$log2FoldChange), ] #825 genes 
down_c1.2_resdata_interaction_padj$Gene <- gsub("\\..*","",as.character(down_c1.2_resdata_interaction_padj$Gene))
#
down_c1.2_resdata_SvsB_padj <- subset(c1.2_resdata_SvsB_padj, log2FoldChange <= 0)
down_c1.2_resdata_SvsB_padj <- down_c1.2_resdata_SvsB_padj[order(down_c1.2_resdata_SvsB_padj$log2FoldChange), ] #25 genes 
down_c1.2_resdata_SvsB_padj$Gene <- gsub("\\..*","",as.character(down_c1.2_resdata_SvsB_padj$Gene))


down_c3.2_resdata_D_SpIfvsIf_padj <- subset(c3.2_resdata_D_SpIfvsIf_padj, log2FoldChange <= 0)
down_c3.2_resdata_D_SpIfvsIf_padj <- down_c3.2_resdata_D_SpIfvsIf_padj[order(down_c3.2_resdata_D_SpIfvsIf_padj$log2FoldChange), ] #164 genes 
down_c3.2_resdata_D_SpIfvsIf_padj$Gene <- gsub("\\..*","",as.character(down_c3.2_resdata_D_SpIfvsIf_padj$Gene))
#
down_c3.2_resdata_SpIf_interaction_padj <- subset(c3.2_resdata_SpIf_interaction_padj, log2FoldChange <= 0)
down_c3.2_resdata_SpIf_interaction_padj <- down_c3.2_resdata_SpIf_interaction_padj[order(down_c3.2_resdata_SpIf_interaction_padj$log2FoldChange), ] #6 genes 
down_c3.2_resdata_SpIf_interaction_padj$Gene <- gsub("\\..*","",as.character(down_c3.2_resdata_SpIf_interaction_padj$Gene))
#
down_c3.2_resdata_SpIfvsIf_padj <- subset(c3.2_resdata_SpIfvsIf_padj, log2FoldChange <= 0)
down_c3.2_resdata_SpIfvsIf_padj <- down_c3.2_resdata_SpIfvsIf_padj[order(down_c3.2_resdata_SpIfvsIf_padj$log2FoldChange), ] #57 genes 
down_c3.2_resdata_SpIfvsIf_padj$Gene <- gsub("\\..*","",as.character(down_c3.2_resdata_SpIfvsIf_padj$Gene))



### merge ensemble to get gene name 
up_c1.2_ensem_BDvsB = merge(ensemble, up_c1.2_resdata_BDvsB_padj, by = "Gene")
down_c1.2_ensem_BDvsB = merge(ensemble, down_c1.2_resdata_BDvsB_padj, by = "Gene")
up_c1.2_ensem_BDvsB$data.source <- paste0("up_BDvsB", up_c1.2_ensem_BDvsB$data.source)
down_c1.2_ensem_BDvsB$data.source <- paste0("down_BDvsB", down_c1.2_ensem_BDvsB$data.source)
#
up_c1.2_ensem_SDvsBD = merge(ensemble, up_c1.2_resdata_SDvsBD_padj, by = "Gene")
down_c1.2_ensem_SDvsBD = merge(ensemble, down_c1.2_resdata_SDvsBD_padj, by = "Gene")
up_c1.2_ensem_SDvsBD$data.source <- paste0("up_SDvsBD", up_c1.2_ensem_SDvsBD$data.source)
down_c1.2_ensem_SDvsBD$data.source <- paste0("down_SDvsBD", down_c1.2_ensem_SDvsBD$data.source)
#
up_c1.2_ensem_interaction = merge(ensemble, up_c1.2_resdata_interaction_padj, by = "Gene")
down_c1.2_ensem_interaction = merge(ensemble, down_c1.2_resdata_interaction_padj, by = "Gene")
up_c1.2_ensem_interaction$data.source <- paste0("up_interaction", up_c1.2_ensem_interaction$data.source)
down_c1.2_ensem_interaction$data.source <- paste0("down_interaction", down_c1.2_ensem_interaction$data.source)
#
up_c1.2_ensem_SvsB = merge(ensemble, up_c1.2_resdata_SvsB_padj, by = "Gene")
down_c1.2_ensem_SvsB = merge(ensemble, down_c1.2_resdata_SvsB_padj, by = "Gene")
up_c1.2_ensem_SvsB$data.source <- paste0("up_SvsB", up_c1.2_ensem_SvsB$data.source)
down_c1.2_ensem_SvsB$data.source <- paste0("down_SvsB", down_c1.2_ensem_SvsB$data.source)


up_c3.2_ensem_D_SpIfvsIf = merge(ensemble, up_c3.2_resdata_D_SpIfvsIf_padj, by = "Gene")
down_c3.2_ensem_D_SpIfvsIf = merge(ensemble, down_c3.2_resdata_D_SpIfvsIf_padj, by = "Gene")
up_c3.2_ensem_D_SpIfvsIf$data.source <- paste0("up_D_SpIfvsIf", up_c3.2_ensem_D_SpIfvsIf$data.source)
down_c3.2_ensem_D_SpIfvsIf$data.source <- paste0("down_D_SpIfvsIf", down_c3.2_ensem_D_SpIfvsIf$data.source)
#
up_c3.2_ensem_SpIf_interaction = merge(ensemble, up_c3.2_resdata_SpIf_interaction_padj, by = "Gene")
down_c3.2_ensem_SpIf_interaction = merge(ensemble, down_c3.2_resdata_SpIf_interaction_padj, by = "Gene")
up_c3.2_ensem_SpIf_interaction$data.source <- paste0("up_SpIf_interaction", up_c3.2_ensem_SpIf_interaction$data.source)
down_c3.2_ensem_SpIf_interaction$data.source <- paste0("down_SpIf_interaction", down_c3.2_ensem_SpIf_interaction$data.source)


### make table containing sig up/down genes, then save tables
all_c1.2_sig_up_down_GENES_BDvsB <- rbind.fill(up_c1.2_ensem_BDvsB,
                                               down_c1.2_ensem_BDvsB)

all_c1.2_sig_up_down_GENES_SDvsBD <- rbind.fill(up_c1.2_ensem_SDvsBD,
                                                down_c1.2_ensem_SDvsBD)

all_c1.2_sig_up_down_GENES_interaction <- rbind.fill(up_c1.2_ensem_interaction,
                                                     down_c1.2_ensem_interaction)

all_c3.2_sig_up_down_GENES_D_SpIfvsIf <- rbind.fill(up_c3.2_ensem_D_SpIfvsIf,
                                                    down_c3.2_ensem_D_SpIfvsIf)

all_c3.2_sig_up_down_GENES_SpIf_interaction <- rbind.fill(up_c3.2_ensem_SpIf_interaction,
                                                          down_c3.2_ensem_SpIf_interaction)


#write.table(all_c1.2_sig_up_down_GENES_BDvsB, file="all_c1.2_sig_up_down_GENES_BDvsB_tximport.txt", row.names = FALSE, sep ="\t")
#write.table(all_c1.2_sig_up_down_GENES_SDvsBD, file="all_c1.2_sig_up_down_GENES_SDvsBD_tximport.txt", row.names = FALSE, sep ="\t")
#write.table(all_c1.2_sig_up_down_GENES_interaction, file="all_c1.2_sig_up_down_GENES_interaction_tximport.txt", row.names = FALSE, sep ="\t")
#
#write.table(all_c3.2_sig_up_down_GENES_D_SpIfvsIf, file="all_c3.2_sig_up_down_GENES_D_SpIfvsIf_tximport.txt", row.names = FALSE, sep ="\t")
#write.table(all_c3.2_sig_up_down_GENES_SpIf_interaction, file="all_c3.2_sig_up_down_GENES_SpIf_interaction_tximport.txt", row.names = FALSE, sep ="\t")








#_______________________________________________________________________________
### Volcano plot #1 (contrast3.2_D_SpIfvsIf; UTIf_DIf_UTSpIf_DSpIf; pairwise)
DT_contrast3.2_D_SpIfvsIf <- contrast3.2_D_SpIfvsIf[, c(2,6)]
DT_contrast3.2_D_SpIfvsIf = as.data.frame(DT_contrast3.2_D_SpIfvsIf)
DT_contrast3.2_D_SpIfvsIf$Gene <- rownames(DT_contrast3.2_D_SpIfvsIf)
DT_contrast3.2_D_SpIfvsIf$Gene <- gsub("\\..*","",as.character(DT_contrast3.2_D_SpIfvsIf$Gene))
DT_contrast3.2_D_SpIfvsIf = merge(ensemble, DT_contrast3.2_D_SpIfvsIf, by = "Gene")
rownames(DT_contrast3.2_D_SpIfvsIf) <- DT_contrast3.2_D_SpIfvsIf$geneName

### Define cutoffs for labeling genes
### Tweak these until you're happy with how many genes are labeled
value4_cpR = subset(DT_contrast3.2_D_SpIfvsIf, padj< 5e-2 & log2FoldChange > 1)
nrow(value4_cpR)
value5_cpR = subset(DT_contrast3.2_D_SpIfvsIf, padj< 5e-2 & log2FoldChange < -1)
nrow(value5_cpR)

### Prepare your volcano plot
cpR_UTIf_DIf_UTSpIf_DSpIf_volcano = ggplot(DT_contrast3.2_D_SpIfvsIf, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("DeSeq2; pairwise; D_SpIfvsIf") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10 adjusted pvalue", x = "log2 fold change") +
  ## Set all dots color to grey
  geom_point(data=DT_contrast3.2_D_SpIfvsIf, colour = "grey") + 
  geom_point(data = value4_cpR, colour = "red", size = 4) +
  geom_point(data = value5_cpR, colour = "blue", size = 4) +
  # Add text label for interesting outliers
  geom_text_repel(data = value4_cpR, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 5, force = 8) +
  geom_text_repel(data = value5_cpR, mapping = aes(log2FoldChange, -log10(padj), label = gene_name), size = 5, force = 8) +
  theme(axis.text.x = element_text(face="bold", size=22),
        axis.text.y = element_text(face="bold", size=22)) +
  theme(axis.title=element_text(size=28), plot.title = element_text(size=28, face = "bold")) 
cpR_UTIf_DIf_UTSpIf_DSpIf_volcano
figure = ggarrange(cpR_UTIf_DIf_UTSpIf_DSpIf_volcano, ncol = 1, nrow = 1)
pdf(file = paste0("deseq_rnaseq_batch_D_SpIfvsIf_volcano_Kristen_tximport", ".pdf"), width = 8, height = 10)
figure
dev.off()

