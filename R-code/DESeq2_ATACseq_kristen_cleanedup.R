#################################################
###clear global environment
rm(list = ls())

### LOAD REQUIRED PACKAGES
BiocManager::install("ggpubr")

library("DESeq2")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("ggpubr")
library("tidyr")

#################################################
### PREP INPUT TABLE

###code for ATAC-seq DEseq2

### set your working dir
setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/atac_seq/")

### read in the tab-delimited count table
countdata_B_BD_S_SD <- read.table("bamCountsWithinPeaks_all.tab", sep="\t", header = FALSE)

###metadata (needed for batched DeSeq2)
ATACmetadata <- read.table("ATACmetadata_DeSeq2_all.txt", header = TRUE)
model.matrix(~batch+genotype+treatment+genotype:treatment, data = ATACmetadata)

### check the table looks correct
head(countdata_B_BD_S_SD)

### create a new column which is the chr and coordinates combined
countdata_B_BD_S_SD$chr <- paste(countdata_B_BD_S_SD$V1,countdata_B_BD_S_SD$V2, sep=":")
countdata_B_BD_S_SD$coord <- paste(countdata_B_BD_S_SD$chr,countdata_B_BD_S_SD$V3, sep="-")

### set countdata$coord to be the rownames
rownames(countdata_B_BD_S_SD) <- countdata_B_BD_S_SD$coord
head(countdata_B_BD_S_SD)

### remove the first four columns (chr, start, stop, regionlabel)
### and the last two columns (coord, chr) since we don't need them anymore
### retaining only the counts
countdata_B_BD_S_SD <- countdata_B_BD_S_SD[, c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)] 
head(countdata_B_BD_S_SD)

### add bam names as column names
colnames(countdata_B_BD_S_SD) <- c("B1","B2","B3","B1D","B2D", "B3D",
                                   "S1","S2","S3","S1D","S2D", "S3D")
head(countdata_B_BD_S_SD)

### convert table to matrix format
countdata_B_BD_S_SD <- as.matrix(countdata_B_BD_S_SD)
head(countdata_B_BD_S_SD)


#################################################
### RUN DESEQ2

### construct a DESeqDataSet
dds2 = DESeqDataSetFromMatrix(countData = countdata_B_BD_S_SD,
                              colData = ATACmetadata,
                              design = ~batch+genotype+treatment+genotype:treatment)


### run deseq2
dds2 <- DESeq(dds2)
resultsNames(dds2)

### call results
contrast1_SDvsBD = results(dds2, contrast=list(c("genotype_II_vs_I", "genotypeII.treatmentDMXAA"))) #SDvsBD

### omit rows with counts of "N/A" 
contrast1_SDvsBD <- na.omit(contrast1_SDvsBD)
head(contrast1_SDvsBD)

### sort results by ascending adjusted pvalue
contrast1_SDvsBD <- contrast1_SDvsBD[order(contrast1_SDvsBD$padj), ]
head(contrast1_SDvsBD)

### report the number of rows with an adjusted pvalue less than 0.05
table(contrast1_SDvsBD$padj<0.05)



### Merge your results table with the normalized count data table. 
###To clarify, this adds additional columns to your results table to reflect the normalized 
###counts for your samples - it does not change the values in your results table in any way.
c1_resdata_SDvsBD <- merge(as.data.frame(contrast1_SDvsBD), 
                           as.data.frame(counts(dds2, normalized = TRUE)), by = "row.names", sort = FALSE)
names(c1_resdata_SDvsBD)[1] <- "ATAC_coord"



### unfiltered tables
c1_resdata_SDvsBD_nf = c1_resdata_SDvsBD
write.table(c1_resdata_SDvsBD_nf, file="ATAC_peaks.txt", row.names = FALSE, sep ="\t")

### pre-deseq2 atac peaks
#write.table(countdata_B_BD_S_SD, file="ATAC_peaks_raw.txt", row.names = TRUE, sep ="\t") #edited final version (columns) in excel




################################################################################
### Code below for plotting


setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/atac_seq/")
### now read in table of nearest gene ATAC peaks
atac_nearest_gene <- read.table("closest_gene_annotated_ATAC_peaks_deseq_SDvsBD.txt",
                                sep = "\t",
                                header = T)
atac_nearest_gene <-as.data.frame(atac_nearest_gene)


setwd("/Users/adamdziulko/Documents/Chuong\ Lab/Projects/t_kristen_mouseSP140/cutrun")
### read in table of nearest gene cutrun peaks
cutrun_nearest_gene <- read.table("pvalsort_closest_altered_gene_annotated_cutrun_peaks.bed",
                                  sep = "\t",
                                  header = T)
cutrun_nearest_gene <- as.data.frame(cutrun_nearest_gene)
cutrun_nearest_gene <- subset(cutrun_nearest_gene, pValue > 3)

### extract genes from cutrun file
cnr_atac_names <- as.data.frame(cutrun_nearest_gene$gene)
cnr_atac_names <- unique.data.frame(cnr_atac_names)
colnames(cnr_atac_names) <- c("gene")
### merge cutrun genes with ATAC gene list to get overlapping nearest genes between the two datasets
cnr_atacseq_nearestgene_overlap <- merge(atac_nearest_gene, cnr_atac_names, by = "gene")

### remove NAs
cnr_atacseq_nearestgene_overlap <- na.omit(cnr_atacseq_nearestgene_overlap)

### make cutoffs (from atac peak dataset) for significant genes. These will be labeled in volcano plot
value4_cnr_atac_nearestgene = subset(cnr_atacseq_nearestgene_overlap, padj< 5e-2 & log2FoldChange > 1)
nrow(value4_cnr_atac_nearestgene)
value5_cnr_atac_nearestgene = subset(cnr_atacseq_nearestgene_overlap, padj< 5e-2 & log2FoldChange < -1)
nrow(value5_cnr_atac_nearestgene)

### make files from atac/cutrun overlap for great analysis
atac_great = subset(cnr_atacseq_nearestgene_overlap, padj< 5e-2 & log2FoldChange > .7)
atac_great = atac_great[,c(-1,-5:-27)]
nearestgene_great = subset(cnr_atacseq_nearestgene_overlap, padj< 5e-2 & log2FoldChange > .7)
nearestgene_great = nearestgene_great[,c(-1:-22,-26,-27)]

### write tables to input into GREAT
#write.table(atac_great, file="great_ATAC_cutrun_nearestgene_overlap.bed", row.names = FALSE, sep ="\t", col.names = FALSE)
#write.table(nearestgene_great, file="great_atac_cutrun_nearestGENE_overlap2.bed", row.names = FALSE, sep ="\t", col.names = FALSE)



#_______________________________________________________________________
### Prepare your volcano plot ()
atac_cutrun_updated_nearestgene_volcano = ggplot(cnr_atacseq_nearestgene_overlap, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("atac/cutrun..nearest gene") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10 adjusted pvalue", x = "log2 fold change") +
  ## Set all dots color to grey
  geom_point(data=cnr_atacseq_nearestgene_overlap, colour = "grey") + 
  geom_point(data = value4_cnr_atac_nearestgene, colour = "red", size = 3) +
  geom_point(data = value5_cnr_atac_nearestgene, colour = "blue", size = 3) +
  geom_text_repel(data = value4_cnr_atac_nearestgene, mapping = aes(log2FoldChange, -log10(padj), label = gene ), size = 5, force = 8) +
  geom_text_repel(data = value5_cnr_atac_nearestgene, mapping = aes(log2FoldChange, -log10(padj), label = gene), size = 5, force = 8) +
  theme(axis.text.x = element_text(face="bold", size=22),
        axis.text.y = element_text(face="bold", size=22)) +
  theme(axis.title=element_text(size=28), plot.title = element_text(size=28, face = "bold")) 
atac_cutrun_updated_nearestgene_volcano
figure = ggarrange(atac_cutrun_updated_nearestgene_volcano, ncol = 1, nrow = 1)
pdf(file = paste0("atac_cutrun_nearest_gene_volcano_Kristen_batched_4", ".pdf"), width = 8, height = 10)
figure
dev.off()

##################################

