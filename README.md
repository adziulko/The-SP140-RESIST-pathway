# The SP140-RESIST pathway

### Scripts and files used in:

Kristen C. Witt, Adam Dziulko, Joohyun An, Filip Pekovic, Arthur Cheng, Grace Liu, Ophelia Vosshall Lee, David J. Turner, Azra Lari, Moritz M. Gaidt, Roberto Chavez,  Stefan A. Fattinger, Preethy Abraham, Harmandeep Dhaliwal, Angus Y. Lee, Dmitri I. Kotov, Laurent Coscoy, Britt Glaunsinger, Eugene Valkov, Edward B. Chuong, Russell E. Vance (2024) "The SP140-RESIST pathway regulates interferon mRNA stability and antiviral immunity".

### GEO:
Accession codes:
- RNA-seq for Ifnar–/– and Sp140–/–Ifnar–/– Bone Marrow-Derived Macrophages (BMMs) is available at GEO accession GSE269761
- RNA-seq for B6 and Sp140–/– BMMs is available at GEO accession GSE269808
- ATAC-Seq for B6 and Sp140–/– BMMs is available at GEO accession GSE269811
- CUT&RUN for HA-SP140 anti-HA is available at GEO accession GSE269315

### Programs used:
- bedtools v2.28.0 (http://bedtools.readthedocs.io/en/latest/)
- deepTools v3.0.1 (https://deeptools.readthedocs.io/en/develop/index.html)
- samtools v1.16.1 (http://www.htslib.org/)
- BBDuk/BBMap v38.05 (https://jgi.doe.gov/data-and-tools/bbtools/)
- salmon v0.13.1 (https://gensoft.pasteur.fr/docs/salmon/1.1.0/salmon.html)
- hisat2 v2.1.0 (https://github.com/DaehwanKimLab/hisat2)
- MACS2 v2.1.1 (https://pypi.org/project/MACS2/)
- BWA v0.7.15 (https://github.com/lh3/bwa)
- DESeq2 v1.38.3 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- singularity v3.1.1 (https://github.com/hpcng/singularity)

### Public databases:
- Cistrome (http://cistrome.org/)
- GREAT  (http://great.stanford.edu/public/html/)
- UCSC (https://genome.ucsc.edu/)
- GENCODE mm10 (https://www.gencodegenes.org/mouse/releases.html)

## Differential expressed genes of: 1) SP140 KO vs WT BMMs & 2) SP140/IFNAR KO vs IFNAR KO BMMs
Started with RNAseq fastq files for two datasets (from WT (B6) or KO BMMs from mice):
1) WT BMMs untreated x3, WT BMMs treated with DMXAA x3, Sp140 KO BMMS untreated x3, Sp140 KO BMMs treated with DMXAA x3 (DMXAA treatment at 100 ug/mL for 8 hours)
2) Ifnar KO BMMs untreated x3, Ifnar KO BMMs treated with DMXAA x3, Sp140 & Ifnar KO BMMS untreated x3, Sp140 & Ifnar KO BMMs treated with DMXAA x3 (DMXAA treatment at 10 ug/mL for 4 hours)

**RNAseq Workflow (for dataset #1):**
1) Trim adapters on raw fastq files using BBDuk
    - [a1_bbduk_array_PE.sbatch](RNA-seq/a1_bbduk_array_PE.sbatch)
2) Take trimmed fastq files and run through Salmon to quantify the expression of transcripts
    - [b1_salmon_PE.sbatch](RNA-seq/b1_salmon_PE.sbatch)
3) Input Salmon quantification files into R, then run DeSeq2 for differential expression (tximport version to convert transcripts to gene level) 
   - [tximport_DESeq2_kristen_cleanedup.R](R-code/tximport_DESeq2_kristen_cleanedup.R)
4) Take trimmed fastq files and run through HISAT2 to map reads to mm10 genome (used mm10 instead of mm39 since there are downstream analysis tools developed for mm10)
    - [c1_hisat2_PE_RNA_TEtran_q10SB_unstrand_mm10.sbatch](RNA-seq/c1_hisat2_PE_RNA_TEtran_q10SB_unstrand_mm10.sbatch)
5) Take BAM files from HISAT2 output and convert to bigwig file (compresses the file/map and makes a more readable format to look at genome map on UCSC genome browser)
   - [d1_deeptools_bam_to_bigwig_mm10_unstranded.sbatch](RNA-seq/d1_deeptools_bam_to_bigwig_mm10_unstranded.sbatch)

**RNAseq Workflow (for dataset #2):**
- Same as above, but use [a2_bbduk_ifnar_array_PE.sbatch](RNA-seq/a2_bbduk_ifnar_array_PE.sbatch), [b2_salmon_ifnar_PE.sbatch](RNA-seq/b2_salmon_ifnar_PE.sbatch), [c2_hisat2_ifnar_PE_RNA_TEtran_q10SB_unstrand_mm10.sbatch](RNA-seq/c2_hisat2_ifnar_PE_RNA_TEtran_q10SB_unstrand_mm10.sbatch), [d2_deeptools_ifnar_bam_to_bigwig_mm10_unstranded.sbatch](RNA-seq/d2_deeptools_ifnar_bam_to_bigwig_mm10_unstranded.sbatch) code

## SP140 chromatin binding
Started with CUT&RUN fastq files for anti-HA or IgG (from Sp140 KO BMMs all treated wtih DMXAA from mice):
1) Sp140 KO BMMs transduced with HA-SP140 x3, Sp140 KO BMMs transduced with non tagged SP140 x3, Sp140 KO BMMs IgG control (DMXAA treatment at 100 ug/mL for 8 hours)

**CUT&RUN Workflow:**
1) Trim adapters on raw fastq files using BBDuk
    - [a_bbduk_array_PE.sbatch](CUT&RUN/a_bbduk_array_PE.sbatch)
2) Take trimmed fastq files and run through BWA-MEM to map reads to mm10 genome 
    - [b_bwaMem_PE_mm10.sbatch](CUT&RUN/b_bwaMem_PE_mm10.sbatch)
3) Take BAM files from BWA-MEM output and convert to broadpeak & bedgraph files using MAC2 
   - [c1_macs2_multi_cNnbN_control_broadP_mm10.sbatch](CUT&RUN/c1_macs2_multi_cNnbN_control_broadP_mm10.sbatch)
4) Same as above, but this time without IgG control (to compare and see if there are any major differences of code; not much difference observed)
   - [c2_macs2_multi_cNnbN_NOcontrol_broadP_mm10.sbatch](CUT&RUN/c2_macs2_multi_cNnbN_NOcontrol_broadP_mm10.sbatch)
5) Take bedgraph (bdg) files from MACS2 output and convert to bigwig file (compresses the file/map and makes a more readable format to look at genome map on UCSC genome browser)
   - [d_bdg_to_bigwig.sbatch](CUT&RUN/d_bdg_to_bigwig.sbatch)

## Accessible chromatin upon Sp140 KO
Started with ATACseq fastq files (from WT (B6) or KO BMMs from mice):
1) WT BMMs untreated x3, WT BMMs treated with DMXAA x3, Sp140 KO BMMS untreated x3, Sp140 KO BMMs treated with DMXAA x3 (DMXAA treatment at 100 ug/mL for 8 hours)

**ATACseq Workflow:**
1) Trim adapters on raw fastq files using BBDuk
    - [a_bbduk_array_PE.sbatch](ATAC-seq/a_bbduk_array_PE.sbatch)
2) Take trimmed fastq files and run through BWA-MEM to map reads to mm10 genome 
    - [b_bwaMem_PE_ATAC_mm10.sbatch](ATAC-seq/b_bwaMem_PE_ATAC_mm10.sbatch)
3) Take BAM files from BWA-MEM output and convert to broadpeak & bedgraph files using MAC2 
   - [c_macs2_atac_NOcontrol_mm10.sbatch](ATAC-seq/c_macs2_atac_NOcontrol_mm10.sbatch)
4) Take BAM files from BWA-MEM output to quantify peaks using bedtools 
   - [e_bamCountTable_all_mm10_ATAC.sbatch](ATAC-seq/e_bamCountTable_all_mm10_ATAC.sbatch)
5) Input bedtools quantification files into R, then run DeSeq2 for differential peak expression 
   - [DESeq2_ATACseq_kristen_cleanedup.R](R-code/DESeq2_ATACseq_kristen_cleanedup.R)   
6) Take bedgraph (bdg) files from MACS2 output and convert to bigwig file (compresses the file/map and makes a more readable format to look at genome map on UCSC genome browser)
   - [d_bdg_to_bigwig_atac.sbatch](ATAC-seq/d_bdg_to_bigwig_atac.sbatch)
  



