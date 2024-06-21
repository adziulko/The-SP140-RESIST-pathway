# The SP140-RESIST pathway

### Scripts and files used in:

Kristen Witt, Adam Dziulko, Joohyun An, Ophelia Vosshall Lee, Grace Liu, Azra Lari, Roberto Chavez, David Turner, Arthur Cheng, Dmitri Kotov, Preethy Abraham, Angus Y. Lee, Harmandeep Dhaliwal, Eugene Valkov, Laurent Coscoy, Britt Glausinger, Edward B. Chuong, Russell E. Vance (2024) "The SP140-RESIST pathway regulates interferon mRNA stability and antiviral immunity".

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

## Differential expressed genes of SP140 KO BMMs & SP140/IFNAR KO vs IFNAR KO
Started with RNAseq fastq files for two datasets (from WT (B6) or KO BMMs from mice):
1) WT BMMs untreated x3, WT BMMs treated with DMXAA x3, Sp140 KO BMMS untreated x3, Sp140 KO BMMs treated with DMXAA x3 (DMXAA treatment at 100 ug/mL for 8 hours)
2) Ifnar KO BMMs untreated x3, Ifnar KO BMMs treated with DMXAA x3, Sp140 & Ifnar KO BMMS untreated x3, Sp140 & Ifnar KO BMMs treated with DMXAA x3 (DMXAA treatment at 10 ug/mL for 4 hours)

**RNAseq Workflow:**
1) [bbduk_PE.sbatch](RNA-seq/a1_bbduk_array_PE.sbatch)
2) [hisat2_PE.sbatch](RNA-seq/b1_salmon_PE.sbatch)
3) [bam_to_bw.sbatch](RNA-seq/c1_hisat2_PE_RNA_TEtran_q10SB_unstrand_mm10.sbatch)
4) [feature_counts.sbatch](RNA-seq/d1_deeptools_bam_to_bigwig_mm10_unstranded.sbatch)
5) [deseq2_genes.R](R-code/tximport_DESeq2_kristen_cleanedup.R)

