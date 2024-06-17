#!/bin/bash

## June 20th 2023
## Example usage:
## sbatch --array 0-x

## General settings
#SBATCH -p short # Partition or queue. In this case, short!
#SBATCH --nodes=1 # Only use a single node
#SBATCH --cpus-per-task=32 # cpus (max = 64)
#SBATCH --mem=64gb # Memory limit (max = 512 gb per node, but only request as needed)
#SBATCH --time=4:00:00 # Time limit hrs:min:sec

## Job name and output
#SBATCH --job-name=macs2
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

## Email Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

## Load modules
module load python/2.7.14/MACS/2.1.1

## Code to execute

inDir=/scratch/Users/addz6536/t_kristen_mouseSP140/atac_seq/c_bwa_bam/mm10 # sorted alignment bams
outDir=/scratch/Users/addz6536/t_kristen_mouseSP140/atac_seq/d_macs2/mm10
MACSGENOME=mm # either the size of mappable genome eg 3e9 or a species abbreviation, hs = human, mm = mouse.

queries_Sp140=($(ls ${inDir}/*.sorted.bam | xargs -n 1 basename | sed 's/.sorted.bam//g' | uniq))

# Call peaks with macs2
pwd; hostname; date

echo "Calling peaks with macs2"
echo "macs2 version: "$(macs2 --version)
echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}



macs2 callpeak --format BAMPE \
--treatment ${inDir}/${queries_Sp140[$SLURM_ARRAY_TASK_ID]}.sorted.bam  \
--name ${queries_Sp140[$SLURM_ARRAY_TASK_ID]}.atac -g ${MACSGENOME} \
--outdir ${outDir} --SPMR -B --broad

echo $(date +"[%b %d %H:%M:%S] Done")

## Explanation of arguments:
# '--treatment <file>' - input file
# '--name <name>' - file name; will be used to name output
# '-B' - whether or not to save extended fragment pileup and local lambda tracks at every bp into a bedGraph file (default = false)
# '--SPMR' - save signal per million reads for fragment pileup profiles; requires -B to be set (default = false)
# '--gsize <int>' - effective genome size; can be 1.0e+9 or 1000000000; can also use shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
# '--format <AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE>' - format of tag file (default = AUTO); use '-f BAMPE' if using paired end data
# '--outdir <outDir>' - output will be written here
# '--keep-dup' - specifies whether macs2 will retain or discard duplicates
