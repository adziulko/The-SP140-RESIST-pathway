#!/bin/bash

## Script for computing bam read count table
## Date: 27 March 2024
##
## Example usage:
## bamOrder=/Shares/CL_Shared/data/atma/cohen2017_chip/bamOrder.txt \
## bed=/Shares/CL_Shared/data/atma/cohen2017_chip/ALL_peaks.narrowPeak \
## outDir=/Shares/CL_Shared/data/atma/cohen2017_chip/bamCounts \
## sbatch computeBamCountTable.q
##
## Note:
## bamOrder is a text file containing the full path and file names of all the bams, in the order they should appear as columns in the table
## You need to create the bamOrder file before running this script
## And you must include the full path to the bamOrder file, as shown above
##
## Here's what the inside of a bamOrder file could look like:
## /Shares/CL_Shared/data/atma/cohen2017_chip/file1.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file2.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file3.bam.sorted
## /Shares/CL_Shared/data/atma/cohen2017_chip/file4.bam.sorted
##

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=1-00:00
#SBATCH --mem=64GB

# Job name and output
#SBATCH -J bamCountTable
#SBATCH -o /Users/%u/slurm_out/slurm-%j.out
#SBATCH -e /Users/%u/slurm_err/slurm-%j.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# load modules
module load bedtools

# Define key variables
bamDir=/scratch/Users/addz6536/t_kristen_mouseSP140/atac_seq/c_bwa_bam/mm10
## peakDir=/Shares/CL_Shared/data/adam/c_BRD4/chip_seq/d_macs2_hs \
bamOrder=bamOrder_all.txt
bed=/scratch/Users/addz6536/t_kristen_mouseSP140/atac_seq/d_macs2/mm10/all_peaks.broadPeak
outDir=/scratch/Users/addz6536/t_kristen_mouseSP140/atac_seq/f_bamCountstable
bams=$(cat ${bamDir}/${bamOrder})

# run the thing
pwd; hostname; date

echo "Starting bedtools..."
echo $(date +"[%b %d %H:%M:%S] Computing matrix of bam read counts...")

bedtools multicov -bams ${bams} -bed ${bed} > $outDir/bamCountsWithinPeaks_all.tab

echo $(date +"[%b %d %H:%M:%S] Done!")
