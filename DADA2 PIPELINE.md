This markdown contains the code used to extract amplicon sequence variants from 16S rRNA sequencing data available as multiplexed fastq files with a separate index for the barcodes.
   

Date:       14.07.18
Function:   Extract exact sequence variants from 16S rRNA sequencing data  
Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk

--------------------------------------------------------------------------

Usage:
======

Prior to using this script, the quality of the fastq files should be examined, primers removed and adapters trimmed. The cleaning of the reads (trimming bad quality bases, removing PhiX) and the forward/reverse reads merging is not performed in the first steps of the analysis on the multiplexed files like in QIIME, but at a later stage using DADA2 in R following DADA2 recommended parameters.

The fastq files are demultiplexed into individual fastq files using QIIME scripts (split_libraries_fastq.py and split_libraries_fastq.py). All the log files are stored and their content summarised using custom python scripts

Reads in the individual fastq files are denoised using DADA2. Exact sequence variants are resolved and the taxonomy assigned using SILVA v132 database.



## ---------------------------------------------------------------------------------------------------
## HEAD OF THE SCRIPT TO RUN ON THE HPC CLUSTER 
## ---------------------------------------------------------------------------------------------------

```{bash}
#$ -S /bin/sh

#$ -cwd

#$ -q HighMemLongterm.q
#$ -l h_vmem=16G
#$ -pe smp 4


# I run qiime1.9.2 with anaconda on the cluster
module load general/miniconda/4.3.21

module load general/R/3.5.0
# R packages are installed in my home directory

# Start qiime1.
source activate qiime1
```


## ---------------------------------------------------------------------------------------------------
## DEMULTIPLEX THE FASTQ FILES BASED ON THE READ INDEX # 
## ---------------------------------------------------------------------------------------------------


All FASTQ files were demultiplexed using QIIME scripts (split_libraries_fastq.py and split_libraries_fastq.py). The split_libraries_fastq.py parameters used are recommended by the Earth Microbiome Project (http://www.earthmicrobiome.org/protocols-and-standards/initial-qiime-processing/)

The following script is used as a loop across directories, each containing the different multiplexed files for a unique sequencing run


```{bash}

cd ~/PATH/TO/DIRECTORY

for dir in `find . -maxdepth 1 -mindepth 1 -type d -printf "%f\n"`; do
  cd ${dir}
  echo ${dir}
  
# The forward reads (R1):
  split_libraries_fastq.py \
  -i *_R1_001.fastq.gz \
  -b *_I1_001.fastq.gz \
  --rev_comp_mapping_barcodes  \
  -o ~/brc_scratch/16S/batch_1_illumina/${dir}_R1/  \
  -m ~/brc_scratch/16S/${dir}.txt \
  --store_demultiplexed_fastq \
  -r 999 \
  -n 999 \
  -q 0  \
  -p 0.0001 \
  -v
  
# The reverse reads (R2):
  split_libraries_fastq.py \
  -i *_R2_001.fastq.gz \
  -b *_I1_001.fastq.gz \
  --rev_comp_mapping_barcodes  \
  -o ~/brc_scratch/16S/batch_1_illumina/${dir}_R2/  \
  -m ~/brc_scratch/16S/${dir}.txt \
  --store_demultiplexed_fastq \
  -r 999 \
  -n 999 \
  -q 0  \
  -p 0.0001 \
  -v

  cd ..

done


```

Check the number of sequences written in the log files. It should be the same for the forward (R1) and reverse (R2) reads.
If some reads have been dropped, the R1 and R2 fastq files may not be in the same order. This can be fixed using the repair.sh script of BBTools (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).


## ---------------------------------------------------------------------------------------------------
## SPLIT THE RESULT OF split_libraries_fastq.py INTO PER-SAMPLE FASTQ FILES # 
## ---------------------------------------------------------------------------------------------------

The output by split_sequence_file_on_sample_ids.py is named sample_name.fastq, and don't have R1 or R2 in the name. The following scripts add the correct names to the sample of a given sequencing run.

```{bash}
cd ~/PATH/TO/DIRECTORY
for dir in `find . -maxdepth 1 -mindepth 1 -type d -printf "%f\n"`; do

  split_sequence_file_on_sample_ids.py \
    -i ${dir}/seqs.fastq \
    --file_type fastq \
    -o fastq_${dir}

  rm -r ${dir}

done


# The read orientation should be added in the name of per-sample fastq files 
for file in *_R1/*.fastq; do   mv "$file" "${file%.fastq}_R1.fastq"; done
for file in *_R2/*.fastq; do   mv "$file" "${file%.fastq}_R2.fastq"; done

source deactivate # QIIME is not used in the following steps
```

The forward and reverse fastq files are grouped in similar folders (the content of SEQUENCING_RUN_R1 and SEQUENCING_RUN_R2 should be transferred in a folder SEQUENCING RUN. I do that manually since it allows me to inspect the output of the previous steps. 

FASTQ files with a very low number of reads should be discarded as they might cause an error in DADA2.
A threshold of 10,000 seems appropriate and can be used. 



## ---------------------------------------------------------------------------------------------------
## DENOISE THE READS AND RESOLVE EXACT SEQUENCE VARIANTS FOR EACH FASTQ USING DADA2
## ---------------------------------------------------------------------------------------------------

The following R script can be ran in order to loop across the different directories and process the reads using DADA2.

https://benjjneb.github.io/dada2/faq.html

Each Miseq run was processed separately in order to account for the differences in sequencing quality in the error modelling step.

```{bash}

R CMD BATCH ~/PATH/TO/DIRECTORY/Dada2_Script_Twins_Loop.R

```

```{R Dada2_Script_Twins_Loop}


require(dada2)
require(ggplot2)

setwd("~/PATH/TO/DIRECTORY")
directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
actual_directory <- getwd()

for(i in directories ){for(j in i){
  print(j)


###########  PREPARE THE DATA   #######################

setwd(j)
extension <- substring(j, 2)  

path <- paste(actual_directory, extension, sep = "", collapse = NULL)

list.files(path) 


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

# Extract sample names (format SAMPLENAME_R1.fastq)
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


###########  EXAMINE READ QUALITY PROFILES   #######################

plotQualityProfile_F <- plotQualityProfile(fnFs[1:16])
ggsave('plotQualityProfile_F', plotQualityProfile_F, device="jpg")

plotQualityProfile_R <- plotQualityProfile(fnRs[1:16])
ggsave('plotQualityProfile_R', plotQualityProfile_R, device="jpg")



###########  FILTER THE FORWARD AND REVERSE READS   #######################

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Forward and reverse reads were trimmed by 5 and 10 bp at their 3’ end, respectively. We also trimmed 10 bases at the 5’ side of reads. 
# This parameter depends on the quality of the FASTQ files

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,240), trimLeft = 10,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE, matchIDs=TRUE) 
head(out)



###########  Learn the Error Rates   #######################

errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors_F <- plotErrors(errF, nominalQ=TRUE)
ggsave("plotErrors_F.jpg", plotErrors_F, device="jpg")
dada2:::checkConvergence(errF)

errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors_R <- plotErrors(errR, nominalQ=TRUE)
ggsave("plotErrors_R.jpg", plotErrors_R, device="jpg")
dada2:::checkConvergence(errR)




########### Dereplicate the filtered fastq file  ###########

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


########### Infer the sequence variants in each sample #####

# By default, dada2 processes each sample separately. In this case, we pooled information across samples since it is shown that it can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]     # Inspecting the dada-class object returned by dada
dadaRs[[1]]


########  Merge the denoised forward and reverse reads #####

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


######## Construct sequence table #########################

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))   # # Inspect distribution of sequence lengths


########  Remove chimeric sequences  #######################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


########  Track reads through the pipeline###############

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

rm(derepFs, derepRs)  # Remove the dereplicated files which are very large
unlink("filtered", recursive=TRUE) # Remove the filtered file which are not used anymore


# Assign the name of the sequencing run to the ASV table in order to use all ASV tables in the following analysis
imagename <- substring(j, 3)
assign(paste(imagename, "seqtab.nochim", sep = '_'), seqtab.nochim) 

# Save the progression
save.image(paste(imagename, 'image', sep = "_", collapse = NULL))

# Return to the main 16S directory before looping to another folder containing the sequencing runs
setwd("PATH/TO/DIRECTORY") 

}}

# Save the final result
save.image("~/PATH/TO/DIRECTORY/Dada2_Seqtables.R")
```

At this stage, the 16S folder is supposed to contain subfolders with each sequencing run details including the FASTQ files, the logfiles, the plots generated during the processing by DADA2, and one R image containing a table 'seq_tab


```{R}
list = ls(pattern = "seqtab.nochim")
final_seqtab <- mergeSequenceTables(paste(list))

# Assign the taxonomy up to the genus level using the Silva 132 databases available from DADA2 website
tt <- assignTaxonomy(final_seqtab, "~/PATH/TO/DIRECTORY/silva_nr_v132_train_set.fa", tryRC=TRUE)

# Add the species
taxa <- addSpecies(tt, "~/PATH/TO/DIRECTORY/silva_species_assignment_v132.fa", allowMultiple=TRUE, tryRC=TRUE)

# Write the taxa and sequence tables as .txt tables
write.table(taxa, "taxa.txt", sep ='\t')
write.table(final_seqtab, "final_seqtab.txt", sep ='\t')

```

The analysis with DADA2 is now complete. 

--------------------------------------------------------------------------

Author:     Robin Mesnage, King's College London, 2018
Contact:    robin.mesnage@kcl.ac.uk
