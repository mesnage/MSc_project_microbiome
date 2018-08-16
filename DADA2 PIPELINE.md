This markdown contains the code used to extract amplicon sequence variants from 16S rRNA sequencing data available as multiplexed fastq files with a separate index for the barcodes.
   

Date:       14.07.18
Function:   Extract exact sequence variants from 16S rRNA sequencing data  
Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk

--------------------------------------------------------------------------

Prior to using this script, the quality of the fastq files should be examined, primers removed and adapters trimmed. The cleaning of the reads (trimming bad quality bases, removing PhiX) and the forward/reverse reads merging is not performed in the first steps of the analysis on the multiplexed files like in QIIME, but at a later stage using DADA2 functions in R.

The fastq files are demultiplexed into individual fastq files using QIIME scripts (split_libraries_fastq.py and split_libraries_fastq.py). 
Reads in the individual fastq files are denoised using DADA2. Exact sequence variants are resolved and the taxonomy assigned using SILVA v132 database.


### HEAD OF THE SCRIPT TO RUN ON THE HPC CLUSTER 
I perform my data analysis on Rosalind High Performance Compute (HPC) cluster, part of BRC at Guy's and St Thomas' NHS Foundation Trust, using 4 threads with a maximum RAM of 64 GB. 
However, this workflow can also be done on a laptop since DADA2 is processing samples independently which means that the memory required is not inflating with the number of samples like in some other pipelines (e.g. Mothur). For instance, it is possible to process 100 sample, each containing around 50,000 reads, overnight on a MacBook Pro 2016 with 16 Go of RAM by multithreading the analysis in DADA2 functions.

I recommend to use a recent version of DADA2 (> 1.7) because taxonomic assignment at the species level has been improved in the more recent versions.


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




### DEMULTIPLEX THE FASTQ FILES BASED ON THE READ INDEX # 

All FASTQ files were demultiplexed using QIIME scripts (split_libraries_fastq.py and split_libraries_fastq.py). 
The split_libraries_fastq.py parameters used are recommended by the Earth Microbiome Project (http://www.earthmicrobiome.org/protocols-and-standards/initial-qiime-processing/)

The following script is used as a loop across directories, each containing the different multiplexed files for a unique sequencing run. This script can be modified to run into a single folder on a single MiSeq run.


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

IMPORTANT: Check the number of sequences written in the log files. It should be the same for the forward (R1) and reverse (R2) reads.
If some reads have been dropped, the R1 and R2 fastq files may not be in the same order. This can be fixed using the repair.sh script of BBTools (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

When running this pipeline on a large dataset, it can be interesting to group the content of all log files together in order to extract summary statistics.
I have written a short Python script that can be run to perform summary statistics on a series of split_libraries_fastq logfiles. The log files produced by split_libraries_fastq.py should be concatenated as a single file.
This can be done at the command line using 'cat split_library_log*.txt > Concat_Demux_Log.txt'

```{Python}
with open("/Users/robin/Logfiles/Split_library_log/Concat_Demux_Log.txt", 'r') as f:
    alignments_log = f.read().splitlines()

import re

# Statistics on the total number of reads processed
p = re.compile(r'Total number of input sequences: (.*)')

total_inputs =[]

for line in alignments_log:
    match = p.search(line)
    if match:
        total_inputs.append(float(match.group(1)))

print('Total number of input sequences:')
print("     Total is", sum(total_inputs))
print("     Average is", sum(total_inputs)/len(total_inputs))
print("     Minimum is", min(total_inputs))
print("     Maximum is", max(total_inputs))


# Statistics on the number of reads written
p = re.compile(r'Total number seqs written *(.*)')

reads_written =[]

for line in alignments_log:
    match = p.search(line)
    if match:
        reads_written.append(float(match.group(1)))

print('')
print('Total number seqs written:')
print("     Total is", sum(reads_written))
print('     Average is', sum(reads_written)/len(reads_written))
print("     minimum is", min(reads_written))
print("     maximum is", max(reads_written))


# Statistics on the number of reads per sample
q = re.compile(r'	(.+?)$')

reads_per_sample =[]
for line in alignments_log:
    match = q.search(line)
    if match:
        reads_per_sample.append(float(match.group(1)))

discarded_samples = [a for a in reads_per_sample if (a <= 10000)]
retained_samples = [a for a in reads_per_sample if (a >= 10000)]

print('')
print('There are', len(discarded_samples), 'fastq with less than 10000 reads')
print('There are', len(retained_samples), 'fastq with more than 10000 reads')

print('')
print('The average number of reads per fastq is', sum(retained_samples)/len(retained_samples))
```


### SPLIT THE RESULT OF split_libraries_fastq.py INTO PER-SAMPLE FASTQ FILES # 

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

FASTQ files with a very low number of reads should be discarded as they might cause an error in DADA2. A threshold of 10,000 seems appropriate and this is what I used. 
Although it can be possible to extract meaningful information with a lower number of reads per sample, the user should be careful seems large differences in the number of reads per sample could introduce some biases when read counts are normalised to a total sum unit.




### FIND EXACT SEQUENCE VARIANTS FOR EACH FASTQ USING DADA2

The following R script can be ran in order to loop across the different directories and process the reads using DADA2.

https://benjjneb.github.io/dada2/faq.html

Each Miseq run was processed separately in order to account for the differences in sequencing quality in the error modelling step.

```{bash}

R CMD BATCH ~/PATH/TO/DIRECTORY/Dada2_Script_Twins_Loop.R

```

```{R Dada2_Script_Twins_Loop}

# THIS IS THE Dada2_Script_Twins_Loop.R FILE

require(dada2)
require(ggplot2) 

# SET A LOOP ACROSS ALL DIRECTORIES (each directory contains FASTQ files of a unique sequencing run with the forward reads '_R1.fastq' and the reverse reads '_R2.fastq'.
setwd("~/PATH/TO/DIRECTORY")
directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
actual_directory <- getwd()

for(i in directories ){for(j in i){
  print(j)


# ENTER THE DIRECTORY TO PREPARE THE DATA 

setwd(j)
extension <- substring(j, 2)  

path <- paste(actual_directory, extension, sep = "", collapse = NULL)

# Sort ensures forward/reverse reads are in same order with the given pattern
forward_reads <- sort(list.files(path, pattern="_R1.fastq"))
reverse_reads <- sort(list.files(path, pattern="_R2.fastq"))

# Extract sample names (format SAMPLENAME_R1.fastq for the forward reads for example)
sample.names <- sapply(strsplit(forward_reads, "_"), `[`, 1)

# Specify the full path to the forward_reads and reverse_reads
forward_reads <- file.path(path, forward_reads)
reverse_reads <- file.path(path, reverse_reads)


# EXAMINE READ QUALITY PROFILES AND SAVE THE PLOTS AS A jpg FILE. 
# This step is not necessary and is more a sanity control. The quality profile of the reads should be thoroughly inspected before starting this analysis using a tool like FASTQC

plotQualityProfile_F <- plotQualityProfile(forward_reads[1:16])
ggsave('plotQualityProfile_F', plotQualityProfile_F, device="jpg")

plotQualityProfile_R <- plotQualityProfile(reverse_reads[1:16])
ggsave('plotQualityProfile_R', plotQualityProfile_R, device="jpg")



# FILTER THE FORWARD AND REVERSE READS
# This is one of the most important step. The presence of bad quality bases will cause problem to lear the error rates. A too stringent trimming will make the reads to short to assign their taxonomy.
# It should also be kept in mind that the reads are not merged yet, and that a minimum overlap of 15-20 bases between the forward and the reverse read is generally recommended to complete the merging step successfully.

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtered_forward_reads <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtered_reverse_reads <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Forward and reverse reads were trimmed by 5 and 10 bp at their 3’ end, respectively. We also trimmed 10 bases at the 5’ side of reads. 
# This parameter depends on the quality of the FASTQ files
# The details of these parameters can be found in the DADA2 pipeline tutorial https://benjjneb.github.io/dada2/tutorial.html

out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, truncLen=c(245,240), trimLeft = 10,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE, matchIDs=TRUE) 
head(out)



#  LEARN THE ERROR RATE
# DADA2 models sequencing errors to resolve the correct DNA sequence of each amplicon (Callahan et al., 2016). 
Errors are modelled using the FASTQ quality scores. In other words, when other methods consider single nucleotide differences between sequences (e.g. C ->T) as a mismatch, DADA2 uses an error model to denoise the data and resolve ‘Amplicon Sequence Variants’ (ASVs).

# Learn the error rates for the forward (errF) and reverse (errR) reads
errF <- learnErrors(filtered_forward_reads, multithread=TRUE)
errR <- learnErrors(filtered_reverse_reads, multithread=TRUE)

# Save the error rate plots
plotErrors_F <- plotErrors(errF, nominalQ=TRUE)
ggsave("plotErrors_F.jpg", plotErrors_F, device="jpg")
plotErrors_R <- plotErrors(errR, nominalQ=TRUE)
ggsave("plotErrors_R.jpg", plotErrors_R, device="jpg")

# The following option allows the user to check the error rate across the different iterations. 
# This can be useful in case the error rate modelling algorithm does not converge properly in order to check that the error rate is decreasing with increasing base-calling quality, and also decreasing after several iterations.
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)


# DEREPLICATE THE FILTERED FASTQ FILES 

derepFs <- derepFastq(filtered_forward_reads, verbose=TRUE)
derepRs <- derepFastq(filtered_reverse_reads, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# INFER THE ASV FOR EACH SAMPLE
# By default, DADA2 processes each sample separately. In this case, we pooled information across samples since it is shown that it can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples.
# Pooling samples implies that memory requirements increase with the number of samples. This 'pool' option should be deactivates if the pipeline is run on a local machine.

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

dadaFs[[1]]     # Inspecting the dada-class object returned by dada for the forward reads
dadaRs[[1]].    # Inspecting the dada-class object returned by dada for the reverse reads


# MERGE THE DENOISED FORWARD AND REVERSE READS AND CONSTRUCT SEQUENCE TABLES
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))   # # Inspect distribution of sequence lengths


# REMOVE CHIMERIC SEQUENCES

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


rm(derepFs, derepRs)  # Remove the dereplicated files which are very large
unlink("filtered", recursive=TRUE) # Remove the filtered file which are not used anymore


# Assign the name of the sequencing run to the ASV table
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

At this stage, the 16S folder is supposed to contain subfolders with each sequencing run details including the FASTQ files, the logfiles, the plots generated during the processing by DADA2, and one R image containing a table 'seq_tab'


```{R}
list = ls(pattern = "seqtab.nochim")
final_seqtab <- mergeSequenceTables(paste(list))

# Assign the taxonomy up to the genus level using the Silva 132 databases available from DADA2 website
tt <- assignTaxonomy(final_seqtab, "~/PATH/TO/DIRECTORY/silva_nr_v132_train_set.fa", tryRC=TRUE)

# Add the species
taxa <- addSpecies(tt, "~/PATH/TO/DIRECTORY/silva_species_assignment_v132.fa", allowMultiple=TRUE, tryRC=TRUE)

#To make visualization easier, we can assign a unique ID to each of our sequence variants. I will create a new category in the taxonomy table called Seq that contains the SV sequence and also a category for the unique ID.
taxa <- as.data.frame(taxa)
taxa$Seq<-rownames(taxa)
taxa$SV_ID<-paste0("SV_", seq(from=1, to=nrow(taxa), by=1))
colnames(final_seqtab)<-taxa[colnames(final_seqtab),]$SV_ID #Our table of SVs currently has sequences for names, switch it to the unique ID you created above
rownames(taxa)<-taxa$SV_ID


# Write the taxa and sequence tables as .txt tables
write.table(taxa, "taxa.txt", sep ='\t')
write.table(final_seqtab, "final_seqtab.txt", sep ='\t')

```

The analysis with DADA2 is now complete. 
The count table 'final_seqtab' contains read counts for each ASV. 
The taxa table contains ASV taxonomy up to the species level

--------------------------------------------------------------------------

Author:     Robin Mesnage, King's College London, 2018
Contact:    robin.mesnage@kcl.ac.uk
