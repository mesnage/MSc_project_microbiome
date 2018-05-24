This markdown contains the results of a reanalysis of the 16S rRNA sequencing data (PRJNA412959) from a study of the effect of glyphosate on the gut microbiome in Sprague Dawley rats. This data was published in the study of Lene Nørby Nielsen and collaborators in the journal Environmental Pollution (233: 364-376) in 2018. We demonstrate here that the use of Exact sequence variants is giving better results than operational taxonomic units.  

```{bash}

# Create directories for the data
cd ~
mkdir Glyphosate_Study

# List accessions numbers for this study
accessions=( SRR6127721 SRR6127722 SRR6127723  SRR6127724  SRR6127725  SRR6127726  SRR6127727  SRR6127728  SRR6127729  SRR6127730  SRR6127731  SRR6127732  SRR6127733  SRR6127734  SRR6127735  SRR6127736  SRR6127737  SRR6127738  SRR6127739  SRR6127740  SRR6127741  SRR6127742  SRR6127743  SRR6127744  SRR6127879  SRR6127720  SRR6127745  SRR6127746  SRR6127747  SRR6127748  SRR6127749  SRR6127750  SRR6127751  SRR6127752  SRR6127753  SRR6127754  SRR6127755  SRR6127756  SRR6127757  SRR6127758  SRR6127759  SRR6127760  SRR6127761  SRR6127762  SRR6127763  SRR6127764  SRR6127765  SRR6127766  SRR6127767  SRR6127768  SRR6127769  SRR6127770  SRR6127771  SRR6127772  SRR6127773  SRR6127774  SRR6127775  SRR6127776  SRR6127777  SRR6127778  SRR6127779  SRR6127780  SRR6127781  SRR6127782  SRR6127783  SRR6127784  SRR6127785  SRR6127786  SRR6127787  SRR6127788  SRR6127789  SRR6127790  SRR6127791  SRR6127792  SRR6127793  SRR6127794  SRR6127795  SRR6127796  SRR6127797  SRR6127798  SRR6127799  SRR6127800  SRR6127801  SRR6127802  SRR6127803  SRR6127804  SRR6127805  SRR6127806  SRR6127807  SRR6127808  SRR6127809  SRR6127810  SRR6127811  SRR6127812  SRR6127813  SRR6127814  SRR6127815  SRR6127816  SRR6127817  SRR6127818  SRR6127819  SRR6127820  SRR6127821  SRR6127822  SRR6127823  SRR6127824  SRR6127825  SRR6127826  SRR6127827  SRR6127828  SRR6127829  SRR6127830  SRR6127831  SRR6127832  SRR6127833  SRR6127834  SRR6127835  SRR6127836  SRR6127837  SRR6127838  SRR6127839  SRR6127840  SRR6127841  SRR6127842  SRR6127843  SRR6127844  SRR6127845  SRR6127846  SRR6127847  SRR6127848  SRR6127849  SRR6127850  SRR6127851  SRR6127852  SRR6127853  SRR6127854  SRR6127855  SRR6127856  SRR6127857  SRR6127858  SRR6127859  SRR6127860  SRR6127861  SRR6127862  SRR6127863  SRR6127864  SRR6127865  SRR6127866  SRR6127867  SRR6127868  SRR6127869  SRR6127870  SRR6127871  SRR6127872  SRR6127873  SRR6127874  SRR6127875  SRR6127876  SRR6127877  SRR6127878 )

# Fetch the fastq files corresponding to this data 
# Warning: Make sure that fastq-dum from the SRA toolkit is in the path
for i in "${accessions[@]}"
do
	fastq-dump -I --split-files $i
done

```


## Process the fastq using dada2 to create the sequence table:
```{r dada2, echo=FALSE}
setwd("~/Glyphosate_Study")

require(dada2)

path <- "~/Glyphosate_Study/" # set to the directory containing the fastq files 
list.files(path) 

# Sort ensures forward/reverse reads are in same order
fn <- sort(list.files(path, pattern="_1.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(fn, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fn <- file.path(path, fn)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fn[1:12])

```

![](images/Figure_1.tiff?raw=true)

```
filt_path <- file.path(path, "filtered") # Place filtered files in filtered subdirectory
filt <- file.path(filt_path, paste0(sample.names, "_filt.fastq"))

# Filter the forward and reverse reads:
out <- filterAndTrim(fn, filt, truncLen=10, trimLeft=15, 
                      maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE, verbose = TRUE)
head(out)

# Learn the Error Rates 
err <- learnErrors(filt, multithread=FALSE)
plotErrors(err, nominalQ=TRUE)
dada2:::checkConvergence(err)
```


