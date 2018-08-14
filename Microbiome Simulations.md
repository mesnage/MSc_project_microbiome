This markdown contains the code used to compare the output of QIIME and DADA2 on some simulated data.

I have compared the accuracy of taxonomic assignments between QIIME (open source software package the most cited for OTU analysis) and DADA2 (the most popular method for ESV analysis) on a simulated PCR amplification of the 16S rRNA gene V3-V4 regions from 50 bacteria strains using the sequencing simulator Grinder.

These have been chosen based on the results of the study 'Subspecies in the global human gut microbiome' by Costea et al., in Mol Syst Biol (2017, 13(12):960). A difficult taxonomy assignment scenario was also added. Four species of the Bacillus cereus group were added in order to test if the pathogenic agent Bacillus anthracis can be detected using 16S rRNA sequencing and if it can be differentiated from closely related members of this taxonomic group. 


### Download the fasta nucleic acid files from the NCBI servers

```{bash}
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/825/GCF_000012825.1_ASM1282v1/GCF_000012825.1_ASM1282v1_genomic.fna.gz  # Bacteroides vulgatus ATCC 8482 NCBI:txid435590
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/162/015/GCF_000162015.1_ASM16201v1/GCF_000162015.1_ASM16201v1_genomic.fna.gz # Faecalibacterium prausnitzii A2-165 NCBI:txid411483
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/185/GCF_000146185.1_ASM14618v1/GCF_000146185.1_ASM14618v1_genomic.fna.gz # [Eubacterium] eligens ATCC 27750 NCBI:txid515620
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/065/GCF_000011065.1_ASM1106v1/GCF_000011065.1_ASM1106v1_genomic.fna.gz # Bacteroides thetaiotaomicron VPI-5482 NCBI:txid226186
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/160/055/GCF_000160055.1_ASM16005v1/GCF_000160055.1_ASM16005v1_genomic.fna.gz # Dialister invisus DSM 15470 NCBI:txid592028
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/209/955/GCA_000209955.1_ASM20995v1/GCA_000209955.1_ASM20995v1_genomic.fna.gz # [Eubacterium] rectale M104/1 (firmicutes) NCBI:txid657317
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/875/GCF_000155875.1_ASM15587v1/GCF_000155875.1_ASM15587v1_genomic.fna.gz # Coprococcus comes ATCC 27758 NCBI:txid470146
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/845/GCF_000012845.1_ASM1284v1/GCF_000012845.1_ASM1284v1_genomic.fna.gz # Parabacteroides distasonis ATCC 8503 NCBI:txid435591
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/225/GCF_000020225.1_ASM2022v1/GCF_000020225.1_ASM2022v1_genomic.fna.gz # Akkermansia muciniphila ATCC BAA-835 NCBI:txid349741
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/015/GCF_000156015.1_ASM15601v1/GCF_000156015.1_ASM15601v1_genomic.fna.gz  # Butyrivibrio crossotus DSM 2876 NCBI:txid511680

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/345/GCF_000225345.1_ASM22534v1/GCF_000225345.1_ASM22534v1_genomic.fna.gz # Roseburia hominis A2-183 NCBI:txid585394
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/153/885/GCF_000153885.1_ASM15388v1/GCF_000153885.1_ASM15388v1_genomic.fna.gz # Eubacterium ventriosum ATCC 27560 NCBI:txid411463
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/025/215/GCF_001025215.1_ASM102521v1/GCF_001025215.1_ASM102521v1_genomic.fna.gz #  Bifidobacterium pseudocatenulatum DSM 20438 NCBI:txid547043
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/195/GCF_000156195.1_ASM15619v1/GCF_000156195.1_ASM15619v1_genomic.fna.gz  # Bacteroides finegoldii DSM 17565 NCBI:txid483215
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/230/275/GCF_000230275.1_ASM23027v1/GCF_000230275.1_ASM23027v1_genomic.fna.gz # Acidaminococcus intestini RyC-MR95 NCBI:txid568816
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/175/GCF_000188175.1_ASM18817v1/GCF_000188175.1_ASM18817v1_genomic.fna.gz # Phascolarctobacterium succinatutens YIT 12067 NCBI:txid62693
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/190/535/GCF_000190535.1_ASM19053v1/GCF_000190535.1_ASM19053v1_genomic.fna.gz # Odoribacter splanchnicus DSM 20712 NCBI:txid709991
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/185/705/GCF_000185705.2_Bilo_wads_3_1_6_V2/GCF_000185705.2_Bilo_wads_3_1_6_V2_genomic.fna.gz # Bilophila wadsworthia 3_1_6  NCBI:txid563192
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/005/GCF_000091005.1_ASM9100v1/GCF_000091005.1_ASM9100v1_genomic.fna.gz  # Escherichia coli O26:H11 str. 11368  NCBI:txid573235
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/915/GCF_000157915.1_ASM15791v1/GCF_000157915.1_ASM15791v1_genomic.fna.gz # Bacteroides coprophilus DSM 18228 = JCM 13818 NCBI:txid547042

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/655/GCF_000156655.1_ASM15665v1/GCF_000156655.1_ASM15665v1_genomic.fna.gz # Holdemanella biformis DSM 3989 NCBI:txid51863
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/253/335/GCF_000253335.1_ASM25333v1/GCF_000253335.1_ASM25333v1_genomic.fna.gz # Streptococcus salivarius CCHSS3 NCBI:txid1048332
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/332/875/GCF_000332875.2_ASM33287v2/GCF_000332875.2_ASM33287v2_genomic.fna.gz # Anaerostipes hadrus DSM 3319 NCBI:txid649757
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/425/GCF_000010425.1_ASM1042v1/GCF_000010425.1_ASM1042v1_genomic.fna.gz # Bifidobacterium adolescentis ATCC 15703 NCBI:txid367928
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/173/975/GCF_000173975.1_ASM17397v1/GCF_000173975.1_ASM17397v1_genomic.fna.gz # [Eubacterium] hallii DSM 3353 NCBI:txid411469
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/015/GCF_000157015.1_ASM15701v1/GCF_000157015.1_ASM15701v1_genomic.fna.gz # Bacteroides fragilis 3_1_12  NCBI:txid457424
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/975/GCF_000144975.1_ASM14497v1/GCF_000144975.1_ASM14497v1_genomic.fna.gz # Burkholderiales bacterium 1_1_47 NCBI:txid469610
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/205/GCF_000155205.1_ASM15520v1/GCF_000155205.1_ASM15520v1_genomic.fna.gz  # Ruminococcus lactaris ATCC 29176 NCBI:txid471875
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/815/GCF_000155815.1_ASM15581v1/GCF_000155815.1_ASM15581v1_genomic.fna.gz # Bacteroides eggerthii DSM 20697 NCBI:txid483216
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/675/GCF_000156675.1_ASM15667v1/GCF_000156675.1_ASM15667v1_genomic.fna.gz # Blautia hansenii DSM 20583 NCBI:txid537007

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/445/GCF_000209445.1_Lach_bact_9_1_43BFAA_V1/GCF_000209445.1_Lach_bact_9_1_43BFAA_V1_genomic.fna.gz # Lachnospiraceae bacterium 9_1_43BFAA  NCBI:txid658088
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/905/GCF_000165905.1_ASM16590v1/GCF_000165905.1_ASM16590v1_genomic.fna.gz # Bifidobacterium bifidum PRL2010 NCBI:txid702459
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/575/GCF_000210575.1_ASM21057v1/GCF_000210575.1_ASM21057v1_genomic.fna.gz # Alistipes shahii WAL 8301 NCBI:txid717959 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/210/655/GCA_000210655.1_ASM21065v1/GCA_000210655.1_ASM21065v1_genomic.fna.gz # Roseburia intestinalis XB6B4  NCBI:txid718255
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/505/GCF_000186505.1_Sutt_wads_3_1_45B_V1/GCF_000186505.1_Sutt_wads_3_1_45B_V1_genomic.fna.gz # Sutterella wadsworthensis 3_1_45B  NCBI:txid742821
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/615/GCF_000195615.1_ASM19561v1/GCF_000195615.1_ASM19561v1_genomic.fna.gz # Bacteroides clarus YIT 12056 NCBI:txid762984
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/253/395/GCF_000253395.1_ASM25339v1/GCF_000253395.1_ASM25339v1_genomic.fna.gz # Streptococcus thermophilus JIM 8232 NCBI:txid1051074
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/078/845/GCF_000007845.1_ASM784v1/GCF_000007845.1_ASM784v1_genomic.fna.gz # Bacillus_anthracis
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/832/605/GCF_000832605.1_ASM83260v1/GCF_000832605.1_ASM83260v1_genomic.fna.gz # Bacillus_mycoides
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/078/825/GCF_000007825.1_ASM782v1/GCF_000007825.1_ASM782v1_genomic.fna.gz # Bacillus_cereus_ATCC_4579

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/161/455/GCF_000161455.1_ASM16145v1/GCF_000161455.1_ASM16145v1_genomic.fna.gz # Bacillus_pseudomycoides
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/478/885/GCF_000478885.1_ASM47888v1/GCF_000478885.1_ASM47888v1_genomic.fna.gz # Adlercreutzia equolifaciens
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/422/625/GCF_000422625.1_ASM42262v1/GCF_000422625.1_ASM42262v1_genomic.fna.gz # Enterorhabdus mucosicola
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/403/355/GCF_000403355.2_Ente_caec_B7_V1/GCF_000403355.2_Ente_caec_B7_V1_genomic.fna.gz # Enterorhabdus caecimuris
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/505/GCF_000008505.1_ASM850v1/GCF_000008505.1_ASM850v1_genomic.fna.gz # Bacillus thuringiensis serovar konkukian
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/154/425/GCF_000154425.1_ASM15442v1/GCF_000154425.1_ASM15442v1_genomic.fna.gz # Coprococcus eutactus ATCC 27759
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/045/GCF_000018045.1_ASM1804v1/GCF_000018045.1_ASM1804v1_genomic.fna.gz # Citrobacter koseri ATCC BAA-895
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/405/555/GCF_001405555.1_13414_6_33/GCF_001405555.1_13414_6_33_genomic.fna.gz # Fusicatenibacter saccharivorans
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/103/815/GCF_900103815.1_IMG-taxon_2657245269_annotated_assembly/GCF_900103815.1_IMG-taxon_2657245269_annotated_assembly_genomic.fna.gz # Lachnospira pectinoschiza
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/106/845/GCF_900106845.1_PRJEB14233/GCF_900106845.1_PRJEB14233_genomic.fna.gz # Romboutsia timonensis
```


### Use Grinder to create 16S data reads

Random 16S rRNA sequencing data is then created using Grinder (Angly et al., Nucleic Acids Research, Volume 40, Issue 12, 1) with the following parameters. We don't use a length bias (-length_bias 0) because the length of the reference genomes should not affect the relative abundance of amplicons. However, we took into account the fact that some genomes can have several 16S rRNA copies so that amplicon libraries are copy biased.

The primers were chosen to cover the V3-V4 regions (forward was CCTACGGGNGGCWGCAG, reverse was GACTACHVGGGTATCTAATCC) according to Thijs et al., Front Microbiol. 2017; 8: 494). We used the default n-mer distribution: 89% of bimeras, 11% trimeras and 0.3% of quadrimera, corresponding to average published values (Quince, et al., 2011). The fraction of chimera was 10% (recommended value). The number of reads generated was 100,000. We selected a read distribution around 250+-15 bp as in (Computational Methods for Next Generation Sequencing Data Analysis, Wiley Series in Bioinformatics). 

```{bash}
# Create the FASTA file containing the forward and reverse primers and store it as 'Primers.fa'
echo -e '>forward\nCCTACGGGNGGCWGCAG\n>reverse\nGACTACHVGGGTATCTAATCC' >Primers.fa

# Unzip the genome fastq files downloaded from the previous section, and concatenate them into a single reference file
gunzip *.fna.gz
cat *genomic.fna > 16S.fna

# Create the 16S rRNA read FASTQ file by simulating a PCR with GRINDER
Grinder -reference_file 16S.fna -forward_reverse Primers.fa -length_bias 0 -unidirectional 1 -rd 550 normal 30 -mo FR -tr 100000 -ql 30 30 -fq 1 -diversity 50

```


### Taxonomic profiling with DADA2

Amplicon Sequence Variants from the FASTQ file created with Grinder are extracted using DADA2 in R. 
This section contains the code necessary to dereplicate the sequences, infer the sequence variants, and assign the taxonomy.
The taxonomy is assigned using the SILVA database (available here https://benjjneb.github.io/dada2/training.html)

```{r Dada2}

library(dada2)

# Dereplicate the fastq file
derep <- derepFastq('grinder-reads.fastq')

# Infer the sequence variants in each sample using the DADA2 algorithm 
dada <- dada(derep, err=inflateErr(tperr1, 3))

# Assign taxonomy using the SILVA database 
taxa <- assignTaxonomy(dada, "~/silva_nr_v132_train_set.fa")
genus.species <- assignSpecies(dada, "~/silva_species_assignment_v132.fa", allowMultiple = TRUE)

# Find the number of reads assigned to each ASV
require(dplyr)
genus.species <- as.data.frame(genus.species)
genus.species$seq <- row.names(genus.species)
seqtab <- makeSequenceTable(dadaFs)
seqtab <- as.data.frame(t(seqtab))
seqtab$seq <- rownames(seqtab)
abundance <- full_join( genus.species, seqtab, by ='seq')

# Calculate relative abundances and save the result table
abundance$rel_abun <- abundance$V1/100000*100
write.table(abundance, 'abundance_data2_grinder.txt', sep = '\t', row.names = FALSE)
```



### Taxonomic binning and profiling with QIIME

OTUs from the FASTQ file created with Grinder are extracted using QIIME1

```{bash}
source activate QIIME1

# Change the name of the FASTQ file because QIIME does not like the symbol '-'
mv grinder-reads.fastq grinderreads.fastq

# Convert the fastq file in a fasta file
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f grinderreads.fastq -o fastaqual

# Create a dummy mapping file otherwise QIIME does not compute the relative abundances
echo -e '#SampleID\tGroup\tInputFastaFileName\nsample1\tX\tgrinderreads.fna' >Mapping_File.txt

# Add QIIME labels for the file using the dummy mapping file 
add_qiime_labels.py -i fastaqual/ -m Mapping_File.txt -c InputFastaFileName -n 1

# Pick OTU using UCLUST
pick_open_reference_otus.py -i combined_seqs.fna -o uclust_openref/ -f

```

### Process the results of the BIOM file created with QIIME using R package phyloseq

The BIOM file created by QIIME is not directly readable with a text editor. Taxonomy and abundances values are extracted using R package Phyloseq.

```{R QIIME}
require(phyloseq)
require(dplyr)

# Read the BIOM file and extract the taxonomy and the OTU number tables
phyloseq_qiime <- import_biom(BIOMfilename = '/Users/robin/test/uclust_openref/otu_table_mc2_w_tax_no_pynast_failures.biom', 
                        parseFunction = parse_taxonomy_greengenes)
otutable_qiime <- as.data.frame(otu_table(phyloseq_qiime ))
taxa_qiime <- as.data.frame(tax_table(phyloseq_qiime ))

# Combine both taxonomy and read number tables in one table
otutable_qiime$OTU <- row.names(otutable_qiime)
taxa_qiime$OTU <- row.names(taxa_qiime)
qiime_abundance <- inner_join(otutable_qiime, taxa_qiime[c(1,2,3,4,5,6,7,11)], by = 'OTU')

# Calculate relative abundances
qiime_abundance$rel_abund <- qiime_abundance$sa1/100000*100

# Save the file as a table
write.csv(qiime_abundance, '/Users/robin/test/uclust_openref/abundance_qiime_grinder.csv', row.names = FALSE)

```

The results of the simulation including the 50 species described above are available:

https://github.com/mesnage/microbiome/blob/master/Simulation_50species_Grinder_QIIME.csv

https://github.com/mesnage/microbiome/blob/master/Simulation_50species_Grinder_DADA2.csv


Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk
--------------------------------------------------------------------------
