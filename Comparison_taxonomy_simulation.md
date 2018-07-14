This markdown contains the code used to compare the output of different types of gut microbiome analysis on some simulated data. 

I have compared the accuracy of taxonomic assignments between QIIME (open source software package the most cited for OTU analysis) and DADA2 (the most popular method for ESV analysis) on a simulated PCR amplification of the 16S rRNA gene V3-V4 regions from 50 bacteria strains using the sequencing simulator Grinder

   
--------------------------------------------------------------------------

## Download the fasta nucleic acid files from the NCBI servers

These have been chosen based on the results of the study 'Subspecies in the global human gut microbiome' by Costea et al., in Mol Syst Biol (2017, 13(12):960). A difficult taxonomy assignment scenario was also added. Four species of the Bacillus cereus group were added in order to test if the pathogenic agent Bacillus anthracis can be detected using 16S rRNA sequencing and if it can be differentiated from closely related members of this taxonomic group. A large number of parameters can have an influence on the quality of taxonomic assignment and a comprehensive benchmarking study is beyond the scope of this preliminary analysis. 


```{bash}
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/825/GCF_000012825.1_ASM1282v1/GCF_000012825.1_ASM1282v1_genomic.fna.gz  # Bacteroides vulgatus ATCC 8482 NCBI:txid435590
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/166/035/GCF_000166035.1_ASM16603v1/GCF_000166035.1_ASM16603v1_genomic.fna.gz # Faecalibacterium cf. prausnitzii KLE1255 NCBI:txid748224
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/162/015/GCF_000162015.1_ASM16201v1/GCF_000162015.1_ASM16201v1_genomic.fna.gz # Faecalibacterium prausnitzii A2-165 NCBI:txid411483
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/185/GCF_000146185.1_ASM14618v1/GCF_000146185.1_ASM14618v1_genomic.fna.gz # [Eubacterium] eligens ATCC 27750 NCBI:txid515620
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/209/855/GCA_000209855.1_ASM20985v1/GCA_000209855.1_ASM20985v1_genomic.fna.gz # Faecalibacterium prausnitzii SL3/3 NCBI:txid657322
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/065/GCF_000011065.1_ASM1106v1/GCF_000011065.1_ASM1106v1_genomic.fna.gz # Bacteroides thetaiotaomicron VPI-5482 NCBI:txid226186
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/160/055/GCF_000160055.1_ASM16005v1/GCF_000160055.1_ASM16005v1_genomic.fna.gz # Dialister invisus DSM 15470 NCBI:txid592028
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/209/955/GCA_000209955.1_ASM20995v1/GCA_000209955.1_ASM20995v1_genomic.fna.gz # [Eubacterium] rectale M104/1 (firmicutes) NCBI:txid657317
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/209/835/GCA_000209835.1_ASM20983v1/GCA_000209835.1_ASM20983v1_genomic.fna.gz # Ruminococcus sp. SR1/5  NCBI:txid657323
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/875/GCF_000155875.1_ASM15587v1/GCF_000155875.1_ASM15587v1_genomic.fna.gz # Coprococcus comes ATCC 27758 NCBI:txid470146
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/845/GCF_000012845.1_ASM1284v1/GCF_000012845.1_ASM1284v1_genomic.fna.gz # Parabacteroides distasonis ATCC 8503 NCBI:txid435591
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/225/GCF_000020225.1_ASM2022v1/GCF_000020225.1_ASM2022v1_genomic.fna.gz # Akkermansia muciniphila ATCC BAA-835 NCBI:txid349741
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/835/GCF_000210835.1_ASM21083v1/GCF_000210835.1_ASM21083v1_genomic.fna.gz # Bacteroides fragilis 638R  NCBI:txid862962
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/015/GCF_000156015.1_ASM15601v1/GCF_000156015.1_ASM15601v1_genomic.fna.gz  # Butyrivibrio crossotus DSM 2876 NCBI:txid511680
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/345/GCF_000225345.1_ASM22534v1/GCF_000225345.1_ASM22534v1_genomic.fna.gz # Roseburia hominis A2-183 NCBI:txid585394
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/153/885/GCF_000153885.1_ASM15388v1/GCF_000153885.1_ASM15388v1_genomic.fna.gz # Eubacterium ventriosum ATCC 27560 NCBI:txid411463
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/025/215/GCF_001025215.1_ASM102521v1/GCF_001025215.1_ASM102521v1_genomic.fna.gz #  Bifidobacterium pseudocatenulatum DSM 20438 NCBI:txid547043
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1/GCF_000016525.1_ASM1652v1_genomic.fna.gz # Methanobrevibacter smithii ATCC 35061 NCBI:txid420247
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/195/GCF_000156195.1_ASM15619v1/GCF_000156195.1_ASM15619v1_genomic.fna.gz  # Bacteroides finegoldii DSM 17565 NCBI:txid483215
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/230/275/GCF_000230275.1_ASM23027v1/GCF_000230275.1_ASM23027v1_genomic.fna.gz # Acidaminococcus intestini RyC-MR95 NCBI:txid568816
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/175/GCF_000188175.1_ASM18817v1/GCF_000188175.1_ASM18817v1_genomic.fna.gz # Phascolarctobacterium succinatutens YIT 12067 NCBI:txid62693
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/190/535/GCF_000190535.1_ASM19053v1/GCF_000190535.1_ASM19053v1_genomic.fna.gz # Odoribacter splanchnicus DSM 20712 NCBI:txid709991
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/185/705/GCF_000185705.2_Bilo_wads_3_1_6_V2/GCF_000185705.2_Bilo_wads_3_1_6_V2_genomic.fna.gz # Bilophila wadsworthia 3_1_6  NCBI:txid563192
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/005/GCF_000091005.1_ASM9100v1/GCF_000091005.1_ASM9100v1_genomic.fna.gz  # Escherichia coli O26:H11 str. 11368  NCBI:txid573235
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/915/GCF_000157915.1_ASM15791v1/GCF_000157915.1_ASM15791v1_genomic.fna.gz # Bacteroides coprophilus DSM 18228 = JCM 13818 NCBI:txid547042
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/655/GCF_000156655.1_ASM15665v1/GCF_000156655.1_ASM15665v1_genomic.fna.gz # Holdemanella biformis DSM 3989 NCBI:txid51863
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/154/365/GCF_000154365.1_ASM15436v1/GCF_000154365.1_ASM15436v1_genomic.fna.gz # [Clostridium] bolteae ATCC BAA-613 NCBI:txid411902
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/253/335/GCF_000253335.1_ASM25333v1/GCF_000253335.1_ASM25333v1_genomic.fna.gz # Streptococcus salivarius CCHSS3 NCBI:txid1048332
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/332/875/GCF_000332875.2_ASM33287v2/GCF_000332875.2_ASM33287v2_genomic.fna.gz # Anaerostipes hadrus DSM 3319 NCBI:txid649757
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/425/GCF_000010425.1_ASM1042v1/GCF_000010425.1_ASM1042v1_genomic.fna.gz # Bifidobacterium adolescentis ATCC 15703 NCBI:txid367928
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/173/975/GCF_000173975.1_ASM17397v1/GCF_000173975.1_ASM17397v1_genomic.fna.gz # [Eubacterium] hallii DSM 3353 NCBI:txid411469
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/015/GCF_000157015.1_ASM15701v1/GCF_000157015.1_ASM15701v1_genomic.fna.gz # Bacteroides fragilis 3_1_12  NCBI:txid457424
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/975/GCF_000144975.1_ASM14497v1/GCF_000144975.1_ASM14497v1_genomic.fna.gz # Burkholderiales bacterium 1_1_47 NCBI:txid469610
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/205/GCF_000155205.1_ASM15520v1/GCF_000155205.1_ASM15520v1_genomic.fna.gz  # Ruminococcus lactaris ATCC 29176 NCBI:txid471875
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/155/815/GCF_000155815.1_ASM15581v1/GCF_000155815.1_ASM15581v1_genomic.fna.gz # Bacteroides eggerthii DSM 20697 NCBI:txid483216
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/156/675/GCF_000156675.1_ASM15667v1/GCF_000156675.1_ASM15667v1_genomic.fna.gz # Blautia hansenii DSM 20583 NCBI:txid537007
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/162/215/GCF_000162215.1_ASM16221v1/GCF_000162215.1_ASM16221v1_genomic.fna.gz # Bacteroides sp. D20  NCBI:txid585543
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/445/GCF_000209445.1_Lach_bact_9_1_43BFAA_V1/GCF_000209445.1_Lach_bact_9_1_43BFAA_V1_genomic.fna.gz # Lachnospiraceae bacterium 9_1_43BFAA  NCBI:txid658088
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/905/GCF_000165905.1_ASM16590v1/GCF_000165905.1_ASM16590v1_genomic.fna.gz # Bifidobacterium bifidum PRL2010 NCBI:txid702459
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/575/GCF_000210575.1_ASM21057v1/GCF_000210575.1_ASM21057v1_genomic.fna.gz # Alistipes shahii WAL 8301 NCBI:txid717959 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/210/655/GCA_000210655.1_ASM21065v1/GCA_000210655.1_ASM21065v1_genomic.fna.gz # Roseburia intestinalis XB6B4  NCBI:txid718255
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/186/505/GCF_000186505.1_Sutt_wads_3_1_45B_V1/GCF_000186505.1_Sutt_wads_3_1_45B_V1_genomic.fna.gz # Sutterella wadsworthensis 3_1_45B  NCBI:txid742821
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/615/GCF_000195615.1_ASM19561v1/GCF_000195615.1_ASM19561v1_genomic.fna.gz # Bacteroides clarus YIT 12056 NCBI:txid762984
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/183/485/GCF_000183485.2_ASM18348v1/GCF_000183485.2_ASM18348v1_genomic.fna.gz # Alistipes sp. HGB5 NCBI:txid908612
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/253/395/GCF_000253395.1_ASM25339v1/GCF_000253395.1_ASM25339v1_genomic.fna.gz # Streptococcus thermophilus JIM 8232 NCBI:txid1051074
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/078/845/GCF_000007845.1_ASM784v1/GCF_000007845.1_ASM784v1_genomic.fna.gz # Bacillus_anthracis
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/078/825/GCF_000007825.1_ASM782v1/GCF_000007825.1_ASM782v1_genomic.fna.gz # Bacillus_cereus_ATCC_4579
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/832/605/GCF_000832605.1_ASM83260v1/GCF_000832605.1_ASM83260v1_genomic.fna.gz # Bacillus_mycoides
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/161/455/GCF_000161455.1_ASM16145v1/GCF_000161455.1_ASM16145v1_genomic.fna.gz # Bacillus_pseudomycoides

```

# Create random reads with BBmap using reference genomes downloaded from NCBI for each species

Random reads are then created using the BBmap program with a coverage of 6 and parameters representative of shotgun metagenomics sequencing data. The fastq from the reads of each 50 species are then concatenated to form a reference fastq which can be used for benchmarking. 


```{bash}

randomreads.sh -Xmx4G ref=GCF_000012825.1_ASM1282v1_genomic.fna out=Bacteroides_vulgatus_ATCC8482.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_vulgatus_ATCC8482.log
randomreads.sh -Xmx4G ref=GCF_000162015.1_ASM16201v1_genomic.fna out=Faecalibacterium_prausnitzii_A2-165.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Faecalibacterium_prausnitzii_A2-165.log
randomreads.sh -Xmx4G ref=GCF_000166035.1_ASM16603v1_genomic.fna out=Faecalibacterium_prausnitzii_KLE1255.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Faecalibacterium_prausnitzii_KLE1255.log
randomreads.sh -Xmx4G ref=GCF_000146185.1_ASM14618v1_genomic.fna out=Eubacterium_eligens.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Eubacterium_eligens.log
randomreads.sh -Xmx4G ref=GCA_000209855.1_ASM20985v1_genomic.fna out=Faecalibacterium_prausnitzii_SL3_3.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Faecalibacterium_prausnitzii_SL3_3.log
randomreads.sh -Xmx4G ref=GCF_000011065.1_ASM1106v1_genomic.fna out=Bacteroides_thetaiotaomicron_VPI_548.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_thetaiotaomicron_VPI_548.log
randomreads.sh -Xmx4G ref=GCF_000160055.1_ASM16005v1_genomic.fna out=Dialister_invisus_DSM_15470.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Dialister_invisus_DSM_15470.log
randomreads.sh -Xmx4G ref=GCA_000209955.1_ASM20995v1_genomic.fna out=Eubacterium_rectale_M104_1.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Eubacterium_rectale_M104_1.log
randomreads.sh -Xmx4G ref=GCA_000209835.1_ASM20983v1_genomic.fna out=Ruminococcus_sp_SR1_5.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Ruminococcus_sp_SR1_5.log
randomreads.sh -Xmx4G ref=GCF_000155875.1_ASM15587v1_genomic.fna out=Coprococcus_comes_ATCC_27758.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Coprococcus_comes_ATCC_27758.log
randomreads.sh -Xmx4G ref=GCF_000012845.1_ASM1284v1_genomic.fna out=Parabacteroides_distasonis_ATCC_8503.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Parabacteroides_distasonis_ATCC_8503.log
randomreads.sh -Xmx4G ref=GCF_000020225.1_ASM2022v1_genomic.fna out=Akkermansia_muciniphila_ATCC_BAA_835.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Akkermansia_muciniphila_ATCC_BAA_835.log
randomreads.sh -Xmx4G ref=GCF_000210835.1_ASM21083v1_genomic.fna out=Bacteroides_fragilis_638R.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_fragilis_638R.log
randomreads.sh -Xmx4G ref=GCF_000156015.1_ASM15601v1_genomic.fna out=Butyrivibrio_crossotus_DSM_2876.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Butyrivibrio_crossotus_DSM_2876.log
randomreads.sh -Xmx4G ref=GCF_000225345.1_ASM22534v1_genomic.fna out=Roseburia_hominis_A2_183.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Roseburia_hominis_A2_183.log
randomreads.sh -Xmx4G ref=GCF_000153885.1_ASM15388v1_genomic.fna out=Eubacterium_ventriosum_ATCC_27560.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Eubacterium_ventriosum_ATCC_27560.log
randomreads.sh -Xmx4G ref=GCF_001025215.1_ASM102521v1_genomic.fna out=Bifidobacterium_pseudocatenulatum_DSM_20438.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bifidobacterium_pseudocatenulatum_DSM_20438.log
randomreads.sh -Xmx4G ref=GCF_000156195.1_ASM15619v1_genomic.fna out=Bacteroides_finegoldii_DSM_17565.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_finegoldii_DSM_17565.log
randomreads.sh -Xmx4G ref=GCF_000230275.1_ASM23027v1_genomic.fna out=Acidaminococcus_intestini_RyC-MR95.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Acidaminococcus_intestini_RyC-MR95.log
randomreads.sh -Xmx4G ref=GCF_000188175.1_ASM18817v1_genomic.fna out=Phascolarctobacterium_succinatutens_YIT_12067.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Phascolarctobacterium_succinatutens_YIT_12067.log
randomreads.sh -Xmx4G ref=GCF_000190535.1_ASM19053v1_genomic.fna out=Odoribacter_splanchnicus_DSM_20712.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Odoribacter_splanchnicus_DSM_20712.log
randomreads.sh -Xmx4G ref=GCF_000185705.2_Bilo_wads_3_1_6_V2_genomic.fna out=Bilophila_wadsworthia_3_1_6.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bilophila_wadsworthia_3_1_6.log
randomreads.sh -Xmx4G ref=GCF_000091005.1_ASM9100v1_genomic.fna out=Escherichia_coli_O26H11_str_11368.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Escherichia_coli_O26H11_str_11368.log
randomreads.sh -Xmx4G ref=GCF_000157915.1_ASM15791v1_genomic.fna out=Bacteroides_coprophilus_DSM_18228.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_coprophilus_DSM_18228.log
randomreads.sh -Xmx4G ref=GCF_000156655.1_ASM15665v1_genomic.fna out=Holdemanella_biformis_DSM_3989.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Holdemanella_biformis_DSM_3989.log
randomreads.sh -Xmx4G ref=GCF_000154365.1_ASM15436v1_genomic.fna out=Clostridium_bolteae_ATCC_BAA613.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Clostridium_bolteae_ATCC_BAA613.log
randomreads.sh -Xmx4G ref=GCF_000253335.1_ASM25333v1_genomic.fna out=Streptococcus_salivarius_CCHSS3.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Streptococcus_salivarius_CCHSS3.log
randomreads.sh -Xmx4G ref=GCF_000332875.2_ASM33287v2_genomic.fna out=Anaerostipes_hadrus_DSM_3319.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Anaerostipes_hadrus_DSM_3319.log
randomreads.sh -Xmx4G ref=GCF_000010425.1_ASM1042v1_genomic.fna out=Bifidobacterium_adolescentis_ATCC_15703.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bifidobacterium_adolescentis_ATCC_15703.log
randomreads.sh -Xmx4G ref=GCF_000173975.1_ASM17397v1_genomic.fna out=Eubacterium_hallii_DSM_3353.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Eubacterium_hallii_DSM_3353.log
randomreads.sh -Xmx4G ref=GCF_000157015.1_ASM15701v1_genomic.fna out=Bacteroides_fragilis_3_1_12.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_fragilis_3_1_12.log
randomreads.sh -Xmx4G ref=GCF_000144975.1_ASM14497v1_genomic.fna out=Burkholderiales_bacterium_1_1_47.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Burkholderiales_bacterium_1_1_47.log
randomreads.sh -Xmx4G ref=GCF_000155205.1_ASM15520v1_genomic.fna out=Ruminococcus_lactaris_ATCC_29176.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Ruminococcus_lactaris_ATCC_29176.log
randomreads.sh -Xmx4G ref=GCF_000155815.1_ASM15581v1_genomic.fna out=Bacteroides_eggerthii_DSM_20697.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_eggerthii_DSM_20697.log
randomreads.sh -Xmx4G ref=GCF_000156675.1_ASM15667v1_genomic.fna out=Blautia_hansenii_DSM_20583.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Blautia_hansenii_DSM_20583.log
randomreads.sh -Xmx4G ref=GCF_000162215.1_ASM16221v1_genomic.fna out=Bacteroides_sp_D20.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_sp_D20.log
randomreads.sh -Xmx4G ref=GCF_000209445.1_Lach_bact_9_1_43BFAA_V1_genomic.fna out=Lachnospiraceae_bacterium_9_1_43BFAA.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Lachnospiraceae_bacterium_9_1_43BFAA.log
randomreads.sh -Xmx4G ref=GCF_000165905.1_ASM16590v1_genomic.fna out=Bifidobacterium_bifidum_PRL2010.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bifidobacterium_bifidum_PRL2010.log
randomreads.sh -Xmx4G ref=GCF_000210575.1_ASM21057v1_genomic.fna out=Alistipes_shahii_WAL_8301.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Alistipes_shahii_WAL_8301.log
randomreads.sh -Xmx4G ref=GCA_000210655.1_ASM21065v1_genomic.fna out=Roseburia_intestinalis_XB6B4.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Roseburia_intestinalis_XB6B4.log
randomreads.sh -Xmx4G ref=GCF_000186505.1_Sutt_wads_3_1_45B_V1_genomic.fna out=Sutterella_wadsworthensis_3_1_45B.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Sutterella_wadsworthensis_3_1_45B.log
randomreads.sh -Xmx4G ref=GCF_000195615.1_ASM19561v1_genomic.fna out=Bacteroides_clarus_YIT_12056.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacteroides_clarus_YIT_12056.log
randomreads.sh -Xmx4G ref=GCF_000183485.2_ASM18348v1_genomic.fna out=Alistipes_sp_HGB5.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Alistipes_sp_HGB5.log
randomreads.sh -Xmx4G ref=GCF_000007845.1_ASM784v1_genomic.fna out=Bacillus_anthracis.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacillus_anthracis.log
randomreads.sh -Xmx4G ref=GCF_000007825.1_ASM782v1_genomic.fna out=Bacillus_cereus_ATCC_4579.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacillus_cereus_ATCC_4579.log
randomreads.sh -Xmx4G ref=GCF_000008505.1_ASM850v1_genomic.fna out=Bacillus_thuringiensis.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacillus_thuringiensis.log
randomreads.sh -Xmx4G ref=GCF_000832605.1_ASM83260v1_genomic.fna out=Bacillus_mycoides.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacillus_mycoides.log
randomreads.sh -Xmx4G ref=GCF_000161455.1_ASM16145v1_genomic.fna out=Bacillus_pseudomycoides.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Bacillus_pseudomycoides.log
randomreads.sh -Xmx4G ref=GCF_000253395.1_ASM25339v1_genomic.fna out=Streptococcus_thermophilus_JIM_8232.fq.gz coverage=6 minlength=100 maxlength=100 length=100 maxq=40 midq=20 minq=6 snprate=0.02 insrate=0.02 subrate=0.02 nrate=0.02 maxns=10 metagenome &> Streptococcus_thermophilus_JIM_8232.log

# Concatenate the fastq files to create the simulated shotgun metagenomics data
cat *.fq.gz > simulated_gut_metagenome.fq.gz

```


## Use Grinder to create 16S data reads

Random 16S rRNA sequencing data is then created using Grinder (Angly et al., Nucleic Acids Research, Volume 40, Issue 12, 1) with the following parameters. We don't use a length bias (-length_bias 0) because the length of the reference genomes should not affect the relative abundance of amplicons. However, we took into account the fact that some genomes can have several 16S rRNA copies so that amplicon libraries are copy biased.
The primers were chosen to cover the V3-V4 regions (forward was CCTACGGGNGGCWGCAG, reverse was GACTACHVGGGTATCTAATCC) according to Thijs et al., Front Microbiol. 2017; 8: 494). We used the default n-mer distribution: 89% of bimeras, 11% trimeras and 0.3% of quadrimera, corresponding to average published values (Quince, et al., 2011). The fraction of chimera was 10% (recommended value). The number of reads generated was 50,000. We selected a read distribution around 250+-15 bp as in (Computational Methods for Next Generation Sequencing Data Analysis, Wiley Series in Bioinformatics). Taxa abundances followed a power law distribution, considered as the best way to describe abundance distribution in prokaryotes (Gans, Wolin- sky, & Dunbar, 2005).

```{bash}

cat *genomic.fna > 16S.fna
Grinder -reference_file 16S.fna -forward_reverse Primers_forward.fa -length_bias 0 -unidirectional 1 -rd 550 normal 30 -mo FR -tr 50000 -ql 30 30 -fq 1 -diversity 50

```

## Taxonomic binning and profiling with MetaPhlAn2

```{bash}

metaphlan2.py simulated_gut_metagenome.fq.gz --input_type fastq --nproc 4 > simgut_accuracy_profile.txt 

# create a single tab-delimited table from these files
merge_metaphlan_tables.py simgut_accuracy_profile.txt > simgut_accuracy_abundance_table.txt


# Get the species only
grep -E "(s__)|(^ID)" simgut_accuracy_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > simgut_abundance_table_species.txt

# Get the number of reads for each taxonomic assignment 
metaphlan2.py simulated_gut_metagenome.fq.gz --input_type fastq --nproc 4 -t rel_ab_w_read_stats > simgut_accuracy_profile_reads.txt 

```



## Taxonomic profiling with DADA2

```{r Dada2, echo=FALSE}

require(dada2)

########### Dereplicate the fastq file  ###########

derepFs <- derepFastq('grinder-reads.fastq', verbose=TRUE)

########## Infer the sequence variants in each sample #####

dadaFs <- dada(derepFs, err=inflateErr(tperr1, 3), multithread=TRUE)

##################  Assign taxonomy  ######################

taxa <- assignTaxonomy(dadaFs, "/Users/robin/ncbi/silva_nr_v132_train_set.fa", multithread=TRUE) # Change the path 
genus.species <- assignSpecies(dadaFs, "/Users/robin/ncbi/silva_species_assignment_v132.fa") # Change the path 
genus.species <- as.data.frame(genus.species)
genus.species$seq <- row.names(genus.species)

seqtab <- makeSequenceTable(dadaFs)
seqtab <- as.data.frame(t(seqtab))
seqtab$seq <- rownames(seqtab)

require(dplyr)
abundance <- full_join( genus.species, seqtab, by ='seq')
write.table(abundance, 'abundance_data_grinder.txt', sep = '\t')

```

## Taxonomic binning and profiling with QIIME

```{bash}
source activate QIIME1

convert_fastaqual_fastq.py -c fastq_to_fastaqual -f grinder-reads.fastq -o fastaqual
pick_open_reference_otus.py -i fastaqual/grinder-reads.fna -o pick_otus  --parallel -aO 4

pick_open_reference_otus.py -i combined_seqs.fna 

```


--------------------------------------------------------------------------


DADA2 detected 55 ESVs which were assigned to 39 bacterial species. By contrast, QIIME only identified 11 species in 52 OTUs. Interestingly, the DNA sequence considered to be representative of the OTU assigned to Bacillus cereus by QIIME aligned with 100% similarity to the 16S rRNA sequence of Bacillus anthracis, which confirms that the use of OTU clustering methods may lead to erroneous results. 

The shotgun metagenomics analysis program MetaPhlan2 detected 46 unique taxa (Holdemanella biformis and Alistipes sp. HGB5 were not detected) in the simulated gut metagenome. Metaphlan2 also detected Eubacterium biforme and Ruminococcus torques although they were not present originally. Overall, metaphlan2 managed to identify most species without ambiguity, even if the identification of the complex taxa belonging to the bacillus species remains ambiguous. A unique feature of shotgun metagenomics is that it allows the identification of bacterial communities at the strain level. Metaphlan2 correctly assigned the taxonomy of 33 unique strains among the 50 strains contained in our simulated gut metagenome. Importantly, no misclassification was found at this taxonomic level perhaps suggesting that identification at the strain level could be the most reliable.


--------------------------------------------------------------------------



Copyright:  (c) Robin Mesnage, King's College London, 2018
Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk
