This markdown contains the code used to make a metagenome-wide association using PopPAnTe, starting from a Phyloseq object.

PopPAnTe is a program to assess association of quantitative data in familial samples. The relation between a response variable and a predictor is measured in a linear-mixed model taking into account the family relationship, which is modelled as a random effect. 

The following code is used to make the analysis in the context of an analysis on Twins with family structure and zygosity modelled as a random effect. Although this code was fit for our purposes, it can be modified and adapted to other types of metagenome-wide association studies.

More details on the use of PopPAnTe can be found in Visconti et al BMC Genomics, 2017, 18:150.


### Create the tables used by PopPAnTe from the Phyloseq object

The phyloseq object contains the OTU/ASV count table, the taxonomy table, and the metadata table. 
These table will be extracted using the following lines of code. A folder will be create in the home directory of the user and this folder will be populated with the different tables necessary to use PopPAnTe by extracting information from the Phyloseq object.

The user needs to specify the covariates and the predictors names that will be used in the linear-mixed models. These should be included in the metadata table

This metadata table should also contain the individuals ID ('Individual_id'), the Family ID ('Family_id'), the Gender ('Gender'), the Zygosity (Zygosity), as well as the sample names ('Sample_name'). They will be used to make the pedigree file used by PopPAnte to estimate family effects (more on PED files: https://gatkforums.broadinstitute.org/gatk/discussion/7696/pedigree-ped-files)



```{r Dada2}

library(dplyr)

load("~/DADA2.RData")
dir.create("~/Poppante/")
covariate_names <- c('BMI') # Insert the column names of your covariates as a vector
predictors_names <- c('BMI') # Insert the column names of your predictors as a vector

# Calculate relative abundances from the count table using the phyloseq object
physeq_prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))

# Extract taxa and abundance table
taxa_table <- as.data.frame(tax_table(physeq_prop))
abundance_table <- as.data.frame(otu_table(physeq_prop))
metadata <- as.data.frame(sample_data(physeq_prop))

#To make visualization easier, we can assign a unique ID to each of our sequence variants. I will create a new category in the taxonomy table called Seq that contains the SV sequence and also a category for the unique ID.
taxa_table$esv<-rownames(taxa_table)
taxa_table$esv_id<-paste0("esv_", seq(from=1, to=nrow(taxa_table), by=1))
colnames(abundance_table)<-taxa_table[colnames(abundance_table),]$esv_id #Our table of SVs currently has sequences for names, switch it to the unique ID you created above
rownames(taxa_table)<-taxa_table$esv_id
  
# We used the rank-based inverse normal transformation to counteract departures from normality. However, different procedures can be used such as the popular centered log-ratio transformation.
for(i in c(1:ncol(abundance_table))){
  abundance_table[,i]<- RNOmni::rankNormal(abundance_table[,i])}


# Add the metadata columns to the abundance table
abundance_table$Sample_name <- row.names(abundance_table)
metadata$Sample_name <- row.names(metadata)

abundance_table <- inner_join(metadata[,c(which(colnames(metadata) == 'Sample_name'),
                                            which(colnames(metadata) == 'Zygosity'),
                                            which(colnames(metadata) == 'Gender'),
                                            which(colnames(metadata) == 'Family_id'),
                                            which(colnames(metadata) == 'Individual_id'))], abundance_table, by = 'Sample_name')

# Create an empty PED file
Status = 0 ; FGender = 'F' ; MGender = 'M' ; FMID = 0 ; StatusFM = -9
PED_empty <- cbind(paste(abundance_table$Family_id), paste(abundance_table$Individual_id), paste(abundance_table$Family_id,3, sep = ""), paste(abundance_table$Family_id,4, sep = ""), paste(abundance_table$Gender), Status, as.character(abundance_table$Zygosity))
Father <- cbind(abundance_table$Family_id, paste(abundance_table$Family_id,3, sep = ""), FMID, FMID, FGender, StatusFM, FMID)
Mother <- cbind(abundance_table$Family_id, paste(abundance_table$Family_id,4, sep = ""), FMID, FMID, MGender, StatusFM, FMID)
rm(Status, MGender, StatusFM, FMID, FGender)
PED_empty <- rbind(PED_empty, Father, Mother) ; rm(Father, Mother)
colnames(PED_empty)[1] <- "Family_id" ; colnames(PED_empty)[2] <- "Individual_id" ; colnames(PED_empty)[7] <- "zygocity"
PED_empty <- as.data.frame(PED_empty, stringsAsFactors = FALSE)


# Fill the PED file with the taxa data used as a response
PED_filled <- left_join(PED_empty, abundance_table[-c(1,2,3,4)], by="Individual_id")
PED_filled <- distinct(PED_filled)
PED_filled$V5[PED_filled$V5 == 'F'] <- 2 ; PED_filled$V5[PED_filled$V5 == 'M'] <- 1
PED_filled$zygocity[PED_filled$zygocity == 'MZ'] <- 1 ; PED_filled$zygocity[PED_filled$zygocity == 'DZ'] <- 2
write.table(PED_filled, '~/Poppante/mydata.ped',  sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Create the covariate file
covariate_columns <- vector()
for (i in 1:length(covariate_names))
{
  covariate_columns <- append(covariate_columns, which(colnames(metadata) == covariate_names[i]))
}
covariate <- PED_filled[c(1,2)]
covariate <- inner_join(covariate, metadata[,c(which(colnames(metadata) == 'Individual_id'), covariate_columns)], by = 'Individual_id')
write.table (covariate, '~/Poppante/covariate.txt', sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE) 

# Create the response file
response <- colnames(PED_filled)[8:length(PED_filled)]
write.table(response, "~/Poppante/response.txt", sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Create the predictor file
predictors_columns <- vector()
for (i in 1:length(covariate_names))
{
  predictors_columns <- append(predictors_columns, which(colnames(metadata) == predictors_names[i]))
}
predictors <- as.data.frame(cbind(PED_filled$Family_id, PED_filled$Individual_id), stringsAsFactors = FALSE)
predictors <- left_join(predictors, metadata[,c(which(colnames(metadata) == 'Individual_id'), predictors_columns)], by="Individual_id")
write.table(predictors, '~/Poppante/predictors.txt', sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Create the map file
map <- predictors_names
write.table(map, '~/Poppante/map.txt', col.names = FALSE, row.names = FALSE, quote = FALSE)

```


### Make the PopPAnTe analysis

Make sure that the PopPAnTe executable is in the folder containing the tables created using the previous code chunk.


```{bash}
wget http://www.twinsuk.ac.uk/wp-content/uploads/2018/08/poppante_20180320-v1_0_2.tar.gz

tar -xvzf poppante_20180320-v1_0_2.tar.gz

java -jar poppante_20180320-v1_0_2.jar -mode association -ped mydata.ped -response response.txt -predictor predictors.txt -map map.txt -covariate covariate.tsv -output results -verbose true

```

Author:     Robin Mesnage
contact:    robin.mesnage@kcl.ac.uk
--------------------------------------------------------------------------

