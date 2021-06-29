#------------------------------------------------------#
#               GENERAL INFORMATION
#------------------------------------------------------#
#
# Script name:
#  Isoform-expression-analysis.R
#
# Purpose of script:
#  This code was made to perform transcript differential expression analysis on aligned RNA-seq data.
#
# Author:
#  Floris Vleugels i6190268
#  f.vleugels@student.maastrichtuniversity.nl
#  florisvleugels@hotmail.com
#
# Date created:
#  28-06-2021
#
# Personal comment:
#   This script requires some manual inputs and mainly acts as a template to aid the RNA-seq isoform expression analysis workflow that was created. 
#
#
# Inputs: 
#  BAM (binary alignment) files
#  Reference genome 
#
# Outputs: 
#  Counts_raw.csv
#  Counts_filtered.csv
#  Counts_junctions.csv
#  LikelihoodTest.csv or ExactTest.csv
#  summaryTest.csv
#  topTags.csv
#  BCVanalysis.PNG
#  MDplots.PNG
#  MDSplot.PNG
#  Volcanoplot.PNG
# 
#---------------------------------------------------------#
#                        INFO
#---------------------------------------------------------#
# The RNA-seq data that was used was downloaded from the SRA. The SRA toolkit was ran in windows CMD and the working directory was made through CMD. 
# The BAM file reference genomes can be found in the BAM file meta data, which can be analyzed in CMD. 
# Correct reference genome link can be selected on: https://hgdownload.soe.ucsc.edu/downloads.html#human
# Further analysis of the BAM files is done in R. 
# SRA toolkit can be downloaded and installed using the SRA github page: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
#
#---------------------------------------------------------#
#                DOWNLOAD SRA DATA IN CMD 
#---------------------------------------------------------#
# Open CMD:
#  press Windows+R, type "CMD" and press "OK"
#
# Create BAM directory:
#  cd \Users\flori\OneDrive\Bureaublad
#  mkdir Results
#
# Open SRA toolkit:
#  cd C:\Users\flori\OneDrive\Bureaublad\sratoolkit.2.11.0-win64\bin
#
# Download SRA data and specify genomic region:
# For example only chromosome 1 of data ERR3256992 is downloaded
#  sam-dump ERR3256992 --aligned-region chr1 --output-file C:\Users\flori\OneDrive\Bureaublad\Results\ERR3256992.bam
#
# Find reference genome build in the BAM file meta data:
#  findstr /n ID ERR3256992
#
#---------------------------------------------------------#
#              (INSTALL &) LOAD PACKAGES
#---------------------------------------------------------#
# list of required packages
s_requiredpackages =
  c(
    "edgeR",
    "ggplot2",
    "ggrepel",
    "Rsubread"
  )
#
#
# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (i in s_requiredpackages) {
  if (!requireNamespace(i, quietly = TRUE))
    BiocManager::install(i, ask = F)  # dependencies = c("Depends", "Imports")
  require(as.character(i), character.only = TRUE)
  #print(i)
}
#
#---------------------------------------------------------#
#                 GET REFERENCE GENOME
#---------------------------------------------------------#
# Set working directory to BAM directory:
setwd("C:/Users/flori/OneDrive/Bureaublad/Results")
#
# Download reference genome build:
download.file(hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
#
#---------------------------------------------------------#
#             Rsubread transcript counting
#---------------------------------------------------------#
# List all BAM files in the BAM directory:
bam <- list.files(pattern = ".bam")
#
# Store featureCounts in object for every BAM file in bam list:
fc = featureCounts(files = bam, annot.ext = "hg38.ncbiRefSeq.gtf",
                              isPairedEnd = F, 
                              useMetaFeatures = TRUE, 
                              GTF.attrType = "transcript_id", 
                              allowMultiOverlap = TRUE, largestOverlap = TRUE, 
                              isGTFAnnotationFile = TRUE, juncCounts = TRUE)
#
# Store raw transcript count data in object:
Counts_raw <- fc$counts
#
# Store junction counts in object: 
Counts_junction <- fc$counts_junction
#
# Export raw count data: 
write.table(Counts_raw, file = "Counts_raw", col.names = NA)
#
# Export junction count data:
write.table(jc1, file = "Counts_junction", col.names = NA)
#
#
#----------------------------------------------------------#
#                 edgeR count processing
#----------------------------------------------------------#
#
# Define study design (groups):
group <- c(1,1,1,2,2,2)
#
# Make DGE object:
DGE <- DGEList(counts = Counts_raw, group = group)
#
# Filter out the low count transcripts:
keep <- filterByExpr(DGE)
DGE <- DGE[keep, , keep.lib.sizes =FALSE]
#
# Show total inclusion/exclusion:
summary(keep)
#
# Normalize the count data:
DGE <- calcNormFactors(DGE)
#
# Estimate dispersions:
#  Replicates:
DGE <- estimateDisp(DGE)
#  No replicate samples:
bcv <- 0.4
#
# Make a BCV plot to visualize dispersions:
plotBCV(DGE)
export
#
# Quality control of normalization through mean-difference plot for all BAM files:
plotMD.DGEList(DGE, column = 1)
abline(h=0, col="red", lty=2, lwd=2)
export
#
# Export processed count table:
write.table(DGE$counts, file = "Counts_filtered.csv", col.names = NA )
#
# Examine samples for outliers and relationships:
plotMDS.DGEList(DGE)
export
#
#-----------------------------------------------------------#
#               Differential expression analysis
#-----------------------------------------------------------#
#
# For replicates perform likelihood ratio test, exact test for no replicates.
#
# Perform exact test or likelihood ratio test:
eT <- exactTest(DGE, dispersion = bcv^2)
fit <- glmFit(DGE)
lrt <- glmLRT(fit)
#
# Top results of exact test or likelihood ratio test: 
tT <- topTags(eT)
tT <- topTags(lrt)
#
# Summarize total differential expression of transcripts:
summary(decideTestsDGE(lrt))
summary(decideTestsDGE(et))
write.table(summary(decideTestsDGE(lrt)), file = "summary.csv", sep = ",", col.names = NA)
#
# Export exact test or likelihood ratio test and top tags:
write.table(lrt, file = "likelihood.csv", sep = ",", col.names = NA)
write.table(eT, file = "ExactTest.csv", sep = ",", col.names = NA)
write.table(tT, file = "topTags.csv", sep = ",", col.names = NA)
#
#----------------------------------------------------------#
#                 Mean difference plots
#----------------------------------------------------------#
#
# Mean-difference plots for likelihood test or exact test:
plotMD.DGELRT(lrt)
plotMD.DGEExact(et)
abline(h=c(-1, 1), col="blue")
export
#
#----------------------------------------------------------#
#                       Volcano plot
#----------------------------------------------------------#
#
# Create data frame for volcano plot, either from likelihood or exact test:
volcdata <- et$table
volcdata <- lrt$table
#
# Add column to data frame for expression status:
volcdata$diffexpressed <- "NO"
volcdata$diffexpressed[volcdata$logFC > 1.0 & volcdata$PValue < 0.05] <- "UP"
volcdata$diffexpressed[volcdata$logFC < -1.0 & volcdata$PValue < 0.05] <- "DOWN"
#
# Define colours for Up Down and No diffexepressed:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
#
# Show transcript names based on log p value:
volcdata$logpv<- -log10(Pvalue)
volcdata$label[test$logpv>43]<- rownames(volcdata)[volcdata$logpv>43]
#
# Construct volcano plot:
ggplot(data=volcdata, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=label))+ 
  geom_point() + 
  theme_minimal() + 
  scale_colour_manual(values = mycolors) +
  geom_text_repel()
#
# Export Volcano plot:
export
#
#
#
#
  
  
  
  
  


























