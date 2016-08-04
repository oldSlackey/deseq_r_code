#Scripts to run DESEQ2 analyses on Mariolas data. Again, the count data are taken from the 3i
#pipeline, so they have been produced using the SubRead SRA (short read aligner...)
#Start/Install Biocinductor
source("http://bioconductor.org/biocLite.R")
#Get help on upgrading bioC
?BiocUpgrade
#Carry out BioConductor Upgrade
biocLite("BiocUpgrade")
#Load DESEQ2 library
biocLite("DESeq2")
biocLite("ggplot2")
biocLite("vsn")
biocLite("RColorBrewer")
biocLite("gplots")
biocLite("biomaRt")
biocLite("dplyr")
biocLite("ggrepel")
biocLite("S4Vectors")
biocLite("GenomicRanges")
library("DESeq2")
library("ggplot2")
library("vsn")
library("RColorBrewer")
library("gplots")
library("biomaRt")
library("dplyr")
library("ggrepel")
library("S4Vectors")
library("GenomicRanges")
library("Matrix")
#Load gene counts as data matrix
gene_counts <- read.csv("list.csv.gene.count.csv.tmp", row.names=1)
gene_counts <- as.matrix(gene_counts)
#Correct sample names in count matrix
samples<-strsplit(colnames(gene_counts),split="_")
unlisted_samples<-unlist(samples)
unlisted_indices<-grep("kal|nbval",unlisted_samples, invert=TRUE)
colnames(gene_counts)<-unlisted_samples[unlisted_indices]

#Then extract the sample names to produce the colData object for DESeq2
colData<-colnames(gene_counts)
colData<-as.data.frame(colData)
#Read in sample description file
exp_design<-read.table("column_data.txt", sep="\t", row.names=1, header=TRUE)
exp_design<- as.data.frame(exp_design)
exp_design$names<-rownames(exp_design)
exp_design2 <- exp_design[order(rownames(exp_design)),]
exp_design2<-as.data.frame(exp_design2)
exp_design<-exp_design2
rm(exp_design2)

#need to add in columns to specify treatment groups
mysplitPeriod<-function(x) {
  list<-strsplit(x, split="\\.")
  return (list[[1]][2])
}
mysplitUnderscore<-function(x) {
  list<-strsplit(x, split="_")
  return (list[[1]][2])
}
exp_design$control<-lapply(exp_design$names,mysplitPeriod)
exp_design$SAMPLE<-as.character(exp_design$SAMPLE)
exp_design$gene<-lapply(exp_design$SAMPLE,mysplitUnderscore)
exp_design$cont_gene<-paste(exp_design$control,exp_design$gene,sep="_")
colData<-exp_design
colData$cont_gene<-sub("SAM_2","SAM_A",colData$cont_gene)
#Ignore all of this, what a palava...
#colData$control<-exp_design$control
#colData$SAMPLE<-exp_design$SAMPLE
#colData$gene<-exp_design$gene
#colData$cont_gene<-exp_design$cont_gene
#row.names(colData)<-colData$colData
#colData$colData<-NULL
# order the samples in countTable the same as in the colTable
sample.order<-as.vector(row.names(colData))
orderedCountTable<-gene_counts[,sample.order]
#Copying Rickys, Pawels and Grahams methods
#Create vector of contrasts although in this case this is useless...
contrasts<-c("SAMAvsCONTA","SAMBvsCONTB","SAMAvsSAMB")
#Subset data and create new objects according to contrasts required
colDataSAMAvsCONTA<-colData[ which(colData$cont_gene=="CON_A" | colData$cont_gene=="SAM_A"),]
SAMAvsCONTA_samples<-row.names(colDataSAMAvsCONTA)
orderedCountTableSAMAvsCONTA<-orderedCountTable[,SAMAvsCONTA_samples]
ddsSAMAvsCONTA <- DESeqDataSetFromMatrix(orderedCountTableSAMAvsCONTA,colData=colDataSAMAvsCONTA,design=~cont_gene)
ddsSAMAvsCONTA <- DESeq(ddsSAMAvsCONTA)

colDataSAMBvsCONTB<-colData[ which(colData$cont_gene=="CON_B" | colData$cont_gene=="SAM_B"),]
SAMBvsCONTB_samples<-row.names(colDataSAMBvsCONTB)
orderedCountTableSAMBvsCONTB<-orderedCountTable[,SAMBvsCONTB_samples]
ddsSAMBvsCONTB <- DESeqDataSetFromMatrix(orderedCountTableSAMBvsCONTB,colData=colDataSAMBvsCONTB,design=~cont_gene)
ddsSAMBvsCONTB <- DESeq(ddsSAMBvsCONTB)

colDataSAMAvsSAMB<-colData[ which(colData$cont_gene=="SAM_A" | colData$cont_gene=="SAM_B"),]
SAMAvsSAMB_samples<-row.names(colDataSAMAvsSAMB)
orderedCountTableSAMAvsSAMB<-orderedCountTable[,SAMAvsSAMB_samples]
ddsSAMAvsSAMB <- DESeqDataSetFromMatrix(orderedCountTableSAMAvsSAMB,colData=colDataSAMAvsSAMB,design=~cont_gene)
ddsSAMAvsSAMB <- DESeq(ddsSAMAvsSAMB)

#Running analysis again with samples N8.SAM and NS17.SAM removed
test <- orderedCountTable
test[1:10,11]
test[1:10,"N8.SAM"]
test<-test[,-11]
test[1:10,16]
test[1:10,"NS17.SAM"]
test<-test[,-16]
colDataTestSAMAvsCONTA<-colData[ which(rownames(colData)!="N8.SAM" | rownames(colData)!="NS17.SAM")]
colDataN8Gone<-colData[-11,]
colDataN8Gone[16,]
colDataN8NS17Gone<-colDataN8Gone[-16,]

colDataExcSAMAvsCONTA<-colDataN8NS17Gone[ which(colDataN8NS17Gone$cont_gene=="CON_A" | colDataN8NS17Gone$cont_gene=="SAM_A"),]
SAMAvsCONTAExc_samples<-row.names(colDataExcSAMAvsCONTA)
orderedCountTableExcSAMAvsCONTA<-test[,SAMAvsCONTAExc_samples]
ddsExcSAMAvsCONTA <- DESeqDataSetFromMatrix(orderedCountTableExcSAMAvsCONTA,colData=colDataExcSAMAvsCONTA,design=~cont_gene)
ddsExcSAMAvsCONTA <- DESeq(ddsExcSAMAvsCONTA)

colDataExcSAMBvsCONTB<-colDataN8NS17Gone[ which(colDataN8NS17Gone$cont_gene=="CON_B"  | colDataN8NS17Gone$cont_gene=="SAM_B"),]
SAMBvsCONTBExc_samples<-row.names(colDataExcSAMBvsCONTB)
orderedCountTableExcSAMBvsCONTB<-test[,SAMBvsCONTBExc_samples]
ddsExcSAMBvsCONTB <- DESeqDataSetFromMatrix(orderedCountTableExcSAMBvsCONTB,colData=colDataExcSAMBvsCONTB,design=~cont_gene)
ddsExcSAMBvsCONTB <- DESeq(ddsExcSAMBvsCONTB)

colDataExcSAMAvsSAMB<-colDataN8NS17Gone[ which(colDataN8NS17Gone$cont_gene=="SAM_A" | colDataN8NS17Gone$cont_gene=="SAM_B"),]
ExcSAMAvsSAMB_samples<-row.names(colDataExcSAMAvsSAMB)
orderedCountTableExcSAMAvsSAMB<-test[,ExcSAMAvsSAMB_samples]
ddsExcSAMAvsSAMB <- DESeqDataSetFromMatrix(orderedCountTableExcSAMAvsSAMB,colData=colDataExcSAMAvsSAMB,design=~cont_gene)
ddsExcSAMAvsSAMB <- DESeq(ddsExcSAMAvsSAMB)

#Running analyses again, but distinguishing between different knockdown genes (A 1 and 2 and B one and 2
#vs their respective controls)
colDataGeneA1SAMvsCON<-colData[ which(colData$SAMPLE..GROUP=="shNTC_A" | colData$SAMPLE..GROUP=="shGENE_A_1"),]
SAMA1vsCONA_samples<-row.names(colDataGeneA1SAMvsCON)
orderedCountTableA1vsCONA<-orderedCountTable[,SAMA1vsCONA_samples]
ddsSAMA1vsCONA <- DESeqDataSetFromMatrix(orderedCountTableA1vsCONA,colData=colDataGeneA1SAMvsCON,design=~SAMPLE)
ddsSAMA1vsCONA <- DESeq(ddsSAMA1vsCONA)

colDataGeneA2SAMvsCON<-colData[ which(colData$SAMPLE..GROUP=="shNTC_A" | colData$SAMPLE..GROUP=="shGENE_A_2" | colData$SAMPLE..GROUP=="shGENE A_2"),]
SAMA2vsCONA_samples<-row.names(colDataGeneA2SAMvsCON)
orderedCountTableA2vsCONA<-orderedCountTable[,SAMA2vsCONA_samples]
ddsSAMA2vsCONA <- DESeqDataSetFromMatrix(orderedCountTableA2vsCONA,colData=colDataGeneA2SAMvsCON,design=~SAMPLE)
ddsSAMA2vsCONA <- DESeq(ddsSAMA2vsCONA)

colDataGeneB1SAMvsCON<-colData[ which(colData$SAMPLE..GROUP=="shNTC_B" | colData$SAMPLE..GROUP=="shGENE_B_1"),]
SAMB1vsCONB_samples<-row.names(colDataGeneB1SAMvsCON)
orderedCountTableB1vsCONB<-orderedCountTable[,SAMB1vsCONB_samples]
ddsSAMB1vsCONB <- DESeqDataSetFromMatrix(orderedCountTableB1vsCONB,colData=colDataGeneB1SAMvsCON,design=~SAMPLE)
ddsSAMB1vsCONB <- DESeq(ddsSAMB1vsCONB)

colDataGeneB2SAMvsCON<-colData[ which(colData$SAMPLE..GROUP=="shNTC_B" | colData$SAMPLE..GROUP=="shGENE_B_2"),]
SAMB2vsCONB_samples<-row.names(colDataGeneB2SAMvsCON)
orderedCountTableB2vsCONB<-orderedCountTable[,SAMB2vsCONB_samples]
ddsSAMB2vsCONB <- DESeqDataSetFromMatrix(orderedCountTableB2vsCONB,colData=colDataGeneB2SAMvsCON,design=~SAMPLE)
ddsSAMB2vsCONB <- DESeq(ddsSAMB2vsCONB)

dds_list<-c(ddsSAMA1vsCONA,ddsSAMA2vsCONA,ddsSAMB1vsCONB,ddsSAMB2vsCONB)
for (dds in dds_list){
  contrast_string<-unique((colData(dds)$SAMPLE))
  contr<-paste0(contrast_string[1],"vs",contrast_string[2])

  #Experimenting with different options
  #test_res <- results(dds)

  # Normalised counts
  disp<-estimateSizeFactors(dds)
  norm_counts<-counts(disp,normalized=TRUE)
  # Dispersion plot
  # plotDispEsts(dds)
  #pdf("DispPlot.pdf")
  #plotDispEsts(dds,main="Gene level")
  #dev.off()

#These are Ricky Cunninghams notes on deseq
## Obtaining the DESeq gene list
# train.cl is a vector that contains the class information for each sample
#temp <- as.data.frame(train.cl)
#rownames(temp) <- colnames(gene.data)# gene.data is the count-table that can be read into R using the read.table function
### DESeq requires integer input, so the count-table must be rounded (Kallisto output is non-integer estimated counts)
#temp.data <- apply(gene.data,2,function(x){round(as.numeric(x))})
#rownames(temp.data) <- rownames(gene.data)
#dds <- DESeqDataSetFromMatrix(temp.data, temp, design = ~ train.cl)
#dds <- DESeq(dds)
#res <- results(dds)  # This res variable contains all relevant results"
#Using Grahams script to add more annotations to Mariolas DE output
# BiomaRt database to use
  db<-"hsapiens_gene_ensembl"
  mart<-"ensembl"
  filt<-"ensembl_gene_id"

# Download transcript biotype annotations and add tot he DESeq2 object
# Create the biomaRt object
  ensembl = useEnsembl(biomart=mart, dataset=db)
# Use this to list available annotations
  martAttr <- listAttributes(ensembl)
# List of annotations to download
#att<-c("ensembl_gene_id","gene_biotype")
  att<-c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","gene_biotype")
# Get the biotypes from biomaRt
  biotypes<-getBM(attributes=att,filters=filt,values=rownames(dds),mart=ensembl)
# Get the annotations into a data frame
  DESeq2Features <- data.frame(ensembl_gene_id = rownames(dds))
  DESeq2Features$ensembl_gene_id <- as.character(DESeq2Features$ensembl_gene_id)
  rowData <- dplyr::left_join(DESeq2Features, biotypes, by = att[1])
  rowData <- as(rowData, "DataFrame")

# Add biotypes to DESQ2 object
  mcols(rowRanges(dds)) <- cbind(mcols(rowRanges(dds)),rowData)
#Specifying contrasts now
#Working with contrasts

  treat <- strsplit(contr,"vs")[[1]][1]
  ref <- strsplit(contr,"vs")[[1]][2]
  comp <- results(dds,contrast=c("SAMPLE",treat,ref))
  # Annotate results with biotypes
  # Add biotypes to results
  comp$gene_biotype <- mcols(dds)$gene_biotype
  # Annotate results with gene names
  comp$external_gene_name <- mcols(dds)$external_gene_name
  # # Annotate results with chromosome
  comp$chromosome_name <- mcols(dds)$chromosome_name
  # # Annotate results with start
  comp$start_position <- mcols(dds)$start_position
  # # Annotate results with end
  comp$end_position <- mcols(dds)$end_position
  #Order by adjusted p-values
  comp_sorted <- comp[order(comp$padj),]
  # Convert into matrix
  comp_sorted_matrix<-as.matrix(comp_sorted)
  # Merge with normalized counts
  comp_sorted_matrix_merged<-merge(comp_sorted_matrix,norm_counts, by = "row.names", all = TRUE)
  outfile <- paste(contr,'.Gene','ds2','csv',sep=".")
  write.csv(comp_sorted_matrix_merged, file=outfile)
  #Order by log-fold change
  comp_lfc1 <- results(dds,contrast=c("SAMPLE",treat,ref),lfcThreshold=1)
  comp_lfc1_sorted <- comp_lfc1[order(comp_lfc1$padj),]
  outfile <- paste(contr,'Gene','ds2','lfc1','csv',sep=".")
  write.csv(comp_lfc1_sorted, file=outfile)
  #Creating plots - ignore for now
  #rmafile <- paste(contr,'Gene','ds2','plotMA','pdf',sep=".")
  #pdf(rmafile)
  #main <- paste('Gene:',contr,sep=" ")
  #plotMA(comp,alpha=0.1,main=main,ylim=c(-2,2))
  #dev.off()
  
  #summary(comp)
  #print(paste('contrast: ',contr,' DONE'))
}
for (dds in dds_list){
  cond<-"cont_gene"
  rld <- rlog(ddsSAMAvsSAMB)
  i<-1
  data <- plotPCA(rld, intgroup = cond, returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  #plot_name = paste("rlog_PCA", i, ".pdf")
  #pdf(plot_name)
  ggplot(data, aes(PC1, PC2, color=group, label=name)) +  geom_point(size=5) + geom_text_repel(box.padding = unit(0.5, 'lines')) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
  #dev.off()
  i = i+1
}
plot(data)
