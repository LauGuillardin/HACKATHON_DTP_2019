# HACKATHON_DTP_2019
## BIOINFORMATICS WORKFLOW
![](https://media.giphy.com/media/OkaMjJ7r7cMhO/giphy.gif)

### DOWNLOAD SRR DATA files from NCBI
 

### CONVERT SRR INTO FASTQ and #SPLIT PE FILES INTO _1 _2 (FORWARD AND REVERSE)


`[enter your path]/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --accession SRR5218239.1 --split-files --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+' --split-spot --gzip`


### QUALITY CONTROL USING FASTQC

`[enter your path]/FastQC/fastqc SRR5218242.1_1.fastq.gz` 

### TRIMMOMATIC

`java -jar [enter your path]/Trimmomatic-0.39/trimmomat-0.39.jar PE SRR5218239.1_1.fastq.gz SRR5218239.1_2.fastq.gz SRR5218239.1_1_TRIMMED_paired.fastq.gz SRR5218239.1_1_TRIMMED_unpaired.fastq.gz SRR5218239.1_2_TRIMMED_paired.fastq.gz SRR5218239.1_2_TRIMMED_unpaired.fastq.gz ILLUMINACLIP:[enter your path]/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

### QUALITY CONTROL USING FASTQC

`[enter your path]/FastQC/fastqc SRR5218242.1_1.fastq.gz` 

### KALLISTO
#### index
`[enter your path]/KALLISTO/kallisto index -i transcripts.idx MTYJ01.1.fsa_nt.gz`
#### quant active
`[enter your path]/KALLISTO/kallisto quant -i transcripts.idx -o output39 -b 100 <(zcat SRR5218239.1_1_TRIMMED_paired.fastq.gz) <(zcat SRR5218239.1_2_TRIMMED_paired.fastq.gz)`
#### quant tun
`[enter your path]/KALLISTO/kallisto quant -i transcripts.idx -o output42 -b 100 <(zcat SRR5218242.1_1_TRIMMED_paired.fastq.gz) <(zcat SRR5218242.1_2_TRIMMED_paired.fastq.gz)`

#### to extract just the counts for differential expression analysis
`cut -f 1,4 [enter your path]/output39/abundance.tsv | sed '1d' > active39.kallisto.counts.txt`

`cut -f 1,4 [enter your path]/output42/abundance.tsv | sed '1d' > tun42.kallisto.counts.txt`

#### sort
`sort -k 1 [enter your path]/active39.kallisto.counts.txt > active39.kallisto.counts.sorted.txt` 

`sort -k 1 [enter your path]/active42.kallisto.counts.txt > active42.kallisto.counts.sorted.txt` 

### DESEQ2 (R)

```R
setwd("[enter your path]")

library(DESeq2)

counts<-read.table("jointcount.sorted.csv", header = F)

counts

colnames(counts)<-c("active39", "active39", "tun42", "tun42")

counts_n<-counts[,c(2,4)]

counts_n[1:2] <- lapply(counts_n[1:2], round)

plot(counts_n)

keep <-rowSums(counts_n)>=10; counts_n<- counts_n[keep,]

coldata<-read.table("data_table.txt",row.names = 1, header = T)

dds<-DESeqDataSetFromMatrix(countData = counts_n, colData = coldata, design = ~condition)

dds$condition <- relevel(dds$condition, ref = "active")

dds<-DESeq(dds)

results_DEA<-results(dds, contrast = c("condition","active", "tun"))

plotMA(results_DEA)

idx<-identify(results_DEA$baseMean, results_DEA$log2FoldChange)

write.table(results_DEA, file = "results_DEA.tsv", sep = "\t", quote = F)


#code for creating a heatmap 
ntd <- normTransform(dds)

BiocManager::install("vsn")
library("vsn")

meanSdPlot(assay(ntd))

select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
