library(Rsubread)

setwd("/net/data.isilon/ag-saez/bq_rramirez/LifeInformatics/E-GEOD-55296")

# Target file
RawTargets = read.table("E-GEOD-55296.sdrf.txt",sep = "\t",header = T,stringsAsFactors = F)

#Download files automatically
FileName = RawTargets$Scan.Name

readfile1 = sort(FileName[grep("_1",FileName)])
readfile2 = sort(FileName[grep("_2",FileName)])

in_ix = grep("SRR1175569",readfile1)
las_ix = length(readfile1)

#Align Counts
BAMFiles = paste(unlist(lapply(strsplit(readfile1,split = "_1"),function(x){
  x[1]})),".subread.bam",sep = "")

#buildindex(basename="hg19_subread",reference="GRCh38.primary_assembly.genome.fa",colorspace=TRUE)

align("hg19_subread", readfile1 = readfile1[in_ix:las_ix], readfile2= readfile2[in_ix:las_ix], input_format = "gzFASTQ", output_file = BAMFiles[in_ix:las_ix])
gene <- featureCounts(BAMFiles, useMetaFeatures=TRUE, annot.inbuilt="hg38", allowMultiOverlap=TRUE)

save(gene, file="genecounts.ro")
