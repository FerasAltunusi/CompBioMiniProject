#Download different packages "sleuth", and "dblyr" 
library(sleuth)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if(args[1]=="True"){
  output_dir<-"test_outputs"
} else {
  output_dir<-"miniProject_Feras_Altunusi"
}

#read in the table describing samples and kallisto output
stab <- read.table(paste(output_dir,"/kallisto/sample_table.tsv",sep=""),header=TRUE,stringsAsFactors = FALSE,sep="\t")

#initialize sleuth
so <- sleuth_prep(stab)

#compare the two conditions
so <- sleuth_fit(so, ~days.post.infection, 'full')

#compare the likelihood ratio test by using "reduced"
so <- sleuth_fit(so, ~1, 'reduced')

#get the differential expression between both conditions (2dpi & 6dpi)
so <- sleuth_lrt(so,'reduced','full')

#extract the test results from the sleuth object 
sleuth_table <- sleuth_results(so,'reduced:full','lrt',show_all = FALSE)

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

#write the info for the significant SNPs to the output file
output <- sleuth_significant[,c(1,4,2,3)]
write.table(output, file=paste(output_dir,"/miniProject.log",sep=""),quote=FALSE,row.names=TRUE,sep="\t",append=TRUE)
