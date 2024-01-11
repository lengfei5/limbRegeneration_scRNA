##########################################################################
##########################################################################
# Project:
# Script purpose: clean GTF file of axolotl
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Feb 24 15:30:59 2021
##########################################################################
##########################################################################
##########################################
# functions used to clean the gtf 
##########################################
rm(list=ls())

extract.gene.transcript.gtf.attribute = function(test)
{
  save = rep(NA, 6)
  test = unlist(strsplit(as.character(test), ';'))
  jj = which(test == ' ')
  if(length(jj) >0 ) test = test[-jj]
  ii.gene = grep('gene_id', test)
  ii.transcript = grep('transcript_id', test)
  ii.other = setdiff(c(1:length(test)), unique(c(ii.gene, ii.transcript)))
  
  gene.id = test[ii.gene]
  gene.id = unlist(strsplit(as.character(gene.id), '[|]'))
  save[1] = gsub('gene_id ', '', paste0(gene.id[grep('AMEX6', gene.id, invert = TRUE)], collapse = ';'))
  
  transcript.id = test[ii.transcript]
  transcript.id = unlist(strsplit(as.character(transcript.id), '[|]'))
  save[2] = gsub('transcript_id ', '', paste0(gene.id[grep('AMEX6', transcript.id, invert = TRUE)], collapse = ';'))
  
  save[3] = paste0('gene_id "', gene.id[grep('AMEX6', gene.id)], '"; transcript_id "', transcript.id[grep('AMEX6', transcript.id)], '";')
  
  save[4] = gene.id[grep('AMEX6', gene.id)]
  save[5] = transcript.id[grep('AMEX6', transcript.id)]
  
  if(length(ii.other)>0)  save[6] = paste0(test[ii.other], collapse = ';')
  
  return(save)
  
}

extract.gene.id.gene.names = function(test)
{
  # test = keep$V9[1]
  save = rep(NA, 4)
  names(save) = c('gene.id', 'gene.symbol.nr', 'gene.symbol.hs', 'others')
  
  test = unlist(strsplit(as.character(test), ';'))
  jj = which(test == ' ')
  if(length(jj) >0 ) test = test[-jj]
  
  ii.id = grep('gene_id', test)
  ii.name = grep('gene_name', test)
  ii.other = setdiff(c(1:length(test)), unique(c(ii.id, ii.name)))
  if(length(ii.other)>0) save[4] = paste0(test[ii.other], collapse = ';')
  
  gene.id = test[ii.id]
  gene.id = gsub('gene_id ', '', gene.id)
  save[1] = gene.id
  gene.names = test[ii.name]
  gene.names = gsub('gene_name ', '', gene.names)
  gene.names = gsub(' ', '', gene.names)
  gene.names = unlist(strsplit(as.character(gene.names), '[|]'))
  name.nr = gene.names[grep('nr', gene.names)]
  name.hs = gene.names[grep('hs', gene.names)]
  if(length(name.nr)>0) save[2] = gsub('nr','', gsub("\\[|\\]", "", name.nr))
  if(length(name.hs)>0) save[3] = gsub('hs','', gsub("\\[|\\]", "", name.hs))
  return(save)
  
}

##########################################
#  import gtf data
##########################################

annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
aa = read.delim(file = paste0(annotDir, 'AmexT_v47.release.gtf'), sep = '\t', header = FALSE)

keep = aa

# remove contigs 
keep = keep[grep('chr', keep$V1), ]

write.table(keep, file = paste0(annotDir, 'AmexT_v47.release_rm.contigs.gtf'), sep = '\t', quote = FALSE, col.names = FALSE, 
            row.names = FALSE)

# select genes and make genn id to gene symbol mapping
keep = keep[which(keep$V3 == 'gene'), ]

mapping = t(sapply(keep$V9, extract.gene.id.gene.names))
mapping = data.frame(keep[, c(1:8)],  mapping, stringsAsFactors = FALSE)


write.table(mapping, file = paste0(annotDir, 'AmexT_v47.release_rm.contigs_geneID_geneSymbol.mapping.txt'), sep = '\t', 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


library("tictoc")
tic()

res = t(sapply(keep$V9, extract.gene.transcript.gtf.attribute))

colnames(res) = c('gene.names', 'transcript.names', 'attributes.gtf', 'gene.id.Amex60DD', 'transcript.id.Amex60DD', 'other.infos')  

toc()


keep$V9 = res[,3]
keep = keep[grep('chr', keep$V1), ]
write.table(keep, file = paste0(annotDir, 'ax6_UCSC_2021_01_26_attributes.cleaned_rm.contigs.gtf'), sep = '\t', quote = FALSE, col.names = FALSE, 
            row.names = FALSE)


res = data.frame(keep[, c(1:7)], res, stringsAsFactors = FALSE)

xx = res[grep('chr', res$V1),  ]



Test.for.loop = FALSE
if(Test.for.loop){
  
  library("tictoc")
  tic()
  for(n in 1:9000){
    #n = 1
    if(n%%1000 == 0) cat(n, '\n')
    
    test = keep$V9[n]
    test = unlist(strsplit(as.character(test), ';'))
    jj = which(test == ' ')
    if(length(jj) >0 ) test = test[-jj]
    ii.gene = grep('gene_id', test)
    ii.transcript = grep('transcript_id', test)
    ii.other = setdiff(c(1:length(test)), unique(c(ii.gene, ii.transcript)))
    
    gene.id = test[ii.gene]
    gene.id = unlist(strsplit(as.character(gene.id), '[|]'))
    keep$gene[n] = gsub('gene_id ', '', paste0(gene.id[grep('AMEX6', gene.id, invert = TRUE)], collapse = ';'))
    
    transcript.id = test[ii.transcript]
    transcript.id = unlist(strsplit(as.character(transcript.id), '[|]'))
    keep$transcript[n] = gsub('transcript_id ', '', paste0(gene.id[grep('AMEX6', transcript.id, invert = TRUE)], collapse = ';'))
    
    keep$new.attr[n] = paste0('gene_id "', gene.id[grep('AMEX6', gene.id)], '"; transcript_id "', transcript.id[grep('AMEX6', transcript.id)], '";')
    
    if(length(ii.other)>0)  keep$others[n] = paste0(test[ii.other], collapse = ';')
    
  }
  
  toc()
  
}
 



