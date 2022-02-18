### ---------------
### Jian-Guo Zhou MD.
###### Date: "Thu Sep 23 10:25:33 2021"
### Email: jianguo.zhou@fau.de

###
### ---------------
rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir=getwd()
Figure_dir=getwd()

# load you data
#1.GSE33774
meta=samples1_n[!(samples1_n$group=="periodontitis"),]
expr=t(datExpr0[rownames(datExpr0) %in% rownames(meta),])
expr=data.frame(expr)
expr=round(expr)
#2.GSE106090
meta=test.clin1[!(test.clin1$Group=="Periodontitis"),]
expr=t(test_n[rownames(test_n) %in% rownames(meta),])
expr=data.frame(expr)
expr=round(expr)
dim(expr)
# [1] 60591    80
dim(meta)
# [1] 80 71


#  comparation: implantitis vs healthy
rownames(meta)
#1.GSE33774
group_list=ifelse(meta$Group =='implantitis','implantitis','healthy')
#2.GSE106090
group_list=ifelse(meta$Group =='Periimplantitis','implantitis','healthy')
table(group_list)
exprSet=na.omit(expr)
# rm(DEG_limma_voom)
# rm(DESeq2_DEG)
# rm(edgeR_DEG)
#source('./functions.R')
source('./functions.R')
### ---------------
###
### Firstly run DESeq2 
###
### ---------------

if(T){
  library(DESeq2)
  
  (colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  tmp_f=file.path(Rdata_dir,'GSE33774-DESeq2-dds.Rdata')
  tmp_f=file.path(Rdata_dir,'GSE106090-DESeq2-dds.Rdata')
  if(!file.exists(tmp_f)){
    dds <- DESeq(dds)
    save(dds,file = tmp_f)
  }
  load(file = tmp_f)
  res <- results(dds, 
                 contrast=c("group_list","implantitis","healthy"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  
  nrDEG=DESeq2_DEG[,c(2,6)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')  
  draw_h_v(exprSet,nrDEG,'DEseq2',group_list,1)
}

### ---------------
###
### Then run edgeR 
###
### ---------------
if(T){
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG =nrDEG 
  nrDEG=edgeR_DEG[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'edgeR',group_list,1)
  
}


### ---------------
###
### Lastly run voom from limma
###
### --------------- 
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('implantitis-healthy'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='implantitis-healthy', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'limma',group_list,1)
  
}

tmp_f=file.path(Rdata_dir,'GSE33774-DEG_results.Rdata')
tmp_f=file.path(Rdata_dir,'GSE106090-DEG_results.Rdata')

if(file.exists(tmp_f)){
  save(DEG_limma_voom,edgeR_DEG, file = tmp_f)
  
}else{
  load(file = tmp_f) 
}



nrDEG1=DEG_limma_voom[,c(1,4)]
colnames(nrDEG1)=c('log2FoldChange','pvalue') 

nrDEG2=edgeR_DEG[,c(1,5)]
colnames(nrDEG2)=c('log2FoldChange','pvalue') 

nrDEG3=DESeq2_DEG[,c(2,6)]
colnames(nrDEG3)=c('log2FoldChange','pvalue')  

mi=unique(c(rownames(nrDEG1),rownames(nrDEG1),rownames(nrDEG1)))
lf=data.frame(lf1=nrDEG1[mi,1],
              lf2=nrDEG2[mi,1],
              lf3=nrDEG3[mi,1])
cor(na.omit(lf))
#         lf1       lf2       lf3
# lf1 1.0000000 0.5170625 0.5236079
# lf2 0.5170625 1.0000000 0.9395551
# lf3 0.5236079 0.9395551 1.0000000
save(nrDEG1,nrDEG2,nrDEG3,lf,file = "lf.DEG.Rdata")
save(DEG_limma_voom,edgeR_DEG, file = tmp_f)
