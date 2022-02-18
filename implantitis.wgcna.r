# GSE106090
# GSE33774
load("/Volumes/B/xjqqrza/implantitis/GSE106090.expr.Rdata")
GSE106090.expr<-expr
load("/Volumes/B/xjqqrza/implantitis/GSE33774.expr.Rdata")
GSE33774.expr<-expr
id1<-intersect(rownames(GSE106090.expr),rownames(GSE33774.expr))
load("/Volumes/B/xjqqrza/PRJEB23709_WGCNA/ICB_melanoma_WGCNA_1/gtf_v34_gene.Rdata")
id.exprs<-gtf_v34_gene[gtf_v34_gene$gene_name %in% id1,]
table(id.exprs$gene_type)
# IG_C_gene                          IG_V_gene                    IG_V_pseudogene 
# 1                                  2                                  1 
# lncRNA             polymorphic_pseudogene               processed_pseudogene 
# 365                                 19                                 27 
# protein_coding                           ribozyme                          TR_V_gene 
# 16202                                  1                                  2 
# transcribed_processed_pseudogene     transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene 
# 48                                 27                                136 
# unprocessed_pseudogene 
# 12 
id.exprs.lnc<-id.exprs[id.exprs$gene_type=="lncRNA",]
id.exprs.pc<-id.exprs[id.exprs$gene_type=="protein_coding",]
GSE106090.lnc<-GSE106090.expr[rownames(GSE106090.expr) %in% id.exprs.lnc$gene_name,]
GSE33774.lnc<-GSE33774.expr[rownames(GSE33774.expr) %in% id.exprs.lnc$gene_name,]
GSE106090.pc<-GSE106090.expr[rownames(GSE106090.expr) %in% id.exprs.pc$gene_name,]
GSE33774.pc<-GSE33774.expr[rownames(GSE33774.expr) %in% id.exprs.pc$gene_name,]
colnames(Targets_all)
colnames(Targets_all)<-c("GEO","FileName","sample_id","sample_name","tissue_type","group" )
GSE106090.clin<-Targets_all[Targets_all$GEO=="GSE106090",]
GSE106090.clin<-data.frame(GSE106090.clin)
rownames(GSE106090.clin)<-GSE106090.clin$sample_id
GSE33774.clin<-Targets_all[Targets_all$GEO=="GSE33774",]
GSE33774.clin<-data.frame(GSE33774.clin)
rownames(GSE33774.clin)<-GSE33774.clin$sample_id
GSE33774.clin<-eSet[["GSE33774_series_matrix.txt.gz"]]@phenoData@data
GSE106090.clin<-eSet[["GSE106090_series_matrix.txt.gz"]]@phenoData@data
GSE178351.clin<-eSet[["GSE178351_series_matrix.txt.gz"]]@phenoData@data
GSE57631.clin<-eSet[["GSE57631_series_matrix.txt.gz"]]@phenoData@data
###2021.09.13#####
library(readxl)
clin_data <- read_excel("clin_data.xlsx")
table(clin_data$Group)
# Healthy Periimplantitis   Periodontitis 
# 14              13              13 
clin_data$group<-ifelse(clin_data$Group=="Periimplantitis",1,
                        ifelse(clin_data$Group=="Periodontitis",2,0))
table(clin_data$group)
# 0  1  2 
# 14 13 13
clin_data$Periimplantitis<-ifelse(clin_data$Group=="Periimplantitis",1,0)
clin_data$Periodontitis<-ifelse(clin_data$Group=="Periodontitis",1,0)
clin_data$gender<-ifelse(clin_data$Gender=="male",1,0)
clin_data$age<-ifelse(clin_data$Age>60,1,0) 
clin_data<-data.frame(clin_data)
rownames(clin_data)<-clin_data$Sample_ID
colnames(clin_data)[1]<-"GEO"
table(clin_data[,1])
train.clin1<-clin_data[clin_data$GEO=="GSE33774",]
train.clin<-train.clin1[,19:22]
test.clin1<-clin_data[!(clin_data$GEO=="GSE33774"),]
test.clin<-test.clin1[,19:22]
train.expr<-GSE33774.pc[,colnames(GSE33774.pc) %in% rownames(train.clin1)]
test.expr<-GSE106090.pc[,colnames(GSE106090.pc) %in% rownames(test.clin1)]
save(train.expr,test.expr,train.clin1,test.clin1,file = "train.test.Rdata")

######TIL.28.signatures#######
library(IOBR)
library(EPIC)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(readxl)
TIL_28 <- read_excel("TISIDB/TIL_28.xlsx")
TIL_28.signature<-NULL
for (i in 1:NCOL(TIL_28)) {
  a<-TIL_28[,i]
  a<-a[,]
  a<-as.matrix(a)
  a<-na.omit(a)
  TIL_28.signature[[i]]<-a
}
names(TIL_28.signature)<-colnames(TIL_28)
save(TIL_28.signature,file = "TISIDB/TIL_28.signature.Rdata")
expr<-GSE33774.expr
expr<-GSE106090.expr
TIL_ssgsea<-calculate_sig_score(pdata   = NULL,
                                        eset            = expr,
                                        signature       = TIL_28.signature,
                                        method          = "ssgsea",
                                        adjust_eset = TRUE,
                                        mini_gene_count = 5)
GSE33774.TIL.28_ssgsea<-TIL_ssgsea
GSE33774.TIL.28_ssgsea<-data.frame(GSE33774.TIL.28_ssgsea)
rownames(GSE33774.TIL.28_ssgsea)<-GSE33774.TIL.28_ssgsea$ID
GSE106090.TIL.28_ssgsea<-TIL_ssgsea
GSE106090.TIL.28_ssgsea<-data.frame(GSE106090.TIL.28_ssgsea)
rownames(GSE106090.TIL.28_ssgsea)<-GSE106090.TIL.28_ssgsea$ID

####correlation  #######
TILs<-GSE33774.TIL.28_ssgsea[,3:30]
data.marker<-train.expr[rownames(train.expr) %in% markers,]

TILs<-GSE106090.TIL.28_ssgsea[,3:30]
data.marker<-test.expr[rownames(test.expr) %in% markers,]


p1<-NULL
for (i in 1:nrow(data.marker)) {
  id.1<-rownames(data.marker)[i]
  x<-t(data.marker[i,])
  for (b in 1:ncol(TILs)) {
    id.2<-colnames(TILs)[b]
    y<-TILs[,b]
    data<-data.frame(x,y)
    colnames(data)<-c(id.1,id.2)
    cor<-cor.test(x, y, method=c("pearson", "kendall", "spearman"))
    p.value<-cor[["p.value"]]
    pearson.p<-cor[["estimate"]][["cor"]]
    p<-data.frame(id.1,id.2,p.value,pearson.p)
    p1<-rbind(p1,p)
    print(i)
    print(b)
  }
}
GSE33774.p1<-p1
GSE106090.p1<-p1
write.csv(GSE33774.p1,file = "GSE33774.p1.csv")
write.csv(GSE106090.p1,file = "GSE106090.p1.csv")


####### plot GSEA ######

load("/Volumes/B/xjqqrza/implantitis/implantitis.brown.DEres2.GSEA.Rdata")
library(readxl)
GSEA_up <- read_excel("GSEA_up.xlsx")
GSEA_down <- read_excel("GSEA_down.xlsx")

GSE33774.DEG_limma_brown.gsea<-(test[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]]@result)
GSE106090.DEG_limma_brown.gsea<-test[["allpathway2gene"]][["GSE106090.DEG_limma_brown"]]@result
write.csv(GSE33774.DEG_limma_brown.gsea,file = "GSE33774.DEG_limma_brown.gsea.csv")
write.csv(GSE106090.DEG_limma_brown.gsea,file = "GSE106090.DEG_limma_brown.gsea.csv")
GSE33774.DEG_limma_brown.gsea.down<-GSE33774.DEG_limma_brown.gsea[GSE33774.DEG_limma_brown.gsea$enrichmentScore<0,]
GSE33774.DEG_limma_brown.gsea.down<-GSE33774.DEG_limma_brown.gsea.down[GSE33774.DEG_limma_brown.gsea.down$pvalue<0.05,]
#GSE33774.DEG_limma_brown.gsea.down<-GSE33774.DEG_limma_brown.gsea.down[GSE33774.DEG_limma_brown.gsea.down$p.adjust<0.05,]
GSE33774.DEG_limma_brown.gsea.down<-GSE33774.DEG_limma_brown.gsea.down[GSE33774.DEG_limma_brown.gsea.down$qvalues<0.2,]
GSE33774.DEG_limma_brown.gsea.up<-GSE33774.DEG_limma_brown.gsea[GSE33774.DEG_limma_brown.gsea$enrichmentScore>0,]
GSE33774.DEG_limma_brown.gsea.up<-GSE33774.DEG_limma_brown.gsea.up[GSE33774.DEG_limma_brown.gsea.up$pvalue<0.05,]
GSE33774.DEG_limma_brown.gsea.up<-GSE33774.DEG_limma_brown.gsea.up[GSE33774.DEG_limma_brown.gsea.up$p.adjust<0.05,]
GSE33774.DEG_limma_brown.gsea.up<-GSE33774.DEG_limma_brown.gsea.up[GSE33774.DEG_limma_brown.gsea.up$qvalues<0.05,]
GSE33774.DEG_limma_brown.gsea.up<-GSE33774.DEG_limma_brown.gsea.up[GSE33774.DEG_limma_brown.gsea.up$setSize>20,]

GSE106090.DEG_limma_brown.gsea.down<-GSE106090.DEG_limma_brown.gsea[GSE106090.DEG_limma_brown.gsea$enrichmentScore<0,]
GSE106090.DEG_limma_brown.gsea.down<-GSE106090.DEG_limma_brown.gsea.down[GSE106090.DEG_limma_brown.gsea.down$pvalue<0.05,]
# GSE106090.DEG_limma_brown.gsea.down<-GSE106090.DEG_limma_brown.gsea.down[GSE106090.DEG_limma_brown.gsea.down$p.adjust<0.05,]
GSE106090.DEG_limma_brown.gsea.down<-GSE106090.DEG_limma_brown.gsea.down[GSE106090.DEG_limma_brown.gsea.down$qvalues<0.2,]
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea[GSE106090.DEG_limma_brown.gsea$enrichmentScore>0,]
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea.up[GSE106090.DEG_limma_brown.gsea.up$pvalue<0.05,]
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea.up[GSE106090.DEG_limma_brown.gsea.up$p.adjust<0.05,]
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea.up[GSE106090.DEG_limma_brown.gsea.up$qvalues<0.05,]
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea.up[GSE106090.DEG_limma_brown.gsea.up$setSize>20,]

co.DESeq.gsea.up<-intersect(GSE106090.DEG_limma_brown.gsea.up$ID,GSE33774.DEG_limma_brown.gsea.up$ID)
co.DESeq.gsea.down<-intersect(GSE106090.DEG_limma_brown.gsea.down$ID,GSE33774.DEG_limma_brown.gsea.down$ID)
GSE106090.DEG_limma_brown.gsea.up<-GSE106090.DEG_limma_brown.gsea.up[GSE106090.DEG_limma_brown.gsea.up$ID %in% co.DESeq.gsea.up,]
GSE106090.DEG_limma_brown.gsea.down<-GSE106090.DEG_limma_brown.gsea.down[GSE106090.DEG_limma_brown.gsea.down$ID %in% co.DESeq.gsea.down,]
library("enrichplot")
up.pathways<-GSEA_up$ID[1:10]
down.pathways<-GSEA_down$ID[1:10]
gseaplot2(test[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], up.pathways, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# gseaplot2(test1[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], up.pathways, color = colorspace::rainbow_hcl(4))
ggsave("GSEA.GSE33774.DEG_limma_brown.up.pdf", width = 12, height = 8)
ggsave("GSEA.GSE33774.DEG_limma_brown.up.png", width = 12, height = 8)
gseaplot2(test[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], down.pathways, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# gseaplot2(test1[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], down.pathways, color = colorspace::rainbow_hcl(4))
ggsave("GSEA.GSE33774.DEG_limma_brown.down.pdf", width = 12, height = 8)
ggsave("GSEA.GSE33774.DEG_limma_brown.down.png", width = 12, height = 8)


gseaplot2(test1[["allpathway2gene"]][["GSE106090.DEG_limma_brown"]], up.pathways.1, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# gseaplot2(test1[["allpathway2gene"]][["GSE106090.DEG_limma_brown"]], up.pathways, color = colorspace::rainbow_hcl(4))
ggsave("GSEA.GSE106090.DEG_limma_brown.up.pdf", width = 12, height = 8)
ggsave("GSEA.GSE106090.DEG_limma_brown.up.png", width = 12, height = 8)
gseaplot2(test1[["allpathway2gene"]][["GSE106090.DEG_limma_brown"]], down.pathways.1, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# gseaplot2(test1[["allpathway2gene"]][["GSE106090.DEG_limma_brown"]], down.pathways, color = colorspace::rainbow_hcl(4))
ggsave("GSEA.GSE106090.DEG_limma_brown.down.pdf", width = 12, height = 8)
ggsave("GSEA.GSE106090.DEG_limma_brown.down.png", width = 12, height = 8)

####### plot GO ######
load("/Volumes/B/xjqqrza/implantitis/implantitis.brown.DEres2.Enrichment.Rdata")
GO.GSE33774.DEG_limma_brown<-test2[["go2gene"]][["GSE33774.DEG_limma_brown"]]@result
GO.GSE106090.DEG_limma_brown<-test2[["go2gene"]][["GSE106090.DEG_limma_brown"]]@result
write.csv(GO.GSE33774.DEG_limma_brown,file = "GO.GSE33774.DEG_limma_brown.csv")
write.csv(GO.GSE106090.DEG_limma_brown,file = "GO.GSE106090.DEG_limma_brown.csv")

# GSE33774.DEG_limma_brown
GO.GSE33774.DEG_limma_brown<-GO.GSE33774.DEG_limma_brown[GO.GSE33774.DEG_limma_brown$pvalue<0.05,]
GO.GSE33774.DEG_limma_brown<-GO.GSE33774.DEG_limma_brown[GO.GSE33774.DEG_limma_brown$qvalue<0.2,]
GO.GSE33774.DEG_limma_brown.cc<-GO.GSE33774.DEG_limma_brown[substr(GO.GSE33774.DEG_limma_brown$ID,1,4)=="GOCC",]
GO.GSE33774.DEG_limma_brown.mf<-GO.GSE33774.DEG_limma_brown[substr(GO.GSE33774.DEG_limma_brown$ID,1,4)=="GOMF",]
GO.GSE33774.DEG_limma_brown.bp<-GO.GSE33774.DEG_limma_brown[substr(GO.GSE33774.DEG_limma_brown$ID,1,4)=="GOBP",]
write.csv(GO.GSE33774.DEG_limma_brown.cc,file = "GO.GSE33774.DEG_limma_brown.cc.csv")
write.csv(GO.GSE33774.DEG_limma_brown.mf,file = "GO.GSE33774.DEG_limma_brown.mf.csv")
write.csv(GO.GSE33774.DEG_limma_brown.bp,file = "GO.GSE33774.DEG_limma_brown.bp.csv")
# GSE106090.DEG_limma_brown
GO.GSE106090.DEG_limma_brown<-GO.GSE106090.DEG_limma_brown[GO.GSE106090.DEG_limma_brown$pvalue<0.05,]
GO.GSE106090.DEG_limma_brown<-GO.GSE106090.DEG_limma_brown[GO.GSE106090.DEG_limma_brown$qvalue<0.2,]
GO.GSE106090.DEG_limma_brown.cc<-GO.GSE106090.DEG_limma_brown[substr(GO.GSE106090.DEG_limma_brown$ID,1,4)=="GOCC",]
GO.GSE106090.DEG_limma_brown.mf<-GO.GSE106090.DEG_limma_brown[substr(GO.GSE106090.DEG_limma_brown$ID,1,4)=="GOMF",]
GO.GSE106090.DEG_limma_brown.bp<-GO.GSE106090.DEG_limma_brown[substr(GO.GSE106090.DEG_limma_brown$ID,1,4)=="GOBP",]
write.csv(GO.GSE106090.DEG_limma_brown.cc,file = "GO.GSE106090.DEG_limma_brown.cc.csv")
write.csv(GO.GSE106090.DEG_limma_brown.mf,file = "GO.GSE106090.DEG_limma_brown.mf.csv")
write.csv(GO.GSE106090.DEG_limma_brown.bp,file = "GO.GSE106090.DEG_limma_brown.bp.csv")

cc<-intersect(GO.GSE33774.DEG_limma_brown.cc$ID,GO.GSE106090.DEG_limma_brown.cc$ID)
mf<-intersect(GO.GSE33774.DEG_limma_brown.mf$ID,GO.GSE106090.DEG_limma_brown.mf$ID)
bp<-intersect(GO.GSE33774.DEG_limma_brown.bp$ID,GO.GSE106090.DEG_limma_brown.bp$ID)
up.cc1<-GO.GSE33774.DEG_limma_brown.cc[GO.GSE33774.DEG_limma_brown.cc$ID %in% cc,]
up.cc1<-up.cc1$ID[1:10]
up.mf1<-GO.GSE33774.DEG_limma_brown.mf[GO.GSE33774.DEG_limma_brown.mf$ID %in% mf,]
up.mf1<-up.mf1$ID[1:10]
up.bp1<-GO.GSE33774.DEG_limma_brown.bp[GO.GSE33774.DEG_limma_brown.bp$ID %in% bp,]
up.bp1<-up.bp1$ID[1:10]

up.cc1<-GO.GSE33774.DEG_limma_brown.cc$ID[1:10]
up.mf1<-GO.GSE33774.DEG_limma_brown.mf$ID[1:10]
up.bp1<-GO.GSE33774.DEG_limma_brown.bp$ID[1:10]

library("enrichplot")
library(ggplot2)
kk<-test2[["go2gene"]][["GSE33774.DEG_limma_brown"]]
dotplot(kk,showCategory=up.cc1)
ggsave("GSE33774.DEG_limma_brown.cc.pdf", width = 12, height = 8)
dotplot(kk,showCategory=up.mf1)
ggsave("GSE33774.DEG_limma_brown.mf.pdf", width = 12, height = 8)
dotplot(kk,showCategory=up.bp1)
ggsave("GSE33774.DEG_limma_brown.bp.pdf", width = 12, height = 8)
kk<-test2[["go2gene"]][["GSE106090.DEG_limma_brown"]]
dotplot(kk,showCategory=up.cc1)
ggsave("GSE106090.DEG_limma_brown.cc.pdf", width = 12, height = 8)
dotplot(kk,showCategory=up.mf1)
ggsave("GSE106090.DEG_limma_brown.mf.pdf", width = 12, height = 8)
dotplot(kk,showCategory=up.bp1)
ggsave("GSE106090.DEG_limma_brown.bp.pdf", width = 12, height = 8)

####kegg#####
GSE33774.DEG_limma_brown.kegg<-(test2[["kegg2gene"]][["GSE33774.DEG_limma_brown"]]@result)
GSE106090.DEG_limma_brown.kegg<-test2[["kegg2gene"]][["GSE106090.DEG_limma_brown"]]@result
write.csv(GSE33774.DEG_limma_brown.kegg,file = "GSE33774.DEG_limma_brown.kegg.csv")
write.csv(GSE106090.DEG_limma_brown.kegg,file = "GSE106090.DEG_limma_brown.kegg.csv")
GSE33774.DEG_limma_brown.kegg<-GSE33774.DEG_limma_brown.kegg[GSE33774.DEG_limma_brown.kegg$pvalue<0.05,]
GSE106090.DEG_limma_brown.kegg<-GSE106090.DEG_limma_brown.kegg[GSE106090.DEG_limma_brown.kegg$pvalue<0.05,]
co.DESeq.kegg<-intersect(GSE106090.DEG_limma_brown.kegg$ID,GSE33774.DEG_limma_brown.kegg$ID)
GSE106090.DEG_limma_brown<-GSE106090.DEG_limma_brown.kegg.up[GSE106090.DEG_limma_brown.kegg.up$ID %in% co.DESeq.kegg.up,]
GSE106090.DEG_limma_brown<-GSE106090.DEG_limma_brown.kegg.down[GSE106090.DEG_limma_brown.kegg.down$ID %in% co.DESeq.kegg.down,]
library("enrichplot")
up.pathways<-kegg_up$ID[1:10]
down.pathways<-kegg_down$ID[1:10]
keggplot2(test[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], up.pathways, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# keggplot2(test1[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], up.pathways, color = colorspace::rainbow_hcl(4))
ggsave("kegg.GSE33774.DEG_limma_brown.up.pdf", width = 12, height = 8)
keggplot2(test[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], down.pathways, color = colorspace::rainbow_hcl(4),pvalue_table = TRUE)
# keggplot2(test1[["allpathway2gene"]][["GSE33774.DEG_limma_brown"]], down.pathways, color = colorspace::rainbow_hcl(4))
ggsave("kegg.GSE33774.DEG_limma_brown.down.pdf", width = 12, height = 8)
ggsave("kegg.GSE33774.DEG_limma_brown.down.png", width = 12, height = 8)
