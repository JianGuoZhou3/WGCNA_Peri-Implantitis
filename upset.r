library(ComplexHeatmap)
# install.packages("UpSetR")
library("UpSetR")
movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
library(readr)
upset <- read_csv("upset.csv")
u1<-union(upset$Closeness,upset$Between)
u2<-union(upset$Degree,upset$Stress)
u3<-union(upset$Radiality,upset$MNC)
u4<-union(upset$MCC,upset$EPC)
u5<-union(upset$EcCentricity,upset$DMNC)
u6<-union(upset$ClusteringCoefficient,upset$BottleNeck)
u7<-union(u1,u2)
u8<-union(u3,u4)
u9<-union(u5,u6)
u10<-union(u7,u8)
u11<-union(u9,u10)
upset1<-data.frame(u11)
upset1<-cbind(upset1,upset1,upset1,upset1)
rownames(upset1)<-upset1$u11
for (i in 1:ncol(upset)) {
  for (b in 1:143) {
    upset1[b,i]<-grepl(upset1[b,16],c(upset[,i]))
  }
}
upset2<-upset1[,1:12]
colnames(upset2)<-colnames(upset)
for (i in 1:ncol(upset2)) {
  upset2[,i]<- ifelse(upset2[,i]=="TRUE",1,0)
}
lt<-data.frame(upset)
m1 = make_comb_mat(lt)
UpSet(upset2)
upset(movies, attribute.plots=list(gridrows = 100, ncols = 1, 
                                   plots = list(list(plot=histogram, x="AvgRating",queries=T),
                                                list(plot = scatter_plot, y = "AvgRating", x = "Watches", queries = T))), 
      sets = c("Action", "Adventure", "Children", "War", "Noir"),
      queries = list(list(query = intersects, params = list("War"), active = T),
                     list(query = intersects, params = list("Noir"))))

upset(upset2, sets = colnames(upset2), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

####### boxplot #######
library(ggsignif)
library(ggplot2)
library(viridis)
library(ggpubr)
markers<-c("ACTB","TLR4","IL1B","STAT3","IL10","ITGAM","ITGB1")
id.sig<-markers
##for Train
# there are one case serverity1 have problem.
train.clin1$Severity1<-ifelse( train.clin1$Severity1=="moderate-severe","moderate",train.clin1$Severity1)
#GSM835253
train.clin1$Severity1<-ifelse( train.clin1$geo_accession=="GSM835253","severe",train.clin1$Severity1)
train.clin1<-train.clin1[train.clin1$Severity1!="none",]
train.clin1$Severity1<-ifelse(train.clin1$Severity1=="(healthy)","healthy",train.clin1$Severity1)
data.marker<-train.expr[rownames(train.expr) %in% markers,]
data.marker<-data.marker[,colnames(data.marker) %in% rownames(train.clin1)]
data.marker_1<-data.marker
# for (i in 1:nrow(data.marker)) {
#   x<-data.marker[i,]
#   x<-as.matrix(x)
#   data.marker_1[i,] <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  # Standardize manually
# }
data.marker_1<-t(scale(t(data.marker)))
sd(data.marker_1[1,])
data.marker1<-data.frame(t(data.marker_1))
data<-cbind(data.marker1,train.clin1)
data$Group = as.factor(data$Group)
data$gender = as.factor(data$gender)
data$age = as.factor(data$age)
data$Severity1 = as.factor(data$Severity1)
## for test
data.marker<-test.expr[rownames(test.expr) %in% markers,]
data.marker_1<-data.marker
# for (i in 1:nrow(data.marker)) {
#   x<-data.marker[i,]
#   x<-as.matrix(x)
#   data.marker_1[i,] <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  # Standardize manually
# }
data.marker_1<-t(scale(t(data.marker)))
sd(data.marker_1[1,])
data.marker1<-data.frame(t(data.marker_1))
test.clin1$Severity1<-ifelse(test.clin1$Severity1=="(healthy)","healthy",test.clin1$Severity1)
data<-cbind(data.marker1,test.clin1)
data$Group = as.factor(data$Group)
# data$Gender = as.factor(data$Gender)
data$age = as.factor(data$age)
data$Severity1 = as.factor(data$Severity1)

#### boxplot#######
expr.sig.plot<-NULL
clins<-c("Group","Gender","Severity1","age")
b=3
for (b in 1:length(clins)) {
  expr.sig.plot<-NULL
  clin.compare<-clins[b]
  for (i in 1:length(id.sig)) {
    marker<-id.sig[i]
    expr.sig.plot1<-data[,c("Sample_ID",clin.compare,marker)]
    expr.sig.plot1$sig<-marker
    expr.sig.plot1<-data.frame(expr.sig.plot1)
    colnames(expr.sig.plot1)<-c("Sample_ID",clin.compare,"expression","sig")
    expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
    a<-expr.sig.plot[expr.sig.plot$sig==marker,]
    p <- ggboxplot(a, x = clin.compare, y = "expression",
                   color = clin.compare, palette = "jco",title=marker,
                   add = "jitter")+ stat_compare_means()#  Add p-value
    # Boxplots are automatically dodged when any aesthetic is a factor
    # wilcox.test/t.test
    ggsave(p, filename = paste0(marker, "_", clin.compare, 
                                #"_test_GSE106090_ggplot2.pdf"), width = 6, height = 6)
                                "_train_GSE33774_ggplot2.pdf"), width = 6, height = 6)
  }
  rm(expr.sig.plot1)
}

###### for train data$Gender ######
data<-data[!(data$Gender=="none"),]
for (b in 2:2) {
  expr.sig.plot<-NULL
  clin.compare<-clins[b]
  for (i in 1:length(id.sig)) {
    marker<-id.sig[i]
    expr.sig.plot1<-data[,c("Sample_ID",clin.compare,marker)]
    expr.sig.plot1$sig<-marker
    expr.sig.plot1<-data.frame(expr.sig.plot1)
    colnames(expr.sig.plot1)<-c("Sample_ID",clin.compare,"expression","sig")
    expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
    a<-expr.sig.plot[expr.sig.plot$sig==marker,]
    p <- ggboxplot(a, x = clin.compare, y = "expression",
                   color = clin.compare, palette = "jco",title=marker,
                   add = "jitter")+ stat_compare_means()#  Add p-value
    # Boxplots are automatically dodged when any aesthetic is a factor
    # wilcox.test/t.test
    ggsave(p, filename = paste0(marker, "_", clin.compare, 
                                "_train_GSE33774_ggplot2.pdf"), width = 6, height = 6)
  }
  rm(expr.sig.plot1)
}
###### for train data$Severity1 ######
data<-data[!(data$Severity1=="none"),]
for (b in 3:3) {
  expr.sig.plot<-NULL
  clin.compare<-clins[b]
  for (i in 1:length(id.sig)) {
    marker<-id.sig[i]
    expr.sig.plot1<-data[,c("Sample_ID",clin.compare,marker)]
    expr.sig.plot1$sig<-marker
    expr.sig.plot1<-data.frame(expr.sig.plot1)
    colnames(expr.sig.plot1)<-c("Sample_ID",clin.compare,"expression","sig")
    expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
    a<-expr.sig.plot[expr.sig.plot$sig==marker,]
    p <- ggboxplot(a, x = clin.compare, y = "expression",
                   color = clin.compare, palette = "jco",title=marker,
                   add = "jitter")+ stat_compare_means()#  Add p-value
    # Boxplots are automatically dodged when any aesthetic is a factor
    # wilcox.test/t.test
    ggsave(p, filename = paste0(marker, "_", clin.compare, 
                                "_train_GSE33774_ggplot2.pdf"), width = 6, height = 6)
  }
  rm(expr.sig.plot1)
}


#######violin plot#########
library(ggplot2)
library(ggsci)
colnames(GSE33774.TIL.28_ssgsea)
# [1] "ID"                              "Index"                          
# [3] "Activated.CD8.T.cell"            "Central.memory.CD8.T.cell"      
# [5] "Effector.memeory.CD8.T.cell"     "Activated.CD4.T.cell"           
# [7] "Central.memory.CD4.T.cell"       "Effector.memeory.CD4.T.cell"    
# [9] "T.follicular.helper.cell"        "Gamma.delta.T.cell"             
# [11] "Type.1.T.helper.cell"            "Type.17.T.helper.cell"          
# [13] "Type.2.T.helper.cell"            "Regulatory.T.cell"              
# [15] "Activated.B.cell"                "Immature..B.cell"               
# [17] "Memory.B.cell"                   "Natural.killer.cell"            
# [19] "CD56bright.natural.killer.cell"  "CD56dim.natural.killer.cell"    
# [21] "Myeloid.derived.suppressor.cell" "Natural.killer.T.cell"          
# [23] "Activated.dendritic.cell"        "Plasmacytoid.dendritic.cell"    
# [25] "Immature.dendritic.cell"         "Macrophage"                     
# [27] "Eosinophil"                      "Mast.cell"                      
# [29] "Monocyte"                        "Neutrophil" 

innate_immmunity<-c("Activated.dendritic.cell","CD56bright.natural.killer.cell","CD56dim.natural.killer.cell",
                    "Eosinophil","Immature.dendritic.cell","Plasmacytoid.dendritic.cell",
                    "Macrophage","Mast.cell","Myeloid.derived.suppressor.cell",
                    "Monocyte","Natural.killer.cell","Neutrophil","Natural.killer.T.cell")
adaptive_immmunity<-c("Activated.B.cell","Immature..B.cell","Memory.B.cell","Activated.CD4.T.cell",
                      "Activated.CD8.T.cell","Central.memory.CD4.T.cell","Central.memory.CD8.T.cell",
                      "Effector.memeory.CD4.T.cell","Effector.memeory.CD8.T.cell","Gamma.delta.T.cell","Regulatory.T.cell",
                      "T.follicular.helper.cell","Type.1.T.helper.cell","Type.17.T.helper.cell",
                      "Type.2.T.helper.cell")
#clin
GSE33774.TIL.28_ssgsea<-GSE33774.TIL.28_ssgsea[order(rownames(GSE33774.TIL.28_ssgsea)),]
train.clin1<-train.clin1[order(rownames(train.clin1)),]
data<-cbind(GSE33774.TIL.28_ssgsea,train.clin1)
GSE106090.TIL.28_ssgsea<-GSE106090.TIL.28_ssgsea[order(rownames(GSE106090.TIL.28_ssgsea)),]
test.clin1<-test.clin1[order(rownames(test.clin1)),]
data<-cbind(GSE106090.TIL.28_ssgsea,test.clin1)
data$Group<-as.factor(data$Group)
expr.sig.plot<-NULL
for (i in 3:length(colnames(GSE33774.TIL.28_ssgsea))) {
  marker<-colnames(GSE33774.TIL.28_ssgsea)[i]
  expr.sig.plot1<-data[,c("ID","Group",marker)]
  expr.sig.plot1$sig<-marker
  expr.sig.plot1<-data.frame(expr.sig.plot1)
  colnames(expr.sig.plot1)<-c("Sample_ID","Group","expression","sig")
  expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
  a<-expr.sig.plot[expr.sig.plot$sig==marker,]
  p <- ggviolin(a, "Group", "expression",fill="Group",
            palette = "jama", add = "boxplot", add.params = list(fill = "white"),
                title=marker)+ stat_compare_means()#  Add p-value
  # Boxplots are automatically dodged when any aesthetic is a factor
  # wilcox.test/t.test
  ggsave(p, filename = paste0(marker, "_", "disease.group", 
                              # "_train_GSE33774_ggplot2.pdf"), width = 6, height = 6)
 "_test_GSE106090_ggplot2.pdf"), width = 6, height = 6)

}
rm(expr.sig.plot1)
GSE33774.TIL.all<-expr.sig.plot
GSE106090.TIL.all<-expr.sig.plot

#####heatmap of cor#######

markers<-c("ACTB","TLR4","IL1B","STAT3","IL10","ITGAM","ITGB1")
### corr
GSE33774.innate<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE33774.innate<-data.frame(GSE33774.innate)
colnames(GSE33774.innate)<-innate_immmunity
rownames(GSE33774.innate)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.innate[b1,i1]<-d1[,4]
  }
}

GSE33774.adaptive<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE33774.adaptive<-data.frame(GSE33774.adaptive)
colnames(GSE33774.adaptive)<-adaptive_immmunity
rownames(GSE33774.adaptive)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.adaptive[b1,i1]<-d1[,4]
  }
}

GSE106090.innate<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE106090.innate<-data.frame(GSE106090.innate)
colnames(GSE106090.innate)<-innate_immmunity
rownames(GSE106090.innate)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.innate[b1,i1]<-d1[,4]
  }
}

GSE106090.adaptive<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE106090.adaptive<-data.frame(GSE106090.adaptive)
colnames(GSE106090.adaptive)<-adaptive_immmunity
rownames(GSE106090.adaptive)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.adaptive[b1,i1]<-d1[,4]
  }
}
### p
GSE33774.innate.p<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE33774.innate.p<-data.frame(GSE33774.innate.p)
colnames(GSE33774.innate.p)<-innate_immmunity
rownames(GSE33774.innate.p)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.innate.p[b1,i1]<-d1[,3]
  }
}

GSE33774.adaptive.p<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE33774.adaptive.p<-data.frame(GSE33774.adaptive.p)
colnames(GSE33774.adaptive.p)<-adaptive_immmunity
rownames(GSE33774.adaptive.p)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.adaptive.p[b1,i1]<-d1[,3]
  }
}

GSE106090.innate.p<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE106090.innate.p<-data.frame(GSE106090.innate.p)
colnames(GSE106090.innate.p)<-innate_immmunity
rownames(GSE106090.innate.p)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.innate.p[b1,i1]<-d1[,3]
  }
}

GSE106090.adaptive.p<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE106090.adaptive.p<-data.frame(GSE106090.adaptive.p)
colnames(GSE106090.adaptive.p)<-adaptive_immmunity
rownames(GSE106090.adaptive.p)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.adaptive.p[b1,i1]<-d1[,3]
  }
}
# Heatmap
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
mat<-GSE33774.innate
mat<-GSE33774.adaptive
mat<-GSE106090.innate
mat<-GSE106090.adaptive

pdf(file="GSE33774.innate_heatmap.pdf",width=10,height=6)
pdf(file="GSE33774.adaptive_heatmap.pdf",width=10,height=6)
pdf(file="GSE106090.innate_heatmap.pdf",width=10,height=6)
pdf(file="GSE106090.adaptive_heatmap.pdf",width=10,height=6)
mat<-round(mat,2)
Heatmap(mat,cluster_rows = FALSE,
        cluster_columns = FALSE,c("blue", "#EEEEEE", "red"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(mat[i, j], x, y)
        })
dev.off()


mat<-GSE33774.innate.p
mat<-GSE33774.adaptive.p
mat<-GSE106090.innate.p
mat<-GSE106090.adaptive.p

pdf(file="GSE33774.innate.p_heatmap.pdf",width=10,height=6)
pdf(file="GSE33774.adaptive.p_heatmap.pdf",width=10,height=6)
pdf(file="GSE106090.innate.p_heatmap.pdf",width=10,height=6)
pdf(file="GSE106090.adaptive.p_heatmap.pdf",width=10,height=6)
mat<-signif(mat,2) # for p value,eg. 1.33-11
mat<-round(mat,2) # for p value,eg. 1.33-11
Heatmap(mat,cluster_rows = FALSE,
        cluster_columns = FALSE,c("white", "gray"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(mat[i, j], x, y)
        })
dev.off()

#######violin plot 2021.12.12 updated#########
library(ggplot2)
library(ggsci)
library(ggpubr)
colnames(GSE33774.TIL.28_ssgsea)
# deleted periodontitis
GSE33774.TIL.28_ssgsea<-GSE33774.TIL.28_ssgsea[order(rownames(GSE33774.TIL.28_ssgsea)),]
train.clin1<-train.clin1[order(rownames(train.clin1)),]
data<-cbind(GSE33774.TIL.28_ssgsea,train.clin1)
data<-data[data$Group!="Periodontitis",]
GSE106090.TIL.28_ssgsea<-GSE106090.TIL.28_ssgsea[order(rownames(GSE106090.TIL.28_ssgsea)),]
test.clin1<-test.clin1[order(rownames(test.clin1)),]
data<-cbind(GSE106090.TIL.28_ssgsea,test.clin1)
data<-data[data$Group!="Periodontitis",]
data$Group<-as.factor(data$Group)
expr.sig.plot<-NULL
for (i in 3:length(colnames(GSE33774.TIL.28_ssgsea))) {
  marker<-colnames(GSE33774.TIL.28_ssgsea)[i]
  expr.sig.plot1<-data[,c("ID","Group",marker)]
  expr.sig.plot1$sig<-marker
  expr.sig.plot1<-data.frame(expr.sig.plot1)
  colnames(expr.sig.plot1)<-c("Sample_ID","Group","expression","sig")
  expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
  a<-expr.sig.plot[expr.sig.plot$sig==marker,]
  p <- ggviolin(a, "Group", "expression",fill="Group",
                palette = "jama", add = "boxplot", add.params = list(fill = "white"),
                title=marker)+ stat_compare_means()#  Add p-value
  # Boxplots are automatically dodged when any aesthetic is a factor
  # wilcox.test/t.test
  ggsave(p, filename = paste0(marker, "_", "disease.group", 
                              #"_train_GSE33774_ggplot2.updated.pdf"), width = 6, height = 6)
                              "_test_GSE106090_ggplot2.updated.pdf"), width = 6, height = 6)
}
  
#####heatmap of cor 2021.12.12 updated#######
####correlation  #######
TILs<-GSE33774.TIL.28_ssgsea[,3:30]
TILs<-TILs[rownames(TILs) %in% train.clin2$Sample_ID,]
data.marker<-train.expr[rownames(train.expr) %in% markers,]
train.clin2<-train.clin1[train.clin1$Group!="Periodontitis",]
data.marker<-data.marker[,colnames(data.marker) %in%  train.clin2$Sample_ID ]
TILs<-GSE106090.TIL.28_ssgsea[,3:30]
data.marker<-test.expr[rownames(test.expr) %in% markers,]
test.clin2<-test.clin1[test.clin1$Group!="Periodontitis",]
TILs<-TILs[rownames(TILs) %in% test.clin2$Sample_ID,]
data.marker<-data.marker[,colnames(data.marker) %in%  test.clin2$Sample_ID ]


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
GSE33774.p1.updated<-p1
GSE106090.p1.updated<-p1
write.csv(GSE33774.p1.updated,file = "GSE33774.p1.updated.csv")
write.csv(GSE106090.p1.updated,file = "GSE106090.p1.updated.csv")

markers<-c("ACTB","TLR4","IL1B","STAT3","IL10","ITGAM","ITGB1")
### corr
GSE33774.innate<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE33774.innate<-data.frame(GSE33774.innate)
colnames(GSE33774.innate)<-innate_immmunity
rownames(GSE33774.innate)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.innate[b1,i1]<-d1[,4]
  }
}

GSE33774.adaptive<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE33774.adaptive<-data.frame(GSE33774.adaptive)
colnames(GSE33774.adaptive)<-adaptive_immmunity
rownames(GSE33774.adaptive)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.adaptive[b1,i1]<-d1[,4]
  }
}

GSE106090.innate<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE106090.innate<-data.frame(GSE106090.innate)
colnames(GSE106090.innate)<-innate_immmunity
rownames(GSE106090.innate)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.innate[b1,i1]<-d1[,4]
  }
}

GSE106090.adaptive<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE106090.adaptive<-data.frame(GSE106090.adaptive)
colnames(GSE106090.adaptive)<-adaptive_immmunity
rownames(GSE106090.adaptive)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.adaptive[b1,i1]<-d1[,4]
  }
}
### p
GSE33774.innate.p<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE33774.innate.p<-data.frame(GSE33774.innate.p)
colnames(GSE33774.innate.p)<-innate_immmunity
rownames(GSE33774.innate.p)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.innate.p[b1,i1]<-d1[,3]
  }
}

GSE33774.adaptive.p<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE33774.adaptive.p<-data.frame(GSE33774.adaptive.p)
colnames(GSE33774.adaptive.p)<-adaptive_immmunity
rownames(GSE33774.adaptive.p)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE33774.p1[GSE33774.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE33774.adaptive.p[b1,i1]<-d1[,3]
  }
}

GSE106090.innate.p<-matrix(nrow =length(markers),ncol = length(innate_immmunity) )
GSE106090.innate.p<-data.frame(GSE106090.innate.p)
colnames(GSE106090.innate.p)<-innate_immmunity
rownames(GSE106090.innate.p)<-markers
for (i in 1:length(innate_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-innate_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.innate.p[b1,i1]<-d1[,3]
  }
}

GSE106090.adaptive.p<-matrix(nrow =length(markers),ncol = length(adaptive_immmunity) )
GSE106090.adaptive.p<-data.frame(GSE106090.adaptive.p)
colnames(GSE106090.adaptive.p)<-adaptive_immmunity
rownames(GSE106090.adaptive.p)<-markers
for (i in 1:length(adaptive_immmunity)) {
  for (b in 1:length(markers)) {
    i1<-adaptive_immmunity[i]
    b1<-markers[b]
    d1<-GSE106090.p1[GSE106090.p1$id.1==b1,]
    d1<-d1[d1$id.2==i1,]
    GSE106090.adaptive.p[b1,i1]<-d1[,3]
  }
}
# Heatmap
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
mat<-GSE33774.innate
mat<-GSE33774.adaptive
mat<-GSE106090.innate
mat<-GSE106090.adaptive

pdf(file="GSE33774.innate_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE33774.adaptive_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE106090.innate_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE106090.adaptive_heatmap.updated.pdf",width=10,height=6)
mat<-round(mat,2)
Heatmap(mat,cluster_rows = FALSE,
        cluster_columns = FALSE,c("blue", "#EEEEEE", "red"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(mat[i, j], x, y)
        })
dev.off()


mat<-GSE33774.innate.p
mat<-GSE33774.adaptive.p
mat<-GSE106090.innate.p
mat<-GSE106090.adaptive.p

pdf(file="GSE33774.innate.p_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE33774.adaptive.p_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE106090.innate.p_heatmap.updated.pdf",width=10,height=6)
pdf(file="GSE106090.adaptive.p_heatmap.updated.pdf",width=10,height=6)
mat<-signif(mat,2) # for p value,eg. 1.33-11
mat<-round(mat,2) # for p value,eg. 1.33-11
Heatmap(mat,cluster_rows = FALSE,
        cluster_columns = FALSE,c("white", "gray"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(mat[i, j], x, y)
        })
dev.off()

#### 3.6 Evaluation of the potential regulatory effect of the hub genes####
# We searched for the Spearman correlations of hub genes with the expression level of 38 immune-related factors in TISIDB database, 
# including immune-inhibitory factors, immune-stimulatory factors, chemokines, and receptors. 
# These correlations are presented as heatmaps (Figure 6A).

##### heatmap for brown ######
load("/Volumes/B/xjqqrza/implantitis/implantitis.brown.DEres2.Rdata")
deg<-DEres2[["GSE33774.DEG_limma_brown"]]
deg<-deg[deg$pvalue<0.05,]
deg<-deg[deg$padj<0.05,]
deg.1<-DEres2[["GSE106090.DEG_limma_brown"]]
deg.1<-deg.1[deg.1$pvalue<0.05,]
deg.1<-deg.1[deg.1$padj<0.05,]
id<-intersect(deg.1$genesymbol,deg$genesymbol)
#GSE33774
library(readr)
brown <- read_csv("brown.csv")
brown<-data.frame(brown)
rownames(brown)<-brown$...1
brown<-brown[,-1]
# GSE106090
brown <- read_csv("brown.test.csv")
brown<-data.frame(brown)
rownames(brown)<-brown$...1
brown<-brown[,-1]
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
pdf(file="GSE33774.brown_heatmap.pdf",width=8,height=12)
pdf(file="GSE106090.brown_heatmap.pdf",width=8,height=12)
mat<-brown[rownames(brown) %in% id,]
mat<-scale(t(mat))
mat<-t(mat)
Heatmap(mat,cluster_rows = T,
        cluster_columns = T,col = col_fun)
dev.off()

#######Functional analysis ########
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  install.packages("org.Hs.eg.db")
if (!requireNamespace("db", quietly = TRUE))
  install.packages("db")
if (!requireNamespace("db", quietly = TRUE))
  install.packages("db")
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

deg<-DEres2[["GSE33774.DEG_limma_brown"]]
deg<-deg[deg$pvalue<0.05,]
deg<-deg[deg$padj<0.2,]
# class(str_split(ensemb_genes,'[.]',simplify = T))
# class(unlist(str_split(ensemb_genes,'[.]')))
# deg$ensembl_id=str_split(ensemb_genes,'[.]',simplify = T)[,1]
library(org.Hs.eg.db)
d=deg[order(deg$log2FoldChange),]
diffg= d[abs(d$log2FoldChange)>0.5,] 
id.fc=diffg[,c(7,1)]
colnames(id.fc)=c("entrezgene_id","logFC")
#########
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
ego_BP <- enrichGO(gene = id.fc$entrezgene_id,
                   #小鼠用这行
                   # OrgDb = org.Mm.eg.db,
                   #人类用这行
                   OrgDb = org.Hs.eg.db,
                   #非模式生物用这行，例如玉米
                   #OrgDb = maize.db,
                   ont = "BP", #或MF或CC
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
                   # qvalueCutoff  = 0.05
) 

dim(ego_BP)
ego_MF <- enrichGO(gene = id.fc$entrezgene_id,
                   #小鼠用这行
                   # OrgDb = org.Mm.eg.db,
                   #人类用这行
                   OrgDb = org.Hs.eg.db,
                   #非模式生物用这行，例如玉米
                   #OrgDb = maize.db,
                   ont = "MF", #或MF或CC
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
                   # qvalueCutoff  = 0.01
) 
dim(ego_MF)
ego_CC <- enrichGO(gene = id.fc$entrezgene_id,
                   #小鼠用这行
                   # OrgDb = org.Mm.eg.db,
                   #人类用这行
                   OrgDb = org.Hs.eg.db,
                   #非模式生物用这行，例如玉米
                   #OrgDb = maize.db,
                   ont = "CC", #或MF或CC
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
                   # qvalueCutoff  = 0.01
) 
dim(ego_CC)
#把富集分析结果保存到文件
save(ego_BP,ego_CC,ego_MF,file = "enrichGO_BP_output_0.05.rdata")
write.csv(ego_BP,"enrichGO_BP_output_0.05.csv",quote = F)
write.csv(ego_CC,"enrichGO_CC_output_0.05.csv",quote = F)
write.csv(ego_MF,"enrichGO_MF_output_0.05.csv",quote = F)
barplot(ego_BP, drop=TRUE, showCategory=12)
ggsave("enrichGO_0.05_clusterProfiler_barplot.pdf", width = 12, height = 8)
barplot(ego_CC, drop=TRUE, showCategory=12)
ggsave("enrichGO_0.05_CC_clusterProfiler_barplot.pdf", width = 12, height = 8)
barplot(ego_MF, drop=TRUE, showCategory=12)
ggsave("enrichGO_0.05_MF_clusterProfiler_barplot.pdf", width = 12, height = 8)
###kegg####
##https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#kegg-analysis
geneList <-deg$logFC
names(geneList) <- deg$gene_id
geneList=sort(geneList,decreasing = T)
deg=unique(deg)
##用logFC大于1的
id.fc=id.fc[order(id.fc$entrezgene_id,id.fc$logFC,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
id.fc=id.fc[!duplicated(id.fc$entrezgene_id),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s

geneList <-id.fc$logFC
names(geneList) <- id.fc$entrezgene_id
geneList=sort(geneList,decreasing = T)
save(geneList,file = "geneList.Rdata")
rm(list = ls())
# KEGG Gene Set Enrichment Analysis
kk <- gseKEGG(geneList, organism = "hsa", keyType = "kegg", exponent = 1,
              nPerm = 10000, minGSSize = 10, maxGSSize = 500,
              pvalueCutoff = 0.25, pAdjustMethod = "BH", verbose = TRUE,
              use_internal_data = FALSE, seed = FALSE, by = "DOSE")
head(kk)
save(kk,file="gseKEGG.2020.06.10.cut.0.5.Rdata")
write.csv(kk,"gseKEGG.2020.06.10.csv")
gseaplot(kk, geneSetID = "hsa04970")
dotplot(kk,showCategory=20)
ggsave("gseKEGG.2020.06.10.pdf", width = 12, height = 8)

browseKEGG(kk, 'hsa04970')
library("pathview")
hsa04722 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04722",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
# 这一条语句就做完了KEGG的GSEA分析。

sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
write.table(sortkk,"gsea_output.txt",sep = "\t",quote = F,col.names = T,row.names = F)

### GSEA可视化 用！！
library("enrichplot")

paths <- c("hsa04970","hsa04972","hsa04740","hsa05206")
gseaplot2(kk, paths)
gseaplot2(kk, paths, color = colorspace::rainbow_hcl(4))



#######violin plot 2021.12.17 updated#########
library(ggplot2)
library(ggsci)
library(ggpubr)
colnames(GSE33774.TIL.28_ssgsea)
# deleted periodontitis
GSE33774.TIL.28_ssgsea<-GSE33774.TIL.28_ssgsea[order(rownames(GSE33774.TIL.28_ssgsea)),]
train.clin1<-train.clin1[order(rownames(train.clin1)),]
data<-cbind(GSE33774.TIL.28_ssgsea,train.clin1)
data<-data[data$Group!="Periodontitis",]
GSE106090.TIL.28_ssgsea<-GSE106090.TIL.28_ssgsea[order(rownames(GSE106090.TIL.28_ssgsea)),]
test.clin1<-test.clin1[order(rownames(test.clin1)),]
data<-cbind(GSE106090.TIL.28_ssgsea,test.clin1)
data<-data[data$Group!="Periodontitis",]
data$Group<-as.factor(data$Group)
expr.sig.plot<-NULL
for (i in 3:length(colnames(GSE33774.TIL.28_ssgsea))) {
  marker<-colnames(GSE33774.TIL.28_ssgsea)[i]
  expr.sig.plot1<-data[,c("ID","Group",marker)]
  expr.sig.plot1$sig<-marker
  expr.sig.plot1<-data.frame(expr.sig.plot1)
  colnames(expr.sig.plot1)<-c("Sample_ID","Group","expression","sig")
  expr.sig.plot<-rbind(expr.sig.plot,expr.sig.plot1)
  a<-expr.sig.plot[expr.sig.plot$sig==marker,]
  # http://yuedong.site/learnR/graphic/%E5%9D%90%E6%A0%87%E8%BD%B4.html
  # # This removes all legends
  # bp + theme(legend.position="none")
  p <- ggviolin(a, "Group", "expression",fill="Group",
                palette = c("#3300CC","#FF3300"), add = "boxplot", add.params = list(fill = "white"),
                title=marker)+ ylim(-0.2, 1.0)+ theme(axis.title.y=element_blank())+ theme(legend.position="none")+ theme(axis.text.y = element_blank())+stat_compare_means(label = "p.signif", method = "wilcox.test",
                                                                    ref.group = "Healthy")                    # Pairwise comparison against reference
  
  # Boxplots are automatically dodged when any aesthetic is a factor
  # wilcox.test/t.test
  ggsave(p, filename = paste0(marker, "_", "disease.group", 
                              # "_train_GSE33774_ggplot2.updated.pdf"), width =1.65, height = 8.74)
                              "_test_GSE106090_ggplot2.updated.pdf"), width =1.65, height = 8.74)
}


#####ggscatterstats######
library(ggstatsplot)
#' citation("ggstatsplot")
#' 
#' Patil, I. (2021). Visualizations with statistical details: The
#' 'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
#' doi:10.21105/joss.03167
#' 
#' A BibTeX entry for LaTeX users is
#' 
#' @Article{,
#'   doi = {10.21105/joss.03167},
#'   url = {https://doi.org/10.21105/joss.03167},
#'   year = {2021},
#'   publisher = {{The Open Journal}},
#'   volume = {6},
#'   number = {61},
#'   pages = {3167},
#'   author = {Indrajeet Patil},
#'   title = {{Visualizations with statistical details: The {'ggstatsplot'} approach}},
#'   journal = {{Journal of Open Source Software}},
#' }
# This function creates a scatterplot with marginal distributions overlaid on the axes and results from statistical tests in the subtitle:
# install.packages('ggside')
TILs<-GSE33774.TIL.28_ssgsea[,3:30]
TILs<-TILs[rownames(TILs) %in% train.clin2$Sample_ID,]
data.marker<-train.expr[rownames(train.expr) %in% markers,]
train.clin2<-train.clin1[train.clin1$Group!="Periodontitis",]
data.marker<-data.marker[,colnames(data.marker) %in%  train.clin2$Sample_ID ]
TILs<-GSE106090.TIL.28_ssgsea[,3:30]
data.marker<-test.expr[rownames(test.expr) %in% markers,]
test.clin2<-test.clin1[test.clin1$Group!="Periodontitis",]
TILs<-TILs[rownames(TILs) %in% test.clin2$Sample_ID,]
data.marker<-data.marker[,colnames(data.marker) %in%  test.clin2$Sample_ID ]



markers<-c("ACTB","TLR4","IL1B","STAT3","IL10","ITGAM","ITGB1")
IPT.a<-c("Plasmacytoid.dendritic.cell","Macrophage","Myeloid.derived.suppressor.cell",
"Natural.killer.T.cell")
### corr
library(ggplot2)
TILsa<-TILs[,colnames(TILs) %in% IPT.a]
for (i in 1:length(markers)) {
  for (b in 1:length(IPT.a)) {
    x=data.frame(t(data.marker[i,]))
    y=(TILsa[,b])
    data<-cbind(x,y)
    colnames(data)<-c("a","b")
    p<-ggscatterstats(
      data  = data,
      x     = b,
      y     = a,
      xlab  = colnames(TILsa)[b],
      ylab  = markers[i]
      )
    ggsave(p, file=paste0("plot_GSE33774_",colnames(TILsa)[b],".vs.",markers[i],".pdf"), width = 8, height = 8)
  }
  }


# data  = ggplot2::msleep
# ggscatterstats(
#   data  = ggplot2::msleep,
#   x     = sleep_rem,
#   y     = awake,
#   xlab  = "",
#   ylab  = "Amount of time spent awake (in hours)",
#   title = ""
# )
##### heatmap for brown 7 markers ######
markers
#GSE33774
library(readr)
brown <- read_csv("brown.csv")
brown<-data.frame(brown)
rownames(brown)<-brown$...1
brown<-brown[,-1]
# GSE106090
brown <- read_csv("brown.test.csv")
brown<-data.frame(brown)
rownames(brown)<-brown$...1
brown<-brown[,-1]
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
pdf(file="GSE33774.brown_heatmap.markers.pdf",width=8,height=4)
pdf(file="GSE106090.brown_heatmap.markers.pdf",width=8,height=4)
mat<-brown[rownames(brown) %in% markers,]
mat<-scale(t(mat))
mat<-t(mat)
Heatmap(mat,cluster_rows = T,
        cluster_columns = T,col = col_fun)
dev.off()


