#### T cells only ####

Idents(T.clust) <- T.clust$annotv1
T.clust3 <- subset(T.clust, idents = c( "Naive CD8 T",
                                        "Memory CD8+ Tem",
                                        "Memory CD8 T_c46",
                                        "Memory CD8+ Trm",
                                        "Naive CD4 T",
                                        "Treg_c13",
                                        "Treg_c27",
                                        "Th1",
                                        "Th17",
                                        "Naive dnT_c40",
                                        "Activated dnT_35",
                                        "Activated dnT_c51",
                                        "gdT",
                                        "CCR6+ gdT",
                                        "T-innate-like_c49",
                                        "NKT_c17",
                                        "NKT_c21",
                                        "NKT_c30",
                                        "Activated NKT_c55"))


T.clust3$annotv1 <- factor(T.clust3$annotv1, levels = c( "Naive CD8 T",
                                                         "Memory CD8+ Tem",
                                                         "Memory CD8 T_c46",
                                                         "Memory CD8+ Trm",
                                                         "Naive CD4 T",
                                                         "Treg_c13",
                                                         "Treg_c27",
                                                         "Th1",
                                                         "Th17",
                                                         "Naive dnT_c40",
                                                         "Activated dnT_35",
                                                         "Activated dnT_c51",
                                                         "gdT",
                                                         "CCR6+ gdT",
                                                         "T-innate-like_c49",
                                                         "NKT_c17",
                                                         "NKT_c21",
                                                         "NKT_c30",
                                                         "Activated NKT_c55"))
Idents(T.clust3) <- T.clust3$annotv1



mytcol <- c( "Naive CD8 T" = "#B3B3FF",
             "Memory CD8+ Tem" = "#0000E6",
             "Memory CD8 T_c46" = "#0000B3",
             "Memory CD8+ Trm" = "#000080",
             "Naive CD4 T" = "#FF7F7F",
             "Treg_c13" = "#FF1A1A",
             "Treg_c27" = "#E60000",
             "Th1" = "#CC0000",
             "Th17" = "#B30000",
             "Naive dnT_c40" = "#FFFF66",
             "Activated dnT_35" = "#F2F200",
             "Activated dnT_c51" = "#D9D900",
             "gdT" = "#FF00A6",
             "CCR6+ gdT" = "#66000D",
             "T-innate-like_c49" = "#8e2f2f",
             "NKT_c17" = "#4DFF4D",
             "NKT_c21" = "#00CC00",
             "NKT_c30" = "#009900",
             "Activated NKT_c55" = "#008000")

DimPlot(T.clust3, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, cols=mytcol,) + 
  xlim(-7,7) + ylim(-3, 10) +ggtitle('T cells')
DimPlot(T.clust3, reduction = "umap", label = FALSE, pt.size = 0.5, repel = TRUE, cols=mytcol,) + 
  xlim(-7,7) + ylim(-3, 10)   #width = 900, height = 500



# Heat map, features
T3.markers <- FindAllMarkers(T.clust3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T.markers1<- T3.markers
T.markers1 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
T.markers1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

par(mar = c(4.5,4,1,1))

#Idents(T.clust3) <- forcats::fct_rev(factor(T.clust3$annot))

DoHeatmap(subset(T.clust3, downsample = 50), group.colors = mytcol, features = top20$gene, size = 3 ) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))



t.features = c("Cd3d", "Cd4", "Cd8a", "Cd8b1", "Trac", "Trdc", "Sell", "Ccr7", #CD4,8,gdt,naive,Dnt
               "Tnfrsf14", "Cxcr6", "Cx3cr1", "Itgax", "Nkg7",
               "Ifng", "Tbx21", "Ccr5", # Th1
               "Foxp3", "Icos", "Tnfrsf4", "Maf", "Rora", "Tnf", "Gzmk", "Rorc", "Ccr6",  #Treg, Th17, T-innate-like
               "Ifit1", "Ifit3", "Ifit203", "Irf4", "Ly6c2", "Itga1", "Klrk1",
               "Klra7", "Klra1", "Klre1", "Ccl3", "Ccl4", "Ccl5", "Xcl1"
)

#### T cells - dotplot (Supp Fig S7A) ####
Idents(T.clust3) <- forcats::fct_rev(factor(T.clust3$annotv1))
DotPlot(T.clust3, features = t.features) + RotatedAxis() + 
  geom_hline(yintercept=10.5, linetype="longdash", color = "red", size = 0.8, alpha= 0.5) + 
  geom_hline(yintercept=15.5, linetype="longdash", color = "red", size = 0.8, alpha= 0.5) + 
  
  RotatedAxis() #saving - wt= 953, ht = 532




############ T cells - proportions (Fig 4C) ########################


Idents(T.clust3) <- T.clust3$orig.ident
split.T.clust<- SplitObject(T.clust3 , split.by = "orig.ident")
split.T.wk8 <- split.T.clust[1:4]
split.T.wk12 <- split.T.clust[5:8]


##### Wk 8 T #####
T.prop.wk8 <- c()
for (i in 1:length(split.T.wk8) ){
  T.prop.wk8<- cbind(T.prop.wk8, as.vector(summary(split.T.wk8[[i]]$annotv1)))
}

T.prop.wk8 <- as.data.frame(T.prop.wk8)
rownames(T.prop.wk8) <- names(summary(T.clust3$annotv1))
colnames(T.prop.wk8) <- names(summary(T.clust3$orig.ident))[1:4]

T.prop.wk8[T.prop.wk8 == 0] <-NA
T.prop.wk8 <- na.omit(T.prop.wk8)

immune.total <- as.data.frame(summary(T.clust3$orig.ident))

T.prop.wk8.pct<-c()
for (i in 1:ncol(T.prop.wk8)){
  T.prop.wk8.pct <- cbind(T.prop.wk8.pct, T.prop.wk8[,i] / immune.total[,1][i])
}

T.prop.wk8.pct <- as.data.frame(T.prop.wk8.pct)
rownames(T.prop.wk8.pct) <- rownames(T.prop.wk8)
colnames(T.prop.wk8.pct) <- names(summary(T.clust3$orig.ident))[1:4]

#sanity check 
colSums(T.prop.wk8.pct) #good.

T.prop.wk8.pct$Celltype <- rownames(T.prop.wk8.pct)

T.prop.wk8.pct.gg<-T.prop.wk8.pct%>%tidyr::gather(Treatment,Count,1:4)
T.prop.wk8.pct.gg$Treatment<-factor(T.prop.wk8.pct.gg$Treatment, 
                                    levels= c('w8_ctrl','w8_lidpad','w8_ab','w8_lysm'))

#T.col<-tcol[which(names(mytcol) %in%  rownames(T.prop.wk8.pct))  ]
T.prop.wk8.pct.gg$Celltype<-factor(rownames(T.prop.wk8.pct), levels = names(mytcol) )

T.prop.wk8.pct.ggplot<-ggplot(T.prop.wk8.pct.gg,aes(fill=Celltype,y=Count,x=Treatment))+
  geom_bar(position=position_fill(reverse = F), stat = 'identity')+
  guides(fill=guide_legend(reverse=F))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Treatment")+
  ylab("Abundance (%)") +  scale_fill_manual(values=c(mytcol))

T.prop.wk8.pct.ggplot + ggtitle("Week 8")

##### WK 12 T#####

T.prop.wk12 <- c()
for (i in 1:length(split.T.wk12) ){
  T.prop.wk12<- cbind(T.prop.wk12, as.vector(summary(split.T.wk12[[i]]$annotv1)))
}

T.prop.wk12 <- as.data.frame(T.prop.wk12)
rownames(T.prop.wk12) <- names(summary(T.clust3$annotv1))
colnames(T.prop.wk12) <- names(summary(T.clust3$orig.ident))[5:8]

T.prop.wk12[T.prop.wk12 == 0] <-NA
T.prop.wk12 <- na.omit(T.prop.wk12)

immune.total <- as.data.frame(summary(T.clust3$orig.ident))

T.prop.wk12.pct<-c()
for (i in 1:ncol(T.prop.wk12)){
  T.prop.wk12.pct <- cbind(T.prop.wk12.pct, T.prop.wk12[,i] / immune.total[,1][i+4])
}

T.prop.wk12.pct <- as.data.frame(T.prop.wk12.pct)
rownames(T.prop.wk12.pct) <- rownames(T.prop.wk12)
colnames(T.prop.wk12.pct) <- names(summary(T.clust3$orig.ident))[5:8]

#sanity check 
colSums(T.prop.wk12.pct) #good.

T.prop.wk12.pct$Celltype <- rownames(T.prop.wk12.pct)

T.prop.wk12.pct.gg<-T.prop.wk12.pct%>%tidyr::gather(Treatment,Count,1:4)
T.prop.wk12.pct.gg$Treatment<-factor(T.prop.wk12.pct.gg$Treatment, 
                                     levels= c('w12_ctrl','w12_lidpad','w12_ab','w12_lysm'))

mytcol<-mytcol[which(names(mytcol) %in%  rownames(T.prop.wk12.pct))  ]
T.prop.wk12.pct.gg$Celltype<-factor(rownames(T.prop.wk12.pct), levels = names(mytcol) )

T.prop.wk12.pct.ggplot<-ggplot(T.prop.wk12.pct.gg,aes(fill=Celltype,y=Count,x=Treatment))+
  geom_bar(position=position_fill(reverse = F), stat = 'identity')+
  guides(fill=guide_legend(reverse=F))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Treatment")+
  ylab("Abundance (%)") +  scale_fill_manual(values=c(mytcol))

T.prop.wk12.pct.ggplot + ggtitle("Week 12")

T.prop.wk12.pct.ggplot + ggtitle("Week 12")



Idents(T.clust3) <- factor(T.clust3$annotv1)