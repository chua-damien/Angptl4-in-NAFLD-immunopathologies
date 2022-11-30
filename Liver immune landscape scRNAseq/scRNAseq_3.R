#### Activation status Fig 5a, S9a ####

# T cell immune activation
Idents(T.clust3) <- T.clust3$annotv1
T.clust3 <- RenameIdents(T.clust3,
                         "Memory CD8+ Tem" = "Activated CD8 T",
                         "Memory CD8 T_c46" = "Activated CD8 T",
                         "Memory CD8+ Trm" = "Activated CD8 T",
                         "Treg_c13" = "Activated CD4 T",
                         "Treg_c27" = "Activated CD4 T",
                         "Th1" = "Activated CD4 T",
                         "Th17" = "Activated CD4 T")
T.clust4 <- subset(T.clust3, idents = c("Naive CD4 T", "Activated CD4 T", "Naive CD8 T", "Activated CD8 T"))
T.clust4$T.act <- Idents(T.clust4)
T.clust4 <- RenameIdents(T.clust4, 
                         "Naive CD4 T" = "CD4", 
                         "Activated CD4 T" = "CD4", 
                         "Naive CD8 T" = "CD8", 
                         "Activated CD8 T" = "CD8")

T.clust4

T.clust4$cd4.8 <- Idents(T.clust4)
split.T.act.clust<- SplitObject(T.clust4, split.by = "cd4.8")
split.CD8.clust <- SplitObject(split.T.act.clust[[1]], split.by = "orig.ident")
split.CD4.clust <- SplitObject(split.T.act.clust[[2]], split.by = "orig.ident")




### CD8
CD8.act.prop <- c()
for (i in 1:length(split.CD8.clust) ){
  CD8.act.prop<- cbind(CD8.act.prop, as.vector(summary(split.CD8.clust[[i]]$T.act)))
}

CD8.act.total <- as.data.frame(summary(split.T.act.clust[[1]]$orig.ident))

CD8.act.prop.pct<-c()
for (i in 1:ncol(CD8.act.prop)){
  CD8.act.prop.pct <- cbind(CD8.act.prop.pct, CD8.act.prop[,i] / CD8.act.total[,1][i])
}

CD8.act.prop.pct <- as.data.frame(CD8.act.prop.pct)
rownames(CD8.act.prop.pct) <- names(summary(split.T.act.clust[[1]]$T.act))
colnames(CD8.act.prop.pct) <- names(summary(split.T.act.clust[[1]]$orig.ident))

CD8.act.prop.pct[CD8.act.prop.pct == 0] <-NA
CD8.act.prop.pct <- na.omit(CD8.act.prop.pct)
CD8.act.prop.pct

#sanity check 
colSums(CD8.act.prop.pct) #good.

CD8.act.prop.pct$Celltype <- rownames(CD8.act.prop.pct)


CD8.act.prop.pct.gg<-CD8.act.prop.pct%>%tidyr::gather(Treatment,Count,1:8)
CD8.act.prop.pct.gg$Treatment<-factor(CD8.act.prop.pct.gg$Treatment, 
                                      levels= c('w8_ctrl','w8_lidpad','w8_ab','w8_lysm',
                                                'w12_ctrl','w12_lidpad','w12_ab','w12_lysm'))
CD8.act.prop.pct.gg$Celltype<-factor(rownames(CD8.act.prop.pct))

col<-c("#0000FF","#B3B3FF")

CD8.act.prop.pct.ggplot<-ggplot(CD8.act.prop.pct.gg,aes(fill=Celltype,y=Count,x=Treatment))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity')+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Treatment groups")+
  ylab("Proportion of CD8+ T cells") +  scale_fill_manual(values=c(col))

#plot
CD8.act.prop.pct.ggplot 


### CD4
CD4.act.prop <- c()
for (i in 1:length(split.CD4.clust) ){
  CD4.act.prop<- cbind(CD4.act.prop, as.vector(summary(split.CD4.clust[[i]]$T.act)))
}

CD4.act.total <- as.data.frame(summary(split.T.act.clust[[2]]$orig.ident))

CD4.act.prop.pct<-c()
for (i in 1:ncol(CD4.act.prop)){
  CD4.act.prop.pct <- cbind(CD4.act.prop.pct, CD4.act.prop[,i] / CD4.act.total[,1][i])
}

CD4.act.prop.pct <- as.data.frame(CD4.act.prop.pct)
rownames(CD4.act.prop.pct) <- names(summary(split.T.act.clust[[1]]$T.act))
colnames(CD4.act.prop.pct) <- names(summary(split.T.act.clust[[1]]$orig.ident))

CD4.act.prop.pct[CD4.act.prop.pct == 0] <-NA
CD4.act.prop.pct <- na.omit(CD4.act.prop.pct)
CD4.act.prop.pct

CD4.act.prop.pct$Celltype <- rownames(CD4.act.prop.pct)


CD4.act.prop.pct.gg<-CD4.act.prop.pct%>%tidyr::gather(Treatment,Count,1:8)
CD4.act.prop.pct.gg$Treatment<-factor(CD4.act.prop.pct.gg$Treatment, 
                                      levels= c('w8_ctrl','w8_lidpad','w8_ab','w8_lysm',
                                                'w12_ctrl','w12_lidpad','w12_ab','w12_lysm'))
CD4.act.prop.pct.gg$Celltype<-factor(rownames(CD4.act.prop.pct))

col<-c("#FF0000","#FF7F7F")

CD4.act.prop.pct.ggplot<-ggplot(CD4.act.prop.pct.gg,aes(fill=Celltype,y=Count,x=Treatment))+
  geom_bar(position=position_fill(reverse    = T), stat = 'identity')+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(size=1,color="black"))+
  xlab("Treatment groups")+
  ylab("Proportion of CD4+ T cells (%)") +  scale_fill_manual(values=c(col))

#plot
CD4.act.prop.pct.ggplot 

  