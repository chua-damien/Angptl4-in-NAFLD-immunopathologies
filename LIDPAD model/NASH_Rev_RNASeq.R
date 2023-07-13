library(GEOquery)
library(DESeq2)
library(gtools)
library(readr)
setwd("C:/Users/Damien/OneDrive - Nanyang Technological University (1)/_PhD/_Projects/_LIDPAD MODEL/NASH_Rev/")


#### Loading datasets #### 

temp_Rev <- read_csv("./Datasets/reversion.csv")
colnames(temp_Rev) <- c("gene.ens", "NAFL_R_1", "NAFL_R_2", "NASH_R_1", "NASH_R_2",
                        "LP_w8_r1", "LP_w8_r2", "Ctrl_w12_r1")


###### ONLY Ctrl_w12, LP_w8, NASH_R

#run together. once. 
LP_comb <- temp_Rev
LP_R <- LP_comb[, c("gene.ens", 
                       "Ctrl_w12_r1", "LP_w8_r1", "LP_w8_r2",
                       "NASH_R_1", "NASH_R_2")]
#arrange dataset accding to conditions

LP_R1<- LP_R
rownames(LP_R1) <- LP_R1[,1]
LP_R1 <- LP_R1[,2:ncol(LP_R1)]
LP_comb1_conds <- data.frame(colnames(LP_R1),
                             conditions = factor(c(rep("w12_ctrl", 1), rep("w8_LP",2),
                                                   rep("NASH_R", 2)),
                                                 levels = c("w12_ctrl", "w8_LP", "NASH_R")) )

### DESeQ object  #### 

library(biobroom)
library(dplyr)
library(ggplot2)

dds_LP_R <- DESeqDataSetFromMatrix(countData = as.matrix(LP_R1), 
                                      colData = as.data.frame(LP_comb1_conds),
                                      design = ~ conditions) 
dds_LP_R <- DESeq(dds_LP_R)
res.LP_R <-results(dds_LP_R)
t.res.LP_R <- tidy.DESeqResults(res.LP_R)
t.res.LP_R <- arrange(t.res.LP_R, p.adjusted)

vds_LP_R <- vst(dds_LP_R, blind = F)
vsn::meanSdPlot(assay(vds_LP_R))


res.LP_R_filt<- results(dds_LP_R,alpha = 0.05, lfcThreshold = log2(1.5), contrast = c("conditions","w8_LP", "w12_ctrl"))
summary(res.LP_R_filt)
sum(res.LP_R_filt$padj < 0.05, na.rm=TRUE)
DEG_R_LPw8vCtrlw12<- res.LP_R_filt@rownames[res.LP_R_filt$padj < 0.05 & !is.na(res.LP_R_filt$padj)]

res.LP_R_filt<- results(dds_LP_R,alpha = 0.05, lfcThreshold = log2(1.5), contrast = c("conditions","NASH_R", "w12_ctrl"))
summary(res.LP_R_filt)
sum(res.LP_R_filt$padj < 0.05, na.rm=TRUE)
DEG_R_NASH_RvCtrlw12<- res.LP_R_filt@rownames[res.LP_R_filt$padj < 0.05 & !is.na(res.LP_R_filt$padj)]

res.LP_R_filt<- results(dds_LP_R,alpha = 0.05, lfcThreshold = log2(1.5), contrast = c("conditions","NASH_R", "w8_LP"))
summary(res.LP_R_filt)
sum(res.LP_R_filt$padj < 0.05, na.rm=TRUE)
DEG_R_NASH_RvLPw8<- res.LP_R_filt@rownames[res.LP_R_filt$padj < 0.05 & !is.na(res.LP_R_filt$padj)]

DEG_all <- union(DEG_R_LPw8vCtrlw12, DEG_R_NASH_RvCtrlw12) #309
DEG_all <- union(DEG_all, DEG_R_NASH_RvLPw8) #422


df_col <- data.frame(cluster= colData(dds_LP_R)[,c("conditions")])
rownames(df_col) <- colData(dds_LP_R)[,c("colnames.LP_R1.")]

annot_color <- c(rep("black",1), rep("red",2), rep("grey", 2))
names(annot_color) <- df_col$cluster
annot_color <- list(cluster = annot_color)

p_R<- pheatmap(assay(vds_LP_R)[DEG_all,c("Ctrl_w12_r1", "LP_w8_r1", "LP_w8_r2",
                                                   "NASH_R_1", "NASH_R_2") ], 
               scale = "row", cluster_rows=T, show_rownames=FALSE, 
               cluster_cols=F, show_colnames = F, 
               annotation_col = df_col,
               annotation_colors = annot_color,
               main = "DEG_all") 


  
#### Gene enrichment ####

p_R<- pheatmap(assay(vds_LP_R)[DEG_all,c("Ctrl_w12_r1", "LP_w8_r1", "LP_w8_r2",
                                         "NASH_R_1", "NASH_R_2") ], 
               scale = "row", cluster_rows=T, show_rownames=FALSE, 
               cluster_cols=F, show_colnames = T, 
               annotation_col = df_col,
               annotation_colors = annot_color,
               cutree_row = 4,
               main = "DEG_all") 


#extract terms 
cl_NASH_R <- cutree(p_R$tree_row, 4)

library(clusterProfiler)
library(AnnotationHub)
library(dplyr)
library(annotables)
library(org.Mm.eg.db)

NASH_R_1 <- data.frame(gene = names(cl_NASH_R)[cl_NASH_R == 1])
NASH_R_1 <- inner_join (NASH_R_1, grcm38, by = c("gene" = "ensgene"))
NASH_R_2 <- data.frame(gene = names(cl_NASH_R)[cl_NASH_R == 2])
NASH_R_3 <- inner_join (NASH_R_2, grcm38, by = c("gene" = "ensgene"))
NASH_R_3 <- data.frame(gene = names(cl_NASH_R)[cl_NASH_R == 3])
NASH_R_4 <- inner_join (NASH_R_3, grcm38, by = c("gene" = "ensgene"))
NASH_R_4 <- data.frame(gene = names(cl_NASH_R)[cl_NASH_R == 4])
NASH_R_4 <- inner_join (NASH_R_4, grcm38, by = c("gene" = "ensgene"))

go_NASH_R_1 <-enrichGO(names(cl_NASH_R)[cl_NASH_R == 1],
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         keyType = "ENSEMBL",
                         readable      = TRUE)

go_NASH_R_2 <-enrichGO(names(cl_NASH_R)[cl_NASH_R == 2],
                           OrgDb         = org.Mm.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           keyType       = 'ENSEMBL',
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)

go_NASH_R_3 <-enrichGO(names(cl_NASH_R)[cl_NASH_R == 3],
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       keyType       = 'ENSEMBL',
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

go_NASH_R_4 <-enrichGO(names(cl_NASH_R)[cl_NASH_R == 4],
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       keyType       = 'ENSEMBL',
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

library(DOSE)

res_LP_R <- results(dds_LP_R, contrast = c("conditions", "w8_LP", "w12_ctrl"))
t.res.LP_R <- tidy.DESeqResults(res_LP_R)
t.res.LP_R <- arrange(t.res.LP_R, p.adjusted)
t.res.LP_R <- inner_join (t.res.LP_R, grcm38, by = c("gene" = "ensgene"))

LPw8_genelist <- data.frame(gene = t.res.LP_R$gene,
                            value = t.res.LP_R$estimate)
LPw8_genelist <- LPw8_genelist[which(LPw8_genelist$gene %in% DEG_all),]
LPw8_genelist <- LPw8_genelist[order(LPw8_genelist$value, decreasing = T),]
LPw8_genelist1 <- LPw8_genelist$value
names(LPw8_genelist1) <- LPw8_genelist$gene

res_NASH_R <- results(dds_LP_R, contrast = c("conditions", "NASH_R", "w12_ctrl"))
t.res.NASH_R <- tidy.DESeqResults(res_NASH_R)
t.res.NASH_R <- arrange(t.res.NASH_R, p.adjusted)
t.res.NASH_R <- inner_join (t.res.NASH_R, grcm38, by = c("gene" = "ensgene"))

NASH_R_genelist <- data.frame(gene = t.res.NASH_R$gene,
                              value = t.res.NASH_R$estimate)
NASH_R_genelist <- NASH_R_genelist[which(NASH_R_genelist$gene %in% DEG_all),]
NASH_R_genelist <- NASH_R_genelist[order(NASH_R_genelist$value, decreasing = T),]
NASH_R_genelist1 <- NASH_R_genelist$value
names(NASH_R_genelist1) <- NASH_R_genelist$gene


#### GSEA - whole gene list not working. #### 


library(DOSE)

res_LP_R <- results(dds_LP_R, contrast = c("conditions", "w8_LP", "w12_ctrl"))
t.res.LP_R <- tidy.DESeqResults(res_LP_R)
t.res.LP_R <- arrange(t.res.LP_R, p.adjusted)
t.res.LP_R <- inner_join (t.res.LP_R, grcm38, by = c("gene" = "ensgene"))


LPw8_genelist2 <- data.frame(gene = t.res.LP_R$symbol,
                               value = t.res.LP_R$estimate)
LPw8_genelist2 <- na.omit(LPw8_genelist2)
LPw8_genelist2 <- LPw8_genelist2[order(LPw8_genelist2$value, decreasing = T),]
LPw8_genelist2 <- LPw8_genelist2[!duplicated(LPw8_genelist2),]
LPw8_genelist3 <- LPw8_genelist2$value
names(LPw8_genelist3) <- LPw8_genelist2$gene
write.csv(LPw8_genelist2, "LPw8_genelist2.csv")



res_NASH_R <- results(dds_LP_R, contrast = c("conditions", "NASH_R", "w12_ctrl"))
t.res.NASH_R <- tidy.DESeqResults(res_NASH_R)
t.res.NASH_R <- arrange(t.res.NASH_R, p.adjusted)
t.res.NASH_R <- inner_join (t.res.NASH_R, grcm38, by = c("gene" = "ensgene"))

#for export
NASH_R_genelist2 <- data.frame(gene = t.res.NASH_R$symbol,
                              value = t.res.NASH_R$estimate)
NASH_R_genelist2 <- na.omit(NASH_R_genelist2)
NASH_R_genelist2 <- NASH_R_genelist2[order(NASH_R_genelist2$value, decreasing = T),]
NASH_R_genelist2 <- NASH_R_genelist2[!duplicated(NASH_R_genelist2),]
NASH_R_genelist3 <- NASH_R_genelist2$value
names(NASH_R_genelist3) <- NASH_R_genelist2$gene
write.csv(NASH_R_genelist2, "NASH_R_genelist2.csv")


#### IMPORT FROM EXTERNAL APPLICATION ####


LPw8_gsea_e_res <- read.csv("./GSEA/Output/LPw8_gsea_results.csv")
rownames(LPw8_gsea_e_res) <- LPw8_gsea_e_res$Geneset


NASH_R_gsea_e_res <- read.csv("./GSEA/Output/NASH_R_gsea_results.csv")
rownames(NASH_R_gsea_e_res) <- NASH_R_gsea_e_res$Geneset

#### Plotting ####          

LPw8_gsea_e_res_filt <- LPw8_gsea_e_res[c("ANTIMICROBIAL HUMORAL RESPONSE",
                                  "DEFENSE RESPONSE TO BACTERIUM",
                                  "CELLULAR COMPONENT ASSEMBLY INVOLVED IN MORPHOGENESIS",
                                  "ACUTE PHASE RESPONSE",
                                  "TISSUE REGENERATION",
                                  "REGENERATION",
                                  "NEGATIVE REGULATION OF TISSUE REMODELING",
                                  "COLLAGEN CATABOLIC PROCESS",
                                  "COLLAGEN FIBRIL ORGANIZATION",
                                  "EXTRACELLULAR MATRIX ASSEMBLY",
                                  "CELL SURFACE TOLL LIKE RECEPTOR SIGNALING PATHWAY",
                                  "ALANINE TRANSPORT",
                                  "ORGANIC HYDROXY COMPOUND CATABOLIC PROCESS",
                                  "ACUTE INFLAMMATORY RESPONSE",
                                  "REGULATION OF EXTRACELLULAR MATRIX ASSEMBLY",
                                  "CELLULAR RESPONSE TO OXYGEN RADICAL",
                                  "POSITIVE REGULATION OF TUMOR NECROSIS FACTOR SUPERFAMILY CYTOKINE PRODUCTION",
                                  "SEQUESTERING OF TRIGLYCERIDE",
                                  "RESPONSE TO TYPE II INTERFERON",
                                  "POSITIVE REGULATION OF PATTERN RECOGNITION RECEPTOR SIGNALING PATHWAY",
                                  "POSITIVE REGULATION OF LIPID KINASE ACTIVITY",
                                  "TUMOR NECROSIS FACTOR SUPERFAMILY CYTOKINE PRODUCTION",
                                  "STEROID CATABOLIC PROCESS",
                                  "POSITIVE REGULATION OF LIPID STORAGE",
                                  "NON CANONICAL WNT SIGNALING PATHWAY",
                                  "INOSITOL PHOSPHATE BIOSYNTHETIC PROCESS",
                                  "POSITIVE REGULATION OF ENDOTHELIAL CELL APOPTOTIC PROCESS",
                                  "NEGATIVE REGULATION OF ENDOTHELIAL CELL PROLIFERATION",
                                  "POSITIVE REGULATION OF ENDOTHELIAL CELL PROLIFERATION",
                                  "VENOUS BLOOD VESSEL DEVELOPMENT",
                                  "CHOLESTEROL STORAGE",
                                  "LIPID STORAGE",
                                  "RESPONSE TO LEPTIN",
                                  
                                  "PLATELET MORPHOGENESIS",
                                  "REGULATION OF PLATELET ACTIVATION",
                                  "POSITIVE REGULATION OF PROTEIN DEPOLYMERIZATION",
                                  "LIPOPROTEIN METABOLIC PROCESS",
                                  "INDOLE CONTAINING COMPOUND METABOLIC PROCESS"
                                  
                                  
                                  ), c("Geneset","NES", "FDR")]



LPw8_gsea_e_res_filt<- LPw8_gsea_e_res_filt[order(LPw8_gsea_e_res_filt$NES, decreasing = T),]
LPw8_gsea_e_res_filt$condition <- "LP_w8"
LPw8_gsea_e_res_filt1<- LPw8_gsea_e_res_filt[order(LPw8_gsea_e_res_filt$NES, decreasing = F),]




NASH_R_gsea_e_res_filt <- NASH_R_gsea_e_res[c("ANTIMICROBIAL HUMORAL RESPONSE",
                                  "DEFENSE RESPONSE TO BACTERIUM",
                                  "CELLULAR COMPONENT ASSEMBLY INVOLVED IN MORPHOGENESIS",
                                  "ACUTE PHASE RESPONSE",
                                  "TISSUE REGENERATION",
                                  "REGENERATION",
                                  "NEGATIVE REGULATION OF TISSUE REMODELING",
                                  "COLLAGEN CATABOLIC PROCESS",
                                  "COLLAGEN FIBRIL ORGANIZATION",
                                  "EXTRACELLULAR MATRIX ASSEMBLY",
                                  "CELL SURFACE TOLL LIKE RECEPTOR SIGNALING PATHWAY",
                                  "ALANINE TRANSPORT",
                                  "ORGANIC HYDROXY COMPOUND CATABOLIC PROCESS",
                                  "ACUTE INFLAMMATORY RESPONSE",
                                  "REGULATION OF EXTRACELLULAR MATRIX ASSEMBLY",
                                  "CELLULAR RESPONSE TO OXYGEN RADICAL",
                                  "POSITIVE REGULATION OF TUMOR NECROSIS FACTOR SUPERFAMILY CYTOKINE PRODUCTION",
                                  "SEQUESTERING OF TRIGLYCERIDE",
                                  "RESPONSE TO TYPE II INTERFERON",
                                  "POSITIVE REGULATION OF PATTERN RECOGNITION RECEPTOR SIGNALING PATHWAY",
                                  "POSITIVE REGULATION OF LIPID KINASE ACTIVITY",
                                  "TUMOR NECROSIS FACTOR SUPERFAMILY CYTOKINE PRODUCTION",
                                  "STEROID CATABOLIC PROCESS",
                                  "POSITIVE REGULATION OF LIPID STORAGE",
                                  "NON CANONICAL WNT SIGNALING PATHWAY",
                                  "INOSITOL PHOSPHATE BIOSYNTHETIC PROCESS",
                                  "POSITIVE REGULATION OF ENDOTHELIAL CELL APOPTOTIC PROCESS",
                                  "NEGATIVE REGULATION OF ENDOTHELIAL CELL PROLIFERATION",
                                  "POSITIVE REGULATION OF ENDOTHELIAL CELL PROLIFERATION",
                                  "VENOUS BLOOD VESSEL DEVELOPMENT",
                                  "CHOLESTEROL STORAGE",
                                  "LIPID STORAGE",
                                  "RESPONSE TO LEPTIN",
                                  
                                  "PLATELET MORPHOGENESIS",
                                  "REGULATION OF PLATELET ACTIVATION",
                                  "POSITIVE REGULATION OF PROTEIN DEPOLYMERIZATION",
                                  "LIPOPROTEIN METABOLIC PROCESS",
                                  "INDOLE CONTAINING COMPOUND METABOLIC PROCESS"
                                  
                                  
                                  ), c("Geneset","NES", "FDR")]


NASH_R_gsea_e_res_filt<- NASH_R_gsea_e_res_filt[order(NASH_R_gsea_e_res_filt$NES, decreasing = T),]
NASH_R_gsea_e_res_filt$condition <- "NASH_R"

gsea_dotplot <- rbind(LPw8_gsea_e_res_filt, NASH_R_gsea_e_res_filt)
gsea_dotplot$FDR[which(gsea_dotplot$FDR ==0 )] <- 0.0004 
gsea_dotplot$log10p <- -log10(gsea_dotplot$FDR)
gsea_dotplot$condition <- factor(gsea_dotplot$condition, levels = c("LP_w8", "NASH_R"))
gsea_dotplot$Geneset <- factor(gsea_dotplot$Geneset, levels = LPw8_gsea_e_res_filt1$Geneset)

library(ggthemes)
ggplot(gsea_dotplot, (aes(x=condition, y=Geneset, color = as.numeric(NES), size=as.numeric(log10p)))) + 
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red")    +
  labs(colour = "NES", size = "-log10(FDR)") +
  theme_bw()



#### NASH_R vs LP_w8 #### 

library(DOSE)
res_LP_R <- results(dds_LP_R, contrast = c("conditions", "NASH_R", "w8_LP"))
t.res.LP_R <- tidy.DESeqResults(res_LP_R)
t.res.LP_R <- arrange(t.res.LP_R, p.adjusted)
t.res.LP_R <- inner_join (t.res.LP_R, grcm38, by = c("gene" = "ensgene"))

NASH_RvLPw8_genelist <- data.frame(gene = t.res.LP_R$gene,
                            value = t.res.LP_R$estimate)
NASH_RvLPw8_genelist <- na.omit(NASH_RvLPw8_genelist)
NASH_RvLPw8_genelist <- NASH_RvLPw8_genelist[order(NASH_RvLPw8_genelist$value, decreasing = T),]
NASH_RvLPw8_genelist <- NASH_RvLPw8_genelist[!duplicated(NASH_RvLPw8_genelist),]
NASH_RvLPw8_genelist1 <- NASH_RvLPw8_genelist$value
names(NASH_RvLPw8_genelist1) <- NASH_RvLPw8_genelist$gene

#### Enrichment Plots #### 

NASH_RvLPw8_GSEdo <- gseGO(NASH_RvLPw8_genelist1,
                           ont = "ALL",
                           OrgDb = org.Mm.eg.db,
                           keyType       = 'ENSEMBL',
                           minGSSize = 8,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "fdr",
                           seed = 8888)


gseaplot2(NASH_RvLPw8_GSEdo, geneSetID = c(1, #acute inflam
                                           216, #
                                           35, 
                                           386, 
                                           222), pvalue_table = T, ES_geom = "line")
