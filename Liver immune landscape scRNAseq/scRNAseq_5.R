


Idents(split.T.act.clust$CD8) <- split.T.act.clust$CD8$orig.ident 
DefaultAssay(split.T.act.clust$CD8) <- "RNA"

CD8T_w8_all_ab_lp <- FindMarkers(split.T.act.clust$CD8, ident.1 = "w8_ab", ident.2 = "w8_lidpad", test.use = "DESeq2", logfc.threshold = 0)
CD8T_w8_all_lp_ctrl <- FindMarkers(split.T.act.clust$CD8, ident.1 = "w8_lidpad", ident.2 = "w8_ctrl", test.use = "DESeq2", logfc.threshold = 0)
CD8T_w12_all_lp_ctrl <- FindMarkers(split.T.act.clust$CD8, ident.1 = "w12_lidpad", ident.2 = "w12_ctrl", test.use = "DESeq2", logfc.threshold = 0)
CD8T_w12_all_ab_lp <- FindMarkers(split.T.act.clust$CD8, ident.1 = "w12_ab", ident.2 = "w12_lidpad", test.use = "DESeq2", logfc.threshold = 0)


write.csv(data.frame("Gene"= rownames(CD8T_w8_all_ab_lp ),
                     "avg_log2FC" = CD8T_w8_all_ab_lp $avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD8T_w8_all_ab_lp.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD8T_w8_all_lp_ctrl),
                     "avg_log2FC" = CD8T_w8_all_lp_ctrl$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD8T_w8_all_lp_ctrl.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD8T_w12_all_lp_ctrl),
                     "avg_log2FC" = CD8T_w12_all_lp_ctrl$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD8_w12_lidpad_ctrl.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD8T_w12_all_ab_lp),
                     "avg_log2FC" = CD8T_w12_all_ab_lp$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD8T_w12_all_ab_lp.csv", row.names=F)



Idents(split.T.act.clust$CD4) <- split.T.act.clust$CD4$orig.ident 
DefaultAssay(split.T.act.clust$CD4) <- "RNA"

CD4T_w8_all_ab_lp <- FindMarkers(split.T.act.clust$CD4, ident.1 = "w8_ab", ident.2 = "w8_lidpad", test.use = "DESeq2", logfc.threshold = 0)
CD4T_w8_all_lp_ctrl <- FindMarkers(split.T.act.clust$CD4, ident.1 = "w8_lidpad", ident.2 = "w8_ctrl", test.use = "DESeq2", logfc.threshold = 0)
CD4T_w12_all_lp_ctrl <- FindMarkers(split.T.act.clust$CD4, ident.1 = "w12_lidpad", ident.2 = "w12_ctrl", test.use = "DESeq2", logfc.threshold = 0)
CD4T_w12_all_ab_lp <- FindMarkers(split.T.act.clust$CD4, ident.1 = "w12_ab", ident.2 = "w12_lidpad", test.use = "DESeq2", logfc.threshold = 0)

write.csv(data.frame("Gene"= rownames(CD4T_w8_all_ab_lp ),
                     "avg_log2FC" = CD4T_w8_all_ab_lp $avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD4T_w8_all_ab_lp.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD4T_w8_all_lp_ctrl),
                     "avg_log2FC" = CD4T_w8_all_lp_ctrl$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD4T_w8_all_lp_ctrl.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD4T_w12_all_lp_ctrl),
                     "avg_log2FC" = CD4T_w12_all_lp_ctrl$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD4_w12_lidpad_ctrl.csv", row.names=F)
write.csv(data.frame("Gene"= rownames(CD4T_w12_all_ab_lp),
                     "avg_log2FC" = CD4T_w12_all_ab_lp$avg_log2FC), "C:/Users/LKC/Desktop/Damien/scNASH/Output/log2fc_genelist/CD4T_w12_all_ab_lp.csv", row.names=F)
