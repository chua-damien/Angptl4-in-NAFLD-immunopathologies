### T cell subtypes pseudotime - Fig 5B, S9B ####


# Pseudotime densities (by spatial)
ds <- list( w8_ctrl = density(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ctrl", 1], na.rm=T),
            w8_lidpad = density(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lidpad", 1], na.rm=T),
            w8_ab = density(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ab", 1], na.rm=T),
            w8_lysm = density(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lysm", 1], na.rm=T))
xlim <- range(c(ds$w8_ctrl$x, ds$w8_lidpad$x, ds$w8_ab$x, ds$w8_lysm$x))
ylim <- range(c(ds$w8_ctrl$y, ds$w8_lidpad$y, ds$w8_ab$y, ds$w8_lysm$y))
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")
polygon(c(min(ds$w8_ctrl$x), ds$w8_ctrl$x, max(ds$w8_ctrl$x)), c(0, ds$w8_ctrl$y, 0),
        col = alpha(brewer.pal(4, "Set1")[1], alpha = .5))
polygon(c(min(ds$w8_lidpad$x), ds$w8_ctrl$x, max(ds$w8_ctrl$x)), c(0, ds$w8_lidpad$y, 0),
        col = alpha(brewer.pal(9, "Set1")[2], alpha = .5))
polygon(c(min(ds$w8_ab$x), ds$w8_ctrl$x, max(ds$w8_ctrl$x)), c(0, ds$w8_ab$y, 0),
        col = alpha(brewer.pal(9, "Set1")[5], alpha = .5))
polygon(c(min(ds$w8_lysm$x), ds$w8_ctrl$x, max(ds$w8_ctrl$x)), c(0, ds$w8_lysm$y, 0),
        col = alpha(brewer.pal(9, "Set1")[6], alpha = .5))
legend("topleft", legend = c("w8_ctrl", "w8_lidpad", "w8_ab", "w8_lysm"), 
       fill = alpha(brewer.pal(9, "Set1")[c(1,2,5,6)], alpha = .5), bty = "n")



# creating ridgeplots for differential progression
CD4_pseudot_list <- list()
for(i in 1:ncol(slingPseudotime(CD4T.clust.ss))){
  df_pseudot <- data.frame( treatment = c(
    rep("w8_ctrl",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ctrl", i])),
    rep("w8_lidpad",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lidpad", i])),
    rep("w8_ab",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ab", i])),
    rep("w8_lysm",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lysm", i])),
    rep("w12_ctrl",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_ctrl", i])),
    rep("w12_lidpad",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_lidpad", i])),
    rep("w12_ab",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_ab", i])),
    rep("w12_lysm",length(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_lysm", i]))),
    
    pseudotime = c(slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ctrl", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lidpad", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_ab", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w8_lysm", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_ctrl", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_lidpad", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_ab", i],
                   slingPseudotime(CD4T.clust.ss)[colData(CD4T.clust.ss)$orig.ident == "w12_lysm", i]
    ) )
  
  df_pseudot <- na.omit(df_pseudot)
  
  df_pseudot$treatment <- factor(df_pseudot$treatment, levels = c("w8_ctrl", "w8_lidpad", "w8_ab", "w8_lysm",
                                                                  "w12_ctrl", "w12_lidpad", "w12_ab", "w12_lysm"))
  
  
  CD4_pseudot_list[[paste0("curve_", i)]] <- df_pseudot
}

SlingshotDataSet(CD4T.clust.ss)
#Curve 1 - Treg_c27
#Curve 2 - Th17
#curve 3 - tnf+ CD4 CTL 

names(CD4_pseudot_list)<- c("lin1_Tregc27",
                            "lin2_Th17",
                            "lin3_Th1")

#plotting ridgeline 
library(ridgeline)
for(i in 1:length(CD4_pseudot_list)){
  print(i)
  ridgeline(CD4_pseudot_list[[i]]$pseudotime, CD4_pseudot_list[[i]]$treatment)
  title(main = paste0(names(CD4_pseudot_list)[i]), xlab = "pseudotime")
  
}

library(ggridges)

for(i in 1:length(CD4_pseudot_list)){
  CD4_pseudot_list[[i]]$treatment <- forcats::fct_rev(CD4_pseudot_list[[i]]$treatment)
}


col <- 
  #Fig3 - w12, w8 separated, pseudotime curves overlaid
  for (i in c("w8", "w1")){
    print(i)
    
    for(j in 1:length(CD4_pseudot_list)){
      print(j)
      print(ggplot(CD4_pseudot_list[[j]][which(substr(CD4_pseudot_list[[j]]$treatment,1,2 ) == i),], aes(x=pseudotime, y=treatment, fill = stat(x))) +
              geom_density_ridges_gradient(scale = 1.5, rel_min_heght = 0.02, gradient_lwd = 1.) + 
              #stat_density_ridges(quantile_lines = T, quantiles =2) +
              ggtitle(paste0(names(CD4_pseudot_list)[j])) +
              scale_x_continuous(limits = c(-2,13), expand = c(0, 0)) +
              scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.25))) +
              scale_fill_viridis_c(name = "pseudotime", option = "C") +
              theme_ridges(font_size = 13, grid = TRUE) +
              theme(axis.title.y = element_blank())
      )
    }
  }

#Fig3 - w12, w8 together, pseudotime curves super overlaid

for(j in 1:length(CD4_pseudot_list)){
  print(j)
  print(ggplot(CD4_pseudot_list[[j]], aes(x=pseudotime, y=treatment, fill = treatment)) +
          geom_density_ridges(scale = 6, rel_min_heght = 0.03, size= 0.25, gradient_lwd = 1.,
                              alpha = 0.5, color= "black") + 
          #stat_density_ridges(quantile_lines = T, quantiles =2) +
          ggtitle(paste0(names(CD4_pseudot_list)[j])) +
          scale_x_continuous(limits = c(-2,13), expand = c(0, 0)) +
          #scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.25))) +
          scale_fill_cyclical(breaks =  c("w8_ctrl", "w8_lidpad", "w8_ab", "w8_lysm",
                                          "w12_ctrl", "w12_lidpad", "w12_ab", "w12_lysm"),
                              values = c("#EEC900", "#886300", "#FFFC33", "#FFFF99",
                                         "#CD3333", "#670000", "#AA5555", "#FFBBBB"),
                              name = "Option", guide = "legend") +
          
          theme_ridges(font_size = 13, grid = TRUE) +
          theme(axis.title.y = element_blank())
  )
}
geom_vline (xintercept=8.2, linetype="dotted", color = "black", size = 1 )+
  geom_vline (xintercept=9.8, linetype="dashed", color = "black", size = 1 )
#lin1: 
#lin2: 
#lin3: ctrl - 8.2, lidpad - 9.8



library(ggridges)
ggplot(CD4_pseudot_list$curve_1_Tregc27, aes(x = pseudotime, y = treatment)) + geom_density_ridges() +
  stat_density_ridges(quantile_lines = T, quantiles =2)+ 
  theme_ridges()

#ridgeline(df_pseudot$pseudotime, df_pseudot$treatment) +
#   title(main = "CD4 T cells", xlab = "pseudotime")
