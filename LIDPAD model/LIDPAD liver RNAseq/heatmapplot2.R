library(gplots)

dat1 = read.csv("masteresti310119.csv")


#make row names

dat2 <- as.matrix(dat1)
rownames(dat1) <- dat2[,1]

#remove column 1,2
dat3 <- dat1[-c(1:2)]

#make col names
cname1 <- c("Week1", "Week4", "Week8", "Week12", "Week16", "week32", "week48")
colnames(dat3) <- cname1
dat3[is.na(dat3)] <- 0

colors = c(seq(-5,-1, length = 25), seq(-0.99,-0.001, length = 25), seq(0, length=1),seq(0.001, 0.99, length = 25), seq(1, 5, length=25))
my_palette <- c(colorRampPalette(c("green","black"))(n=49),"#f8f8f8", colorRampPalette(c("black", "red"))(n=50))
plot(hclust(dist(dat3)))

#no dendrogram (according to excel order)
heatmap.2(as.matrix(dat3[dat1,cname1]), Colv = FALSE, Rowv = FALSE, trace = "none", margin=c(8,9), breaks = colors, col = my_palette)

#dendrogram
heatmap.2(as.matrix(dat3[,cname1]), Colv = FALSE, dendrogram = "row", trace = "none", margin=c(8,9), breaks = colors, col = my_palette)

#both dendrograms
heatmap.2(as.matrix(dat3[,cname1]), dendrogram = "both", trace = "none", margin=c(8,9), breaks = colors, col = my_palette)

#dendrogram color
heatmap.2(as.matrix(dat3[,cname1]), Colv = FALSE, dendrogram = "row", trace = "none", margin=c(8,9), breaks = colors, col = my_palette, RowSideColors = as.character(ct2))

#retrieve rownames
x <- heatmap.2(as.matrix(dat3[,cname1]), Colv = FALSE, dendrogram = "row", trace = "none", margin=c(8,9), breaks = colors, col = my_palette)
y <- rev(rownames(dat3)[x$rowInd])
write.table(y, "Geneclustfilt.txt", row.names = FALSE)

#dendrograms grouping
hc.rows <- hclust(dist(dat3))
ct <- cutree(hc.rows, h=10)
rect.hclust(hc.rows, h=10)
table(ct)
tableclust <-  data.frame(dat3,ct)
ct2 <- cutree(hc.rows, h=7)
rect.hclust(hc.rows, h=7)
tableclust2 <- data.frame(tableclust, ct2)
tableclust3 <- data.frame(rownames(dat3),tableclust2)
write.table(tableclust3, "dendrofull.txt", row.names = FALSE)
