library(ggplot2)
library(cowplot)

eig <- read.table("Documents/FIBR_temporal/pca/FIBR_thin_plink.eigenvec")
meta <- read.csv("Documents/FIBR_temporal/pca/meta_data.csv")

final <- data.frame(cbind(meta, eig))

final$POP <- as.factor(final$POP)
final$RIVER <- as.factor(final$RIVER)
final$YEAR <- as.factor((final$YEAR))
final_sort <- final[order(final$POP, final$RIVER),]

table(final_sort$RIVER)
table(final_sort$POP)
table(final_sort$YEAR)

sum_text <- data.frame(stream= c("CA", "TA", "LL", "UL"), x=c(0.08,0.15,-0.1,-0.015), y=c(-0.1,0.05,-0.07,0.05))
sum_text <- data.frame(stream= c("CA", "TA", "LL", "UL", "GH"), x=c(0.08,0.15,-0.1,-0.015,0.05), y=c(-0.1,0.05,-0.07,0.05,-0.01))


pca_plot <- ggplot(data=final_sort,aes(x=V3, y=V4)) +
  geom_point(aes(color=RIVER, shape = YEAR, size=2))+
  xlab("PC1 (4.47%)")+
  ylab("PC2 (3.24%)")+
 scale_colour_manual(values = c(GH = "#009E73", CA = "#FFB90F", TA = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
   theme_bw()+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 20)
  )+
  geom_text(data = sum_text,
            mapping = aes(x = x, y = y, label = stream), size=6,colour=c(CA="#FFB90F",TA="#FF7F00",LL="#5D478B",UL="#CD0000", GH="#009E73"))


pca_plot

pdf("Documents/FIBR_temporal/am_nat_revisions/pca/pca_revisions.pdf", height=5, width=7)
pca_plot
dev.off()

pca_plot_3_4 <- ggplot(data=final_sort,aes(x=V5, y=V6)) +
  geom_point(aes(color=RIVER, shape = YEAR))+
  xlab("PC3 (3.11%)")+
  ylab("PC4 (2.02%)")+
  scale_colour_manual(values = c(GH = "#009E73", CA = "#FFB90F", TA = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
  theme_bw()+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

pca_plot_3_4


eigval <- read.table("Documents/FIBR_temporal/pca/FIBR.thin.eigenval")
head(eigval)
eigval$number <- 1:10  

scree_plot <- ggplot(eigval, aes(x = factor(number), y = V1)) +
  geom_col() +
  labs(x = "Dimension", y = "Eigenvalue") +
  theme_bw()

supplement_plot <- plot_grid(pca_plot_3_4, scree_plot,  ncol = 2,  rel_widths = c(2, 1), labels = c("(A)", "(B)"), label_x = c(0.02, -0.05),label_y = 0.98)
supplement_plot

png("Documents/FIBR_temporal/am_nat_revisions/pca/supplement_pca.png", res = 300, height=4, width=8, units="in")
supplement_plot
dev.off()
