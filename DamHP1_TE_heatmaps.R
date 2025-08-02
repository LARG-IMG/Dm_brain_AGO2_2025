library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(EnrichedHeatmap)

# 2025-06-24
# Script for making heatmaps of TE log2(DamHP1/Dam) scores


# Load bedgraphs

bgdir <- "Path to directory with TE bedgraph files"
bgs <- lapply(dir(bgdir, full.names = T), import.bedGraph)
names(bgs) <- c("AgoMut3", "AgoMut1", "Agomut2", "Control1", "Control2")
bgs <- bgs[-2]
# Make a matrix with average RPMs per TE
mtx <- sapply(bgs, function(bg){
  heh <- as.list(split(bg, seqnames(bg)))
  sapply(heh, function(x) mean(x$score))
})

# Heatmap of plain RPM values
heatmap(
  mtx,
  Colv = NA,                # don’t reorder columns
  col = colorRampPalette(c("blue","white","red"))(100),
  margins = c(5,8),
  xlab = "Samples", ylab = "TEs",
  main = "TEs RPMs (row-wise z-scores)"
)

# Convert to Z-scores (row-wise)
zmtx <- t(apply(mtx, 1, function(x) (x - mean(x))/sd(x)))
# Heatmap
pdf("plots/DamHP1_agom_TE_enrich.pdf", width = 6, height = 12)
heatmap(
  zmtx,
  Rowv = NA,
  revC = T,
  scale = "none",           # we already scaled
  col = colorRampPalette(c("blue","white","red"))(100),
  margins = c(3,7),
  cexCol = 1,
  main = "TE RPMs (row-wise z-scores)"
)
dev.off()

# Average RPMs for replcates and make a log2FC heatmap

mtxc <- sapply(list(c(1,2), c(3,4)), function(x) rowMeans(mtx[,x]))
colnames(mtxc) <- c("AgoMut", "Control")
mtxc <- mtxc[order(mtxc[,2] - mtxc[,1], decreasing = T), ]

heatmap(
  mtxc,
  Rowv = NA,
  Colv = NA,
  #revC = T,
  scale = "none",           # we already scaled
  col = colorRampPalette(c("blue", "white","red2"))(50),
  margins = c(3,7),
  cexCol = 1,
  main = "TE RPMs"
)

# max value in the matrix
m <- max(abs(mtxc))

my_colors <- colorRampPalette(c("navy", "white", "firebrick2"))(100)

# 3) Define breaks so that zero falls exactly between colors
breaks <- seq(-m, m, length.out = length(my_colors) + 1)

# 4) Plot with pheatmap (it draws a legend by default)
pheatmap(
  mtxc,
  color          = my_colors,
  breaks         = breaks,
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 5,
  fontsize_col   = 8,
  main           = "Avg log2(DamHP1/Dam) per TE",
  legend         = TRUE,      # default, but explicit here
  name = 'log2(DamHP1/Dam)',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1)))
)

mtx1 <- mtxc[,2] - mtxc[,1]
pdf("plots/DamHP1_agom_TE_diff.pdf", width = 5, height = 10)
pheatmap(
  as.matrix(mtx1),
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 6,
  fontsize_col   = 8,
  cellwidth = 20,
  main           = "Avg log2(DamHP1/Dam) Control - Ago2Mut",
  legend         = TRUE,      # default, but explicit here
  name = 'log2(DamHP1/Dam)',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1)))
)
dev.off()

mtx1e <- enframe(mtx1) %>% setNames(c("TE", "diff"))
mtx1e <- merge(mtx1e, te3.hier, by.x = 'TE', by.y = 'id', all.x = T)
mtx1e <- mtx1e %>% 
  select(-Alias)
library(openxlsx)
write.xlsx(mtx1e, "tables/agomut_damhp1_TE_dif.xlsx")
load("RData/te.hier.added.envs.RData", verb = T)

avgbg <- function(bg1, bg2){
  heh <- mcolAsRleList(bg1, 'score')
  huh <- mcolAsRleList(bg2, 'score')
  hyh <- lapply(1:length(heh), function(x) log2(((2^heh[[x]]) + (2^huh[[x]]))/2))
  names(hyh) <- names(huh)
  as(RleList(hyh), "GRanges")
}

agomut <- avgbg(bgs[[1]], bgs[[2]])
contro <- avgbg(bgs[[3]], bgs[[4]])

export.bedGraph(agomut, "bed/E2AGOM_DamHP1_TE_log2Dam.bedgraph")
export.bedGraph(contro, "bed/E2C_DamHP1_TE_log2Dam.bedgraph")
