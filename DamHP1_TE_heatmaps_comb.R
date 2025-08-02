library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(EnrichedHeatmap)

# 2025-07-02
# Script for making heatmaps of TE log2(DamHP1/Dam) scores
load("RData/te.hier.added.envs.RData", verb = T)

# Load bedgraphs

bgdir <- "Path to your bedgraph (combined) files"
bgs <- lapply(dir(bgdir, full.names = T, pattern = 'TE'), import.bedGraph)
names(bgs) <- c("E1 Control", "E1 Lam KD", "E2 Ago2 Mut", "E2 Control")

# Make a matrix with average RPMs per TE
mtx <- sapply(bgs, function(bg){
  heh <- as.list(split(bg, seqnames(bg)))
  sapply(heh, function(x) mean(x$score))
})
mtx <- mtx[!grepl("ste", rownames(mtx),ignore.case = T),]

# Heatmap of plain RPM values for E1
heatmap(
  mtx[,1:2],
  Colv = NA,                # don’t reorder columns
  scale = 'none',
  col = colorRampPalette(c("blue","white","red"))(100),
  margins = c(5,8),
  ylab = "Lamin KD",
  main = "TEs RPMs"
)

# Heatmap of plain RPM values for E2
heatmap(
  mtx[,3:4],
  Colv = NA,                # don’t reorder columns
  scale = 'none',
  col = colorRampPalette(c("blue","white","red"))(100),
  margins = c(5,8),
  ylab = "TEs",
  main = "Ago2 Mutation"
)

mtxe1 <- mtx[,1:2]
mtxe2 <- mtx[,3:4]

# Heatmap
pdf("plots/DamHP1_agom_comb_TE_enrich.pdf", width = 6, height = 12)
heatmap(
  mtxe2,
  Rowv = NA,
  revC = T,
  scale = "none",           # we already scaled
  col = colorRampPalette(c("blue","white","red"))(100),
  margins = c(3,7),
  cexCol = 1,
  main = "TE RPMs (row-wise z-scores)"
)
dev.off()



# max value in the matrix
m <- max(abs(mtxe1))

my_colors <- colorRampPalette(c("navy", "white", "firebrick2"))(100)

# 3) Define breaks so that zero falls exactly between colors
breaks <- seq(-m, m, length.out = length(my_colors) + 1)

# 4) Plot with pheatmap (it draws a legend by default)
pheatmap(
  mtxe1,
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
  main           = "log2(DamHP1/Dam) per TE",
  legend         = TRUE,      # default, but explicit here
  name = 'log2(DamHP1/Dam)',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1)))
)

mtxe1diff <- mtxe1[,2] - mtxe1[,1]
mtxe1diff <- sort(mtxe1diff, decreasing = T)

pheatmap(
  as.matrix(mtxe1diff),
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 6,
  fontsize_col   = 8,
  cellwidth = 20,
  main           = "HP1 enrichment in Lam KD - Control",
  legend         = TRUE,      # default, but explicit here
  name = 'log2(DamHP1/Dam)',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1))),
  filename = "plots/DamHP1_E1LamKD_TE_diff.pdf",
  width = 5,
  height = 10
)


dfe1 <- as.data.frame(mtxe1) %>% setNames(c("Control", "LamKD")) %>% 
  rownames_to_column("TE") %>%   mutate(diff = LamKD - Control)
dfe1 <- merge(dfe1, te3.hier, by.x = 'TE', by.y = 'id', all.x = T)
dfe1 <- dfe1 %>% 
  dplyr::select(-Alias)
library(openxlsx)
write.xlsx(dfe1, "tables/lamkd_damhp1_TE_dif.xlsx")


mtxe2diff <- mtxe2[,1] - mtxe2[,2]
mtxe2diff <- sort(mtxe2diff)

pheatmap(
  as.matrix(mtxe2diff),
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 6,
  fontsize_col   = 8,
  cellwidth = 20,
  main           = "HP1 enrichment in Ago2 Mut - Control",
  legend         = TRUE,      # default, but explicit here
  name = 'log2(DamHP1/Dam)',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1))),
  filename = "plots/DamHP1_e2Ago2Mut_TE_diff.pdf", width = 5, height = 10
)


dfe2 <- as.data.frame(mtxe2[,2:1]) %>% setNames(c("Control", "Ago2Mut")) %>% 
  rownames_to_column("TE") %>%   mutate(diff = Ago2Mut - Control)
dfe2 <- merge(dfe2, te3.hier, by.x = 'TE', by.y = 'id', all.x = T)
dfe2 <- dfe2 %>% 
  dplyr::select(-Alias)
write.xlsx(dfe2, "tables/ago2mut_damhp1_TE_dif.xlsx")
