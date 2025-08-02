library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(EnrichedHeatmap)
library(pheatmap)

# 2025-07-05
# Script for making heatmaps of TE log2(DamHP1/Dam) scores
load("RData/te.hier.added.envs.RData", verb = T)

# Load bedgraphs

bgdir <- "Directory with your bedgraph files"
bgs <- lapply(dir(bgdir, full.names = T), import.bedGraph)
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
  col = colorRampPalette(c("white","red"))(100),
  margins = c(5,8),
  ylab = "Lamin KD",
  main = "TEs RPMs"
)

# Heatmap of plain RPM values for E2
heatmap(
  mtx[,4:3],
  Colv = NA,                # don’t reorder columns
  scale = 'none',
  col = colorRampPalette(c("white","red"))(100),
  margins = c(5,8),
  ylab = "TEs",
  main = "Ago2 Mutation"
)

mtxe1 <- mtx[,1:2]
mtxe2 <- mtx[,4:3]

# Heatmap
pdf("~/R/plots/Dam_agom_comb_TE_enrich.pdf", width = 6, height = 12)
heatmap(
  mtxe2[order(mtxe2[,2]),],
  Rowv = NA,
  Colv = NA,
  revC = F,
  scale = "none",           # we already scaled
  col = colorRampPalette(c("white","red"))(50),
  margins = c(3,7),
  cexCol = 1,
  main = "TE RPMs"
)
dev.off()



# max value in the matrix
m <- max(abs(mtxe1))

my_colors <- colorRampPalette(c("blue", "white", "firebrick2"))(50)

# 3) Define breaks so that zero falls exactly between colors
breaks <- seq(-m, m, length.out = length(my_colors) + 1)

# 4) Plot with pheatmap (it draws a legend by default)
pheatmap(
  mtxe1,
  color          = my_colors,
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 5,
  fontsize_col   = 8,
  main           = "Dam RPM per TE",
  legend         = TRUE,      # default, but explicit here
  name = 'Dam RPM',
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1)))
)

# mtxe1diff <- log2((mtxe1[,2] + 1)/(mtxe1[,1] + 1))
mtxe1diff <- log2((mtxe1[,2])/(mtxe1[,1]))[mtxe1[,2] > 0]
mtxe1diff <- sort(mtxe1diff, decreasing = F)
m <- max(abs(mtxe1diff))
my_colors <- colorRampPalette(c("blue", "white", "firebrick2"))(50)
breaks <- seq(-m, m, length.out = length(my_colors) + 1)

pheatmap(
  as.matrix(mtxe1diff),
  breaks = breaks,
  color = my_colors,
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 6,
  fontsize_col   = 8,
  cellwidth = 20,
  main           = "Dam enrichment in Lam KD - Control",
  legend         = TRUE,      # default, but explicit here
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1))),
  filename = "plots/Dam_E1LamKD_TE_diff.pdf",
  width = 5, height = 10
)


dfe1 <- as.data.frame(mtxe1[mtxe1[,2] > 0,]) %>% setNames(c("Control", "LamKD")) %>% 
  rownames_to_column("TE") %>%   mutate(log2_KD_Control = log2((LamKD)/(Control)))
dfe1 <- merge(dfe1, te3.hier, by.x = 'TE', by.y = 'id', all.x = T)
dfe1 <- dfe1 %>% 
  dplyr::select(-Alias)
library(openxlsx)
write.xlsx(dfe1, "~/Work/projects/misc/shev/lamkd_dam_TE_dif.xlsx")


mtxe2diff <- log2((mtxe2[,2])/(mtxe2[,1]))[mtxe2[,2] > 0]
mtxe2diff <- sort(mtxe2diff, decreasing = T)
m <- max(abs(mtxe2diff))
breaks <- seq(-m, m, length.out = length(my_colors) + 1)
my_colors <- colorRampPalette(c("blue", "white", "firebrick2"))(50)

pheatmap(
  as.matrix(mtxe2diff),
  breaks = breaks,
  color = my_colors,
  cluster_rows = FALSE,
  cutree_rows = 2,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize       = 10,
  fontsize_row   = 6,
  fontsize_col   = 8,
  cellwidth = 20,
  main           = "Dam enrichment in Ago2 Mut vs Control",
  legend         = TRUE,      # default, but explicit here
  legend_labels  = c(paste0("-", round(m,1)), "0", paste0("+", round(m,1))),
  filename = "plots/Dam_E2Ago2Mut_TE_diff.pdf",
  width = 5, height = 10
)

dfe2 <- as.data.frame(mtxe2[mtxe2[,2] > 0,]) %>% setNames(c("Control", "Ago2Mut")) %>% 
  rownames_to_column("TE") %>%   mutate(log2_Mut_Control = log2((Ago2Mut)/(Control)))
dfe2 <- merge(dfe2, te3.hier, by.x = 'TE', by.y = 'id', all.x = T)
dfe2 <- dfe2 %>% 
  dplyr::select(-Alias)

write.xlsx(dfe2, "tables/ago2mut_dam_TE_dif.xlsx")
