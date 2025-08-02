library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)
library(ChIPseeker)
library(cowplot)
library(rtracklayer)
library(EnrichedHeatmap)
library(reshape2)

# Plot HP1 enrichment over TE insertions in Ago2 KD and control

# Load TE insertions longer than 100 bp
load("RData/Mcclintock_dm3_insertions.RData", verbose = T)

load("RData/ago2_dm3_tracks_multi.RData", verbose = T)
load("RData/lam_dm3_tracks_multi.RData", verbose = T)

tespl <- as.list(split(m.a.gr, m.a.gr$name))

enh_ver_point <- resize(enh_ver, width = 1, fix = 'center')

hp1_mdg3_mean <- lapply(agotracks, function(bg){
  mat <- normalizeToMatrix(bg, tespl[['mdg3']], value_column = 'score',
                           extend = 1000, w = 50, mean_mode = 'w0')
  return(colMeans(mat))
})

hp1_mdg3_mean_df <- do.call(rbind, hp1_mdg3_mean)

hp1_mdg3 <- t(hp1_mdg3_mean_df) %>% as.data.frame() %>% mutate(x = 1:57) %>%  melt(id.vars = "x")

pdf("plots/uq_bw_mean_TEs.pdf")
ggplot(hp1_mdg3, aes(x = x, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("Mean Enrichment on mdg3")
dev.off()


hp1_copia_mean <- lapply(agotracks, function(bg){
  mat <- normalizeToMatrix(bg, tespl[['copia']], value_column = 'score',
                           extend = 1000, w = 50, mean_mode = 'w0')
  return(colMeans(mat))
})

hp1_copia_mean_df <- do.call(rbind, hp1_copia_mean)

hp1_copia <- t(hp1_copia_mean_df) %>% as.data.frame() %>% mutate(x = 1:67) %>%  melt(id.vars = "x")

pdf("plots/uq_bw_mean_TEs.pdf")
ggplot(hp1_copia, aes(x = x, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("Mean Enrichment on copia")
dev.off()

gypsies <- tespl[['gypsy']]
gypsies <- resize(gypsies, 1, fix = 'start')
hp1_gypsy_mean <- lapply(agotracks, function(bg){
  mat <- normalizeToMatrix(bg, gypsies, value_column = 'score',
                           extend = c(1000, 10000), w = 100, mean_mode = 'w0')
  return(colMeans(mat))
})

hp1_gypsy_mean_df <- do.call(rbind, hp1_gypsy_mean)

hp1_gypsy <- t(hp1_gypsy_mean_df) %>% as.data.frame() %>% mutate(x = seq(-1000, 10000, by = 100)[-11]) %>%  melt(id.vars = "x")


ggplot(hp1_gypsy, aes(x = x, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  xlab("bp from TE insertion start")
  ggtitle("Mean Enrichment on gypsy")
  
  
###############
# plot for all tes
##################
pdf("TE_enrichment_plots.pdf", width = 8, height = 6)

# Loop over each TE insertion set in 'tespl'
for (te_name in names(tespl)) {
  
  # Extract and resize the GRanges for this TE
  te_ranges <- tespl[[te_name]]
  te_ranges <- resize(te_ranges, width = 1, fix = "start")
  
  # For each BigWig track in 'agotracks', normalize around the TE start
  te_mean_list <- lapply(agotracks, function(bg) {
    mat <- normalizeToMatrix(
      bg, 
      te_ranges, 
      value_column   = "score",
      extend         = c(1000, 10000), 
      w              = 100, 
      mean_mode      = "w0"
    )
    colMeans(mat)
  })
  
  # Combine into a single matrix and convert to data.frame
  te_mean_df <- do.call(rbind, te_mean_list)
  te_plot_df <- t(te_mean_df) %>%
    as.data.frame() %>%
    mutate(
      x = seq(-1000, 10000, by = 100)[-11]  # drop the zero‐bin if needed
    ) %>%
    melt(id.vars = "x")
  
  # Build the ggplot
  p <- ggplot(te_plot_df, aes(x = x, y = value, col = variable)) +
    geom_line() +
    theme_bw() +
    xlab("bp from TE insertion start") +
    ylab("Mean normalized score") +
    ggtitle(paste("Mean Enrichment on", te_name))
  
  # Print to the current PDF page
  print(p)
}

dev.off()
  
  

hp1_jockey_mean <- lapply(agotracks, function(bg){
  mat <- normalizeToMatrix(bg, tespl[['jockey']], value_column = 'score',
                           extend = 1000, w = 50, mean_mode = 'w0')
  return(colMeans(mat))
})

hp1_jockey_mean_df <- do.call(rbind, hp1_jockey_mean)

hp1_jockey <- t(hp1_jockey_mean_df) %>% as.data.frame() %>% mutate(x = 1:51) %>%  melt(id.vars = "x")


ggplot(hp1_jockey, aes(x = x, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("Mean Enrichment on jockey")
