library(GenomicRanges)
library(rtracklayer)
library(openxlsx)
library(tidyverse)
library(EnrichedHeatmap)
library(reshape2)

lam_brain <- read.xlsx("~/Downloads/Lam_brain_domains.xlsx")
lam_brain$domain.start <- lam_brain$domain.start + 1
lam_brain_gr <- makeGRangesFromDataFrame(lam_brain)
dm3dm6 <- import.chain("dm3ToDm6.over.chain")
lb_lo <- liftOver(lam_brain_gr, dm3dm6 )
lam_brain_dm6 <- unlist(lb_lo)
lam_brain_dm6 <- lam_brain_dm6[width(lam_brain_dm6) > 300]
system("ls")
table(seqnames(unlist(lb_lo)))
table(seqnames(lam_brain_gr))
lbd_point <- resize(lam_brain_dm6, width = 1, fix = 'center')


damfilall <- dir("~/Work/projects/misc/shev/normalized_bigwig", full.names = T)

e2ago2c <- import.bw(damfilall[8])
seqlevels(e2ago2c) <- paste0("chr", seqlevels(e2ago2c))

e2ago2m1 <- import.bw(damfilall[5])
seqlevels(e2ago2m1) <- paste0("chr", seqlevels(e2ago2m1))
e2ago2m2 <- import.bw(damfilall[6])
seqlevels(e2ago2m2) <- paste0("chr", seqlevels(e2ago2m2))

m_lb_agoc <- normalizeToMatrix(e2ago2c, lbd_point, value_column = "score", extend = 5000,
                               w = 50, mean_mode = "w0")

mlbagoc_df <- colMeans(m_lb_agoc) %>% enframe(value = "Ago_control1") %>% 
  mutate(coord = seq(-5000, 5000, by = 50)[-101])

m_lb_agom2 <- normalizeToMatrix(e2ago2m2, lbd_point, value_column = "score", extend = 5000,
                               w = 50, mean_mode = "w0")

mlbagom2_df <- colMeans(m_lb_agom2) %>% enframe(value = "Ago_mut2") %>% 
  mutate(coord = seq(-5000, 5000, by = 50)[-101])

e1lamkd1 <- import.bw(damfilall[3])
seqlevels(e1lamkd1) <- paste0("chr", seqlevels(e1lamkd1))

m_lb_lamm <- normalizeToMatrix(e1lamkd1, lbd_point, value_column = "score", extend = 5000,
                               w = 50, mean_mode = "w0")

mlblamm_df <- colMeans(m_lb_lamm) %>% enframe(value = "Lam_mut1") %>% 
  mutate(coord = seq(-5000, 5000, by = 50)[-101])
ggplot(mlbagoc_df, aes(x = coord))

names(damfilall) <- c("Lam_control1",
                      "Lam_control2",
                      "LamKD1",
                      "LamKD2",
                      "Ago2_mut1",
                      "Ago2_mut2",
                      "Ago2_control1",
                      "Ago2_control2")

agg_lst <- lapply(names(damfilall), function(x){
  bopo <- import.bw(damfilall[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, lbd_point, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                 agg_lst)

agg_df_e1_m <- melt(agg_df[,-1] %>% select(coord, starts_with("Lam")), id.vars = "coord")

ggplot(agg_df_e1_m, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()

agg_df_e2_m <- melt(agg_df[,-1] %>% select(coord, starts_with("Ago")), id.vars = "coord")

ggplot(agg_df_e2_m, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()

rand_df <- data.frame(
  "chr" = sample(c("chr2L", "chr2R", "chr3L", "chr3R", "chrX"),
                 3500, replace = T),
  "start" = sample(1:20e6, 3500)
) %>% mutate(end = start)
rand_gr <- makeGRangesFromDataFrame(rand_df)

agg_lst2 <- lapply(names(damfilall), function(x){
  bopo <- import.bw(damfilall[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, rand_gr, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df2 <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                 agg_lst2)

agg_df_m2 <- melt(agg_df2[,-1], id.vars = "coord")

ggplot(agg_df_m2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("Random ranges")

m_g_agoc_k3 <- kmeans(m_g_agoc, centers = 3)
m_g_agom2 <- normalizeToMatrix(e2ago2m2, genesies, value_column = "score",
                               extend = 5000, w = 50, mean_mode = "w0")


######################
# Unique mapping
###################
damfiluq <- dir("~/Work/projects/misc/shev/normalized_bigwig_uniq/", full.names = T)

names(damfiluq) <- c("Lam_control1",
                      "Lam_control2",
                      "LamKD1",
                      "LamKD2",
                      "Ago2_mut1",
                      "Ago2_mut2",
                      "Ago2_control1",
                      "Ago2_control2")

# Plot enrichment on Lam domains in brains

agg_lst_uq <- lapply(names(damfiluq), function(x){
  bopo <- import.bw(damfiluq[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, lbd_point, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df_uq <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                 agg_lst_uq)

agg_df_uq_m_e1 <- melt(agg_df_uq[,-1] %>% select(coord, starts_with("Lam")), id.vars = "coord")

pdf("plots/e1_lam_brain_domains_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_m_e1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E1 Lam domains in Brain")
dev.off()

agg_df_uq_m_e2 <- melt(agg_df_uq[,-1] %>% select(coord, starts_with("Ago")), id.vars = "coord")

pdf("plots/e2_lam_brain_domains_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_m_e2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E2 Lam domains in Brain")
dev.off()


# Plot enrichment on random ranges
agg_lst_uq_rand <- lapply(names(damfiluq), function(x){
  bopo <- import.bw(damfiluq[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, rand_gr, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df_uq_rand <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                    agg_lst_uq_rand)

agg_df_uq_rand_m_e1 <- melt(agg_df_uq_rand[,-1] %>% select(coord, starts_with("Lam")), id.vars = "coord")

pdf("plots/e1_random_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_rand_m_e1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E1 Random Ranges")
dev.off()

agg_df_uq_rand_m_e2 <- melt(agg_df_uq_rand[,-1] %>% select(coord, starts_with("Ago")), id.vars = "coord")

pdf("plots/e2_random_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_rand_m_e2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E2 Random Ranges")
dev.off()

# Plot on inter LADs
ilads <- gaps(lam_brain_dm6)
seqlevels(te_dm6) <- paste0("chr", seqlevels(te_dm6))
ilads_no_te <- subsetByOverlaps(ilads, te_dm6, invert = T, ignore.strand = T)
ilads <- resize(ilads[width(ilads) > 9000], 1, fix = 'center')
ilads2 <- resize(ilads_no_te, 1, fix = 'center')

agg_lst_uq_ilad <- lapply(names(damfiluq), function(x){
  bopo <- import.bw(damfiluq[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, ilads2, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df_uq_ilad <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                         agg_lst_uq_ilad)

agg_df_uq_ilad_m_e1 <- melt(agg_df_uq_ilad[,-1] %>% select(coord, starts_with("Lam")), id.vars = "coord")

pdf("plots/e1_ilad2_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_ilad_m_e1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E1 interLADs no repeats")+
  ylab("Mean log2(Dam-HP1/Dam)")
dev.off()

agg_df_uq_ilad_m_e2 <- melt(agg_df_uq_ilad[,-1] %>% select(coord, starts_with("Ago")), id.vars = "coord")

pdf("plots/e2_ilad2_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_ilad_m_e2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E2 interLADs no repeats")+
  ylab("mean log2(Dam-HP1/Dam)")
dev.off()

# Now with removing heterochromatin from ilads
euc.coords <- GRanges(
  seqnames = Rle(c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")),
  ranges = IRanges(
    start = c(1, 5712495, 1, 4174279, 103614),
    end = c(22000000, 25259177, 22906900, 32074278, 23020964)
  )
)

ilads3 <- subsetByOverlaps(ilads_no_te, euc.coords, ignore.strand = T)
ilads3 <- resize(ilads3, 1, fix = 'center')

agg_lst_uq_ilad <- lapply(names(damfiluq), function(x){
  bopo <- import.bw(damfiluq[[x]])
  seqlevels(bopo) <- paste0("chr", seqlevels(bopo))
  matt <- normalizeToMatrix(bopo, ilads3, value_column = "score", extend = 5000,
                            w = 50, mean_mode = "w0")
  dff <- colMeans(matt) %>% enframe(value = x) %>% 
    mutate(coord = seq(-5000, 5000, by = 50)[-101])
})

agg_df_uq_ilad <- Reduce(function(x,y) merge(x, y, by = c("name", "coord")),
                         agg_lst_uq_ilad)

agg_df_uq_ilad_m_e1 <- melt(agg_df_uq_ilad[,-1] %>% select(coord, starts_with("Lam")), id.vars = "coord")

pdf("plots/e1_ilad3_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_ilad_m_e1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E1 interLADs no repeats euc")+
  ylab("Mean log2(Dam-HP1/Dam)")
dev.off()

agg_df_uq_ilad_m_e2 <- melt(agg_df_uq_ilad[,-1] %>% select(coord, starts_with("Ago")), id.vars = "coord")

pdf("plots/e2_ilad3_mean_enrich.pdf", width = 6, height = 6)
ggplot(agg_df_uq_ilad_m_e2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ggtitle("E2 interLADs no repeats euc")+
  ylab("mean log2(Dam-HP1/Dam)")
dev.off()
