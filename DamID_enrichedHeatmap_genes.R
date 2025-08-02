# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
library(GenomicRanges)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(reshape2)
library(cowplot)


tss <- resize(transcripts(txdb), width = 1, fix = 'start')
genesies <- resize(genes(txdb), width = 1, fix = 'start')
# names(genesies) <- genesies$gene_id
seqlevels(damnorm2) <- paste0("chr", seqlevels(damnorm2))
mat1 = normalizeToMatrix(damnorm2, tss, value_column = "score", 
                         extend = 5000, mean_mode = "w0", w = 50)
EnrichedHeatmap(mat1[rowSums(mat1) != 0,], name = "DamHP1")
temid <- resize(te_dm6[te_dm6$class == "DNA"], width = 1, fix = "center")
seqlevels(temid) <- paste0("chr", seqlevels(temid))
mat2 <- normalizeToMatrix(damnorm2, temid, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat2[rowSums(mat2) != 0,])
sum(rowSums(mat2) == 0)

temid_line <- resize(te_dm6[te_dm6$class == "LINE"], width = 1, fix = "center")
seqlevels(temid_line) <- paste0("chr", seqlevels(temid_line))
mat3 <- normalizeToMatrix(damnorm2, temid_line, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat3[rowSums(mat3) != 0,])

temid_ltr <- resize(te_dm6[te_dm6$class == "LTR"], width = 1, fix = "center")
seqlevels(temid_ltr) <- paste0("chr", seqlevels(temid_ltr))
mat4 <- normalizeToMatrix(damnorm2, temid_ltr, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat4[rowSums(mat4) != 0,])

damfilall <- dir("normalized_bigwig", full.names = T)

e2ago2c <- import.bw(damfilall[8])
seqlevels(e2ago2c) <- paste0("chr", seqlevels(e2ago2c))

e2ago2m1 <- import.bw(damfilall[5])
seqlevels(e2ago2m1) <- paste0("chr", seqlevels(e2ago2m1))
e2ago2m2 <- import.bw(damfilall[6])
seqlevels(e2ago2m2) <- paste0("chr", seqlevels(e2ago2m2))

m_g_agoc <- normalizeToMatrix(e2ago2c, genesies, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
m_g_agoc_k3 <- kmeans(m_g_agoc, centers = 3)
m_g_agom2 <- normalizeToMatrix(e2ago2m2, genesies, value_column = "score",
                               extend = 5000, w = 50, mean_mode = "w0")
all(rownames(m_g_agoc) == rownames(m_g_agom2))
all(rownames(m_g_agom2) == names(m_g_agoc_k3$cluster))

m_g_agom1 <- normalizeToMatrix(e2ago2m1, genesies, value_column = "score",
                               extend = 5000, w = 50, mean_mode = "w0")

mgack_means <- by(m_g_agoc, m_g_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pc1 <- ggplot(mgack_means, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-1, 1))+
  ggtitle("Ago2+ TSS by genes")
pdf("plots/Ago2_control_kmeans3_TSS.pdf", width = 5, height = 5)
pc1
dev.off()

mgamk_means1 <- by(m_g_agom1, m_g_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pm1 <- ggplot(mgamk_means1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-1, 1))+
  ggtitle("Ago2- (rep1) TSS by genes")

mgamk_means2 <- by(m_g_agom2, m_g_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pm2 <- ggplot(mgamk_means2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-1, 1))+
  ggtitle("Ago2- (rep2) TSS by genes")

pm12c <- plot_grid(pc1, pm1, pm2, nrow = 1)
pdf("plots/Ago2_mut_cont_rep12_TSS_kmeans3.pdf", height = 5, width = 15)
pm12c
dev.off()


cluster_top_id <- names(m_g_agoc_k3$cluster)[m_g_agoc_k3$cluster == 3]

# Convert l1 FBGn ID to ENTREZ
ct_id_conv <- bitr(cluster_top_id, 
                         fromType = "FLYBASE", 
                         toType = "ENTREZID", 
                         OrgDb = org.Dm.eg.db)
# GO for biological processes
go_r_hp1enr_bp <- enrichGO(gene = ct_id_conv$ENTREZID,
                       OrgDb = org.Dm.eg.db,
                       ont = "BP",  # You can also use "MF" or "CC"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)

dotplot(go_r_hp1enr_bp, showCategory = 10)+
  ggtitle("Biological Process")

#######################
# Mean enrichment TEs
######################
te_proper <- te_dm6[te_dm6$class %in% c("DNA", "LINE", "LTR") & seqnames(te_dm6) %in% c("2L", "2R", "3L", "3R", "4", "X", "Y")]
te_centers <- resize(te_proper, width = 1, fix = 'center')
seqlevels(te_centers) <- paste0("chr", seqlevels(te_centers))
m_te_agoc <- normalizeToMatrix(e2ago2c, te_centers, value_column = "score",
                              extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(m_te_agoc, row_km = 2)
m_te_agoc_k3 <- kmeans(m_te_agoc, centers = 2)
m_te_agom2 <- normalizeToMatrix(e2ago2m2, te_centers, value_column = "score",
                               extend = 5000, w = 50, mean_mode = "w0")
all(rownames(m_te_agoc) == rownames(m_te_agom2))
all(rownames(m_te_agom2) == names(m_te_agoc_k3$cluster))

m_te_agom1 <- normalizeToMatrix(e2ago2m1, te_centers, value_column = "score",
                               extend = 5000, w = 50, mean_mode = "w0")

mgack_means <- by(m_te_agoc, m_te_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pc1 <- ggplot(mgack_means, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-0.5, 2))+
  ggtitle("Ago2+ TEs")
pdf("plots/Ago2_control_kmeans3_TE_center.pdf", width = 5, height = 5)
pc1
dev.off()

mgamk_means1 <- by(m_te_agom1, m_te_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pm1 <- ggplot(mgamk_means1, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-0.5, 2))+
  ggtitle("Ago2- (rep1) TEs")

mgamk_means2 <- by(m_te_agom2, m_te_agoc_k3$cluster, colMeans) %>% do.call(rbind, .) %>% 
  t() %>% as.tibble() %>% mutate("coord" = seq(-5000, 5000, by = 50)[-101]) %>% 
  melt(id.vars = "coord")
pm2 <- ggplot(mgamk_means2, aes(x = coord, y = value, col = variable))+
  geom_line()+
  theme_bw()+
  ylim(c(-0.5, 2))+
  ggtitle("Ago2- (rep2) TEs")

pm12c <- plot_grid(pc1, pm1, pm2, nrow = 1)
pdf("plots/Ago2_mut_cont_rep12_TEs_kmeans2.pdf", height = 5, width = 15)
pm12c
dev.off()















EnrichedHeatmap(m_g_agoc, row_split = factor(m_g_agoc_k3$cluster, levels = c(3,1,2)))+
EnrichedHeatmap(m_g_agom2, row_split = factor(m_g_agoc_k3$cluster, levels = c(3,1,2)))

mat6 <- normalizeToMatrix(e2ago2c, temid_line, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
peep <- EnrichedHeatmap(mat6[rowSums(mat6) != 0,])

mat10 <- normalizeToMatrix(e2ago2c, temid_ltr, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat10[rowSums(mat10) != 0,])


mat7 <- normalizeToMatrix(e2ago2m, temid, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat7[rowSums(mat7) != 0,])

mat8 <- normalizeToMatrix(e2ago2m, temid_line, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat8[rowSums(mat8) != 0,])

mat9 <- normalizeToMatrix(e2ago2m, temid_ltr, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat9[rowSums(mat9) != 0,])

mat11 <- normalizeToMatrix(e2ago2m, tss, value_column = "score",
                           extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat11[rowSums(mat11) != 0,])

mat12 <- normalizeToMatrix(e2ago2c, tss, value_column = "score",
                           extend = 5000, w = 50, mean_mode = "w0")

pp_ltr_c <- EnrichedHeatmap(mat10[rowSums(mat10) != 0,], column_title = "Ago2+",
                            axis_name = c(-5000, "LTR center", 5000), name = "log2(DamHP1/Dam)",
                            top_annotation = HeatmapAnnotation(
                              enrich = anno_enriched(
                                ylim = c(0, 1.05),
                                axis_param = list(side = "left"))))
pp_ltr_m <- EnrichedHeatmap(mat9[rowSums(mat9) != 0,], column_title = "Ago2-",
                            axis_name = c(-5000, "LTR center", 5000), 
                            top_annotation = HeatmapAnnotation(
                              enrich = anno_enriched(ylim = c(0, 1.05))
                            ),
                            show_heatmap_legend = F)
pdf("plots/damhp1_ago2_ltr_heatmap.pdf", width = 6, height = 12)
pp_ltr_c+pp_ltr_m
dev.off()

pp_tss_c <- EnrichedHeatmap(mat12[rowSums(mat12) != 0,], column_title = "Ago2+",
                            axis_name = c(-5000, "TSS", 5000), top_annotation = HeatmapAnnotation(
                              enrich = anno_enriched(
                                ylim = c(-0.4, 0),
                                axis_param = list(side = "left"))),
                            name = "log2(DamHP1/Dam)")
pp_tss_m <- EnrichedHeatmap(mat11[rowSums(mat11) != 0,], column_title = "Ago2-",
                            axis_name = c(-5000, "TSS", 5000), top_annotation = HeatmapAnnotation(
                              enrich = anno_enriched(
                                ylim = c(-0.4, 0))),
                            show_heatmap_legend = F)
pdf("plots/damhp1_ago2_tss_heatmap.pdf", width = 6, height = 12)
pp_tss_c+pp_tss_m
dev.off()

partition = paste0("cluster", kmeans(mat11, centers = 3)$cluster)

mat12_k3 <- kmeans(mat12, centers = 3)

EnrichedHeatmap(mat10[rowSums(mat10) != 0,], name = "log2(DamHP1/Dam)", row_km = 3,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                  column_title = "Ago2+")+
  EnrichedHeatmap(mat9[rowSums(mat9) != 0,], name = "log2(DamHP1/Dam)",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                  column_title = "Ago2-")

boop <- EnrichedHeatmap(mat12[rowSums(mat12) != 0,], name = "log2(DamHP1/Dam)", row_km = 3,
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                column_title = "Ago2+")+
  EnrichedHeatmap(mat11[rowSums(mat11) != 0,], name = "log2(DamHP1/Dam)",
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                  column_title = "Ago2-")
