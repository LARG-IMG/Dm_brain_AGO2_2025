# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
library(GenomicRanges)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)


tss <- resize(transcripts(txdb), width = 1, fix = 'start')
genesies <- resize(genes(txdb), width = 1, fix = 'start')
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

damfilall <- dir("~/Work/projects/misc/shev/normalized_bigwig", full.names = T)

e2ago2c <- import.bw(damfilall[8])
seqlevels(e2ago2c) <- paste0("chr", seqlevels(e2ago2c))

mat5 <- normalizeToMatrix(e2ago2c, temid, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat5[rowSums(mat5) != 0,])

mat6 <- normalizeToMatrix(e2ago2c, temid_line, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat6[rowSums(mat6) != 0,])

mat10 <- normalizeToMatrix(e2ago2c, temid_ltr, value_column = "score",
                          extend = 5000, w = 50, mean_mode = "w0")
EnrichedHeatmap(mat10[rowSums(mat10) != 0,])

e2ago2m <- import.bw(damfilall[6])
seqlevels(e2ago2m) <- paste0("chr", seqlevels(e2ago2m))

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
