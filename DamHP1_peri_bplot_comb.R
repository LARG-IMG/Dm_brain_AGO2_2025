library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(EnrichedHeatmap)

# 2025-07-04
# Script for making boxplots for log2(DamHP1/Dam) scores in control vs KD/Mut

# Load bedgraphs

bgdir <- "Directory with your begraph (combined) files"
bgs <- lapply(dir(bgdir, full.names = T, pattern = 'dm6'), import.bedGraph)
names(bgs) <- c("E1 Control", "E1 Lam KD", "E2 Ago2 Mut", "E2 Control")
bgs <- bgs[c(1,2,4,3)]

pericoor <- GRanges(seqnames = c("2L", "2R", "3L", "3R", "X"),
                    ranges = IRanges(start = c(22001009, 1, 22962476, 1, 22628490),
                                     end = c(23513712, 5398184, 28110227, 4552934, 23542271)))

bgperi <- lapply(bgs, function(bg) subsetByOverlaps(bg, pericoor, ignore.strand = T)$score)

pdf("~/R/plots/DamHP1_pericen_bp.pdf", width = 5, height = 6)
boxplot(bgperi[1:2], outline = F)
boxplot(bgperi[3:4], outline = F)
dev.off()

wilcox.test(bgperi[[4]], bgperi[[1]], alternative = 'less')
plot(density(bgperi[[2]]))
