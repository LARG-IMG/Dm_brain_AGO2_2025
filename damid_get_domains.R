source("~/R/damid_hmm_functions.R")
library(rtracklayer)

damfilall <- dir("~/Work/projects/misc/shev/normalized_bigwig", full.names = T)
damfiluq <- dir("~/Work/projects/misc/shev/normalized_bigwig_uniq", full.names = T)

exps <- basename(damfilall) %>% gsub("_Dam.*", "", .) %>% unique()

damid_log2_all <- lapply(exps, function(expo){
  reps <- damfilall[grepl(expo, damfilall)]
  rep1 <- import.bw(reps[1])
  rep2 <- import.bw(reps[2])
  hmm1 <- run_hmm_on_granges(rep1, n_states = 3,
                             collapse_enriched = T, by_chr = F)
  hmm2 <- run_hmm_on_granges(rep2, n_states = 3,
                             collapse_enriched = T, by_chr = F)
  hmm_12 <- GenomicRanges::intersect(hmm1, hmm2, ignore.strand = T)
  return(hmm_12)
})

names(damid_log2_all) <- exps

damid_dom_uq <- lapply(exps, function(expo){
  reps <- damfiluq[grepl(expo, damfiluq)]
  rep1 <- import.bw(reps[1])
  rep2 <- import.bw(reps[2])
  hmm1 <- run_hmm_on_granges(rep1, n_states = 3,
                             collapse_enriched = T, by_chr = F)
  hmm2 <- run_hmm_on_granges(rep2, n_states = 3,
                             collapse_enriched = T, by_chr = F)
  hmm_12 <- GenomicRanges::intersect(hmm1, hmm2, ignore.strand = T)
  return(hmm_12)
})

names(damid_dom_uq) <- exps

sapply(damid_log2_all, length)
sapply(damid_log2_all, function(x) sum(width(x)))

sapply(damid_dom_uq, length)
sapply(damid_dom_uq, function(x) sum(width(x)))

dir.create("~/Work/projects/misc/shev/bed")
for (expo in exps){
  export.bed(damid_log2_all[[expo]],
             file.path("~/Work/projects/misc/shev/bed",
                       paste0(expo, "_DamHP1_HMM3_domains_all.bed")))
  export.bed(damid_dom_uq[[expo]],
             file.path("~/Work/projects/misc/shev/bed",
                       paste0(expo, "_DamHP1_HMM3_domains_uq.bed")))
}
