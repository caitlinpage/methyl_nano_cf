library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(plyranges)
library(tidyr)
library(stringr)
library(dplyr)
library(OpTiles)
library(edgeR)

seqs <- seqnames(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[1:22]

cg_sites <- lapply(seqs, function(x) {
  cbind(data.frame(Biostrings::matchPattern("CG", BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38[[x]])),
        seqnames = x
  )}) %>% bind_rows()

process_wgbs_short <- function(bismark, bed = NULL) {
  read <- bismark
  names(read) <- c("read_id", "strand", "seqnames", "start", "meth_status")
  read$end <- read$start

  read_overlap <- read
  read_overlap[, meth_status := ifelse(meth_status == "Z", "M", "U")]
  read_overlap[, meth_pattern := paste0(meth_status, collapse = ""), by = read_id]
  read_overlap[, num_meth := str_count(meth_pattern, pattern = "M")]
  read_overlap[, total := str_length(meth_pattern)]
  read_overlap[, alpha := num_meth/total]
  read_overlap[, genom_positions := paste0(start, collapse = ","), by = read_id]

  read_overlap <- read_overlap %>% .[order(.$read_id),]

  read_overlap[, read_start := min(start), by = read_id]
  read_overlap[, read_end := max(start), by = read_id]

  read_overlap
}


bis <- fread("/researchers/caitlin.page/cf_nano/fastq_tothill_emseq_emseq_cancer/CpG_OT_PRJ230223_LPRJ230223_S3_L004_R1_001_bismark_bt2_pe.txt")
a1 <- process_wgbs(bis)
a1 <- data.frame(a1)
a1 <- a1 %>% filter(seqnames %in% seqs)
a1$sample <- "cancer_0223"
saveRDS(a1, "/researchers/caitlin.page/cf_nano/r_output/emseq_cancer_0223.rds")

bis <- fread("/researchers/caitlin.page/cf_nano/fastq_tothill_emseq_emseq_cancer/CpG_OT_PRJ230224_LPRJ230224_S4_L004_R1_001_bismark_bt2_pe.txt")
a4 <- process_wgbs(bis)
a4 <- data.frame(a4)
a4 <- a4 %>% filter(seqnames %in% seqs)
a4$sample <- "cancer_0224"
saveRDS(a4, "/researchers/caitlin.page/cf_nano/r_output/emseq_cancer_0224.rds")

bis <- fread("/researchers/caitlin.page/cf_nano/fastq_tothill_emseq_emseq_cancer/CpG_OT_PRJ230225_LPRJ230225_S5_L004_R1_001_bismark_bt2_pe.txt")
a2 <- process_wgbs(bis)
a2 <- data.frame(a2)
a2 <- a2 %>% filter(seqnames %in% seqs)
a2$sample <- "cancer_0225"
saveRDS(a2, "/researchers/caitlin.page/cf_nano/r_output/emseq_cancer_0225.rds")

bis <- fread("/researchers/caitlin.page/cf_nano/fastq_tothill_emseq_emseq_cancer/CpG_OT_PRJ230454_LPRJ230454_S1_L003_R1_001_bismark_bt2_pe.txt")
a3 <- process_wgbs(bis)
a3 <- data.frame(a3)
a3 <- a3 %>% filter(seqnames %in% seqs)
a3$sample <- "cancer_0454"
saveRDS(a3, "/researchers/caitlin.page/cf_nano/r_output/emseq_cancer_0454.rds")


cols <- c("seqnames", "read_start", "read_end", "strand", "read_id", "sample", "meth_pattern", "num_meth", "total", "alpha", "genom_positions")

all_samples <- rbind(a1[,cols],
                     a2[,cols],
                     a3[,cols],
                     a4[,cols])
all_samples$start <- all_samples$read_start
all_samples$end <- all_samples$read_end
all_samples <- distinct(all_samples)
saveRDS(all_samples, "/researchers/caitlin.page/cf_nano/r_output/emseq_all_cancer.rds")


