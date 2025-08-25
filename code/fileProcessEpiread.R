library(plyranges)
library(dplyr)
library(tidyr)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

cg_sites <- cbind(data.frame(Biostrings::matchPattern("CG", BSgenome.Hsapiens.UCSC.hg38[["chr22"]])), seqnames = "chr22") %>%
  mutate(pos = paste0(seqnames, "-", start)) %>%
  relocate(pos, seqnames)
cg_sites$index <- 1:nrow(cg_sites)

modkit_files <- list.files("/researchers/caitlin.page/cf_nano/modkit_output")
modkit_files <- paste0("/researchers/caitlin.page/cf_nano/modkit_output/", modkit_files)

healthy_samples <- modkit_files[grep("filt_ec", modkit_files)][c(1:2,6:9)]
luad_samples <- modkit_files[grep("filt_ec", modkit_files)][c(3:5,10:12)]
modkit_colnames <- fread(modkit_files[1])
modkit_colnames <- colnames(modkit_colnames)

luad_headers <- gsub("/researchers/caitlin.page/cf_nano/modkit_output/filt_ec_", "", luad_samples)
luad_headers <- gsub(".tsv", "", luad_headers)
luad_epireads <- list()

for(i in 1:length(luad_samples)) {
  make_epiread <- read.table(luad_samples[i], header = FALSE)
  colnames(make_epiread) <- modkit_colnames
  make_epiread$seqnames <- make_epiread$chrom
  make_epiread$start <- make_epiread$alignment_start
  make_epiread$end <- make_epiread$alignment_end
  make_epiread <- make_epiread[grep("..CG.", make_epiread$query_kmer),]

  make_epiread <- make_epiread %>% group_by(read_id) %>% mutate(num_cg = n()) %>% ungroup()
  overlaps <- find_overlaps(as_granges(make_epiread), as_granges(cg_sites)) %>% data.frame() %>%
    group_by(read_id) %>% summarise(min_index = min(index)) %>% ungroup() %>% data.frame()
  make_epiread <- make_epiread %>% mutate(index_cg = overlaps[match(.$read_id, overlaps$read_id), "min_index"]) %>% data.frame()

  make_epiread <- make_epiread %>% group_by(read_id, seqnames, start, end, chrom, alignment_start, alignment_end, read_length, num_cg, index_cg) %>%
    mutate(call_code2 = ifelse(call_code == "m", "C", "T")) %>%
    summarise(meth_pattern = paste(call_code2, collapse = "")) %>% .[order(.$alignment_start),] %>% ungroup()

  cpg_positions <- data.frame(cpg_positions = as.character())
  for(x in 1:nrow(make_epiread)) {
    cpg_positions <- rbind(cpg_positions, paste0(seq(from = make_epiread$index_cg[x], length = make_epiread$num_cg[x]), collapse = ","))

  }
  colnames(cpg_positions) <- "cpg_positions"
  make_epiread <- cbind(make_epiread, cpg_positions)

  genom_positions <- data.frame(genom_positions = as.character())
  for(y in 1:nrow(make_epiread)) {
    genom_positions <- rbind(genom_positions, paste0(cg_sites[seq(from = make_epiread$index_cg[y], length = make_epiread$num_cg[y]), "start"], collapse = ","))

  }
  colnames(genom_positions) <- "genom_positions"
  make_epiread <- cbind(make_epiread, genom_positions)

  make_epiread <- make_epiread %>%
    mutate(num_meth = str_count(make_epiread$meth_pattern, "C"),
           total = str_length(make_epiread$meth_pattern),
           alpha = num_meth/total)

  luad_epireads[[luad_headers[i]]] <- make_epiread
}

# beta value

beta_format <- luad_epireads[1] %>% separate_longer_position(cols = meth_pattern, width = 1) %>%
  group_by(read_id) %>%
  mutate(cpg_pos = index_cg:c(index_cg + n() -1)) %>% ungroup() %>%
  mutate(genom_pos = cg_sites[match(.$cpg_pos, cg_sites$index), "start"])

beta_format <- beta_format %>% group_by(cpg_pos) %>%
  mutate(site_pattern = paste(meth_pattern, collapse = "")) %>% ungroup() %>% distinct(cpg_pos, genom_pos, site_pattern) %>%
  mutate(num_meth = str_count(.$site_pattern, "C"), num_unmeth = str_count(.$site_pattern, "T"))
beta_format <- data.frame(beta_format)

beta_vals <- cg_sites
beta_vals <- beta_vals %>% mutate(meth = beta_format[match(.$index, beta_format$cpg_pos), "num_meth"],
                                  unmeth = beta_format[match(.$index, beta_format$cpg_pos), "num_unmeth"])
beta_vals[is.na(beta_vals)] <- 0

beta_vals$beta <- beta_vals$meth / c(beta_vals$meth + beta_vals$unmeth)
