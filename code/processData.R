process_modkit <- function(modkit_fread, path_to_bed_file) {
  reads <- modkit_fread
  reads <- reads$ref_position + 1 # convert to 1 based
  reads <- reads[, ref_position := ifelse(strand == "-", ref_position - 1, ref_position)] # convert to 1 based
  reads <- reads[, read_and_pos := paste0(reads$read_id, ":", reads$ref_position)]
  reads$start <- reads$ref_position
  reads$end <- reads$ref_position

  bed <- fread(path_to_bed_file)

  alpha_and_beta2 <- find_overlaps(as_granges(reads), as_granges(bed)) %>% data.frame() %>%
    filter(call_prob > 0.66) %>% # solution for alpha and beta match
    group_by(read_id) %>%
    .[order(.$alignment_start),] %>%
    mutate(call_code = ifelse(call_code == "m", "M", "U")) %>%
    mutate(meth_pattern = paste(call_code, collapse = "")) %>%
    mutate(num_meth = str_count(meth_pattern, pattern = "M"),
           total = str_length(meth_pattern),
           alpha = num_meth/total) %>%
    mutate(genom_positions = paste0(ref_position, collapse = ",")) %>%
    ungroup()
  alpha_and_beta2 <- alpha_and_beta2 %>% data.frame()
  alpha_and_beta2
}

modkit_to_pat <- function(alpha_and_beta) {
  alpha_and_beta2_pat <- alpha_and_beta2 %>%
    mutate(start_pos = sapply(str_split(genom_positions, ","), head, 1),
           end_pos = sapply(str_split(genom_positions, ","), tail, 1)) %>%
    filter(ref_position == start_pos) %>%
    distinct(read_id, seqnames, start_pos, meth_pattern, alpha, end_pos, beta) %>%
    group_by(seqnames, start_pos, meth_pattern, alpha, end_pos, beta) %>%
    summarise(n=n()) %>%
    ungroup()

  alpha_and_beta2_pat <- alpha_and_beta2_pat %>%
    mutate(start_pos = as.integer(start_pos),
           end_pos = as.integer(end_pos))
  alpha_and_beta2_pat
}

process_wgbs <- function(bismark, bed = NULL) {
  read <- bismark
  names(read) <- c("read_id", "strand", "seqnames", "start", "meth_status")
  read$end <- read$start

  # if (is.null(bed) == TRUE) {
  #    read_overlap <- read
  # } else {
  #  read_overlap <- find_overlaps(as_granges(read), as_granges(bed)) %>% data.frame()
  # }
  read_overlap <- read
  read_overlap[, meth_status := ifelse(meth_status == "Z", "M", "U")]
  read_overlap[, meth_pattern := paste0(meth_status, collapse = ""), by = read_id]
  read_overlap[, num_meth := str_count(meth_pattern, pattern = "M")]
  read_overlap[, total := str_length(meth_pattern)]
  read_overlap[, alpha := num_meth/total]
  read_overlap[, genom_positions := paste0(start, collapse = ","), by = read_id]

  read_overlap <- read_overlap %>% .[order(.$read_id),]

  betas <- betas_for_wgbs(read_overlap)
  betas <- data.frame(betas)
  read_overlap <- read_overlap %>% group_by(read_id) %>%
    mutate(read_start = min(start), read_end = max(end),
           genom_positions = paste0(start, collapse = ",")) %>%
    ungroup()
  read_overlap <- read_overlap %>% data.frame() %>%
    mutate(beta = betas[match(.$start, betas$start), "beta"])

  read_overlap
}

betas_for_wgbs <- function(bismark_processed) {
  res <- bismark_processed %>% group_by(start, meth_status) %>%
    .[order(.$start),] %>%
    mutate(meth_count_pos = n(), unmeth_count_pos = n()) %>%
    mutate(meth_count_pos = ifelse(meth_status == "M", meth_count_pos, NA),
           unmeth_count_pos = ifelse(meth_status == "U", unmeth_count_pos, NA)) %>%
    ungroup() %>%
    distinct(seqnames, start, end, meth_count_pos, unmeth_count_pos) %>%
    group_by(start) %>%
    .[order(.$meth_count_pos),] %>%
    fill(meth_count_pos, .direction = "down") %>%
    fill(unmeth_count_pos, .direction = "up") %>%
    replace(is.na(.), 0) %>%
    mutate(beta = meth_count_pos/(meth_count_pos + unmeth_count_pos)) %>%
    ungroup() %>%
    distinct(seqnames, start, end, meth_count_pos, unmeth_count_pos, beta)

  res
}

pat_for_wgbs <- function(process_res) {
  pat <- process_res %>%
    mutate(start_pos = sapply(str_split(genom_positions, ","), head, 1),
           end_pos = sapply(str_split(genom_positions, ","), tail, 1)) %>%
    filter(start == start_pos) %>%
    distinct(read_id, seqnames, start_pos, meth_pattern, alpha, end_pos, beta) %>%
    group_by(seqnames, start_pos, meth_pattern, alpha, end_pos, beta) %>%
    summarise(n=n()) %>%
    ungroup()  %>%
    mutate(start_pos = as.integer(start_pos),
           end_pos = as.integer(end_pos))

  pat
}

rrbs_plot <- function(long_data, read_ids, cpgs = cg_sites) {
  long_data %>% filter(read_id %in% read_ids) %>%
    distinct(read_id, read_start, read_end, alpha) %>%
    ggplot() +
    geom_segment(aes(x = read_start, xend = read_end, y = alpha, yend = alpha), alpha = 0.2) +
    geom_point(data = filter(long_data, read_id %in% read_ids),
               aes(x = start, y = alpha, colour = meth_status)) +
    geom_point(data = distinct(filter(long_data, read_id %in% read_ids), start, beta),
               aes(x = start, y = beta), shape = 4, size = 3) +
    geom_segment(data = filter(cg_sites, start >= min(filter(long_data, read_id %in% read_ids)$start),
                               end <= max(filter(long_data, read_id %in% read_ids)$end)),
                 aes(x = start, xend = start, y = -0.15, yend = -0.05), colour = "green")
}
