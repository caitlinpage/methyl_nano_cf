library(plyranges)
library(dplyr)
library(ggplot2)

# plot a region for 1 sample

## find high density region
epiread_format %>% mutate(row_num = 1:n()) %>% .[1:1000,] %>%
  ggplot(aes(x = row_num, y = gap_to_next_read)) +
  geom_point()

epiread_format[110:125,] %>%
  ggplot(aes(x = start, y = alpha, colour = read_id)) +
  geom_segment(aes(x = start, xend = end, y = alpha, yend = alpha))


# Plot a gene promoter region across all samples

pair_overlaps(as_granges(overlap_bins_reads), as_granges(ensdb_genes_filt_prom)) %>% data.frame() %>% dplyr::filter(gene_name == "MGAT3") %>%
  ggplot(aes(x = read_start, y = alpha, colour = sample, group = sample)) +
  geom_segment(aes(x = read_start, xend = read_end, y = alpha, yend = alpha)) +
  labs(title = "MGAT3") +
  theme(legend.position = "none") +
  facet_wrap(~type, nrow = 2)

## by cancer group

pair_overlaps(as_granges(overlap_bins_reads), as_granges(ensdb_genes_filt_prom)) %>% data.frame() %>% dplyr::filter(gene_name == "CYTH4") %>%
  ggplot(aes(x = read_start, y = alpha, colour = sample, group = sample)) +
  geom_segment(aes(x = read_start, xend = read_end, y = alpha, yend = alpha), linewidth = 1) +
  labs(title = "CYTH4") +
  facet_wrap(~cancer_group, nrow = 3) +
  theme_classic() +
  theme(legend.position = "none")

# 1 position per read joined by a line to see the trends more easily

pair_overlaps(as_granges(overlap_bins_reads), as_granges(ensdb_genes_filt_prom)) %>% data.frame() %>% dplyr::filter(gene_name == "CYTH4") %>%
  ggplot(aes(x = read_start, y = alpha, colour = sample, group = sample)) +
  geom_line() +
  labs(title = "CYTH4") +
  facet_wrap(~cancer_group, nrow = 3) +
  theme_classic() +
  theme(legend.position = "none")
