suppressMessages({
  library(tidyverse)
  library(progress)
})

theme_set(
    theme_minimal()
)

options(scipen=15000000)
options(dplyr.summarise.inform = FALSE)

read_data <- function(filename) {
  read_tsv(filename, col_types  = cols()) |>
  filter(FDR < 0.01 & abs(IncLevelDifference) > 0.2) |>
  mutate(event_in_kd = if_else(IncLevelDifference < 0, "Inclusion", "Skipping"))  |>
  mutate(my_id = str_c(upstreamEE,upstreamES,exonStart_0base,exonEnd,downstreamEE,downstreamES,event_in_kd, sep = '_'))
}

my_data_u251 <- read_data('input_skipcount/u251_SE.MATS.JC.txt') |> mutate(assay = 'U251')
my_data_u343 <- read_data('input_skipcount/u343_SE.MATS.JC.txt')  |> mutate(assay = 'U343')
my_data_new <- read_data('input_skipcount/reCTRLvskdDOXY_SE.MATS.JC.txt')  |> mutate(assay = 'reCTRLvskdDOXY')

my_data <- tibble() |>
  bind_rows(my_data_u251) |>
  bind_rows(my_data_u343) |>
  bind_rows(my_data_new) |>
  filter(event_in_kd == "Skipping")

ids_in_two_or_more <-
  my_data |>
  count(my_id) |>
  filter(n > 1) |>
  pull(my_id)

my_data <- my_data |>
  filter(my_id %in% ids_in_two_or_more)

exon_count_data <-
  read_tsv('input/ExonCountData-Table 1.tsv', show_col_types = FALSE) |>
  janitor::clean_names() |>
  mutate(chr_start_stop_gene_name_enst_id = 
    str_replace(chr_start_stop_gene_name_enst_id,
    "([^-]*-[^-]*-[^-]*-[^-]*).*", '\\1'))

my_output <- "all_nolastexon2"
my_clip_data <- 
  exon_count_data |>
  select(-'rna_seq_a', -'rna_seq_b') |>
  select(-'casc3_eif4a3_xl_a', -'casc3_eif4a3_xl_b') |>
  select(-'rnps1_eif4a3_xl_a', -'rnps1_eif4a3_xl_b') |>
  # select(-'casc3_eif4a3_chx_xl_a', -'casc3_eif4a3_chx_xl_b') |>
  # select(-'rnps1_eif4a3_chx_xl_a', -'rnps1_eif4a3_chx_xl_a') |>
  select(-'magoh_eif4a3', -'casc3_eif4a3', -'rnps1_eif4a3',
         -'magoh_casc3', -'magoh_rnps1', -'casc3_eif4a3', -'magoh_casc3',-'magoh_rnps1') |>
  separate(
    chr_start_stop_gene_name_enst_id, 
    into = c('chr','start','stop','gene_name'), sep = '-', convert = TRUE) |>
  distinct(chr,start,stop, .keep_all = TRUE) |>
  pivot_longer(cols = -c('chr','start','stop','gene_name'), names_to = 'assay', values_to = 'counts') |>
  group_by(chr, start, stop, gene_name) |> 
  summarise(counts_sum = sum(counts, na.rm = TRUE)) |>
  ungroup() |>
  rename(counts = 'counts_sum') |>
  filter(counts > 0)
  
my_rnaseq_data <- 
  exon_count_data |>
  select('chr_start_stop_gene_name_enst_id', 'rna_seq_a', 'rna_seq_b') |>
  mutate(rna_seq_a = rna_seq_a + rna_seq_b) |>
  select(-'rna_seq_b') |>
  separate(
    chr_start_stop_gene_name_enst_id, 
    into = c('chr','start','stop','gene_name'), sep = '-', convert = TRUE) |>
  distinct(chr,start,stop, .keep_all = TRUE) |>
  pivot_longer(cols = -c('chr','start','stop','gene_name'), names_to = 'assay', values_to = 'counts') |>
  group_by(chr, start, stop, gene_name) |> 
  summarise(counts_sum = sum(counts, na.rm = TRUE)) |>
  ungroup() |>
  rename(counts_rnaseq = 'counts_sum') |>
  filter(counts_rnaseq > 0) |>
  select(-gene_name)

my_data_wcounts <-
  my_data |>
  mutate(start = exonStart_0base, stop = exonEnd) |>
  mutate(start = if_else(strand == '+', upstreamES, downstreamES )) |>
  mutate(stop = if_else(strand == '+', upstreamEE, downstreamEE)) |>
  left_join(my_clip_data, by = c('chr', 'start', 'stop')) |>
  left_join(my_rnaseq_data, by = c('chr', 'start', 'stop')) |>
  mutate(counts = if_else(is.na(counts), 0, counts)) |>
  mutate(counts_rnaseq = if_else(is.na(counts), 0, counts)) |>
  filter(counts_rnaseq > 0)

sum_counts <-
  my_clip_data |>
  group_by(gene_name) |>
  summarise(sum_counts = sum(counts), n_exons = n())

gene_strand <-
  read_tsv('gencode.v29.annotation.genes.gtf',
         col_names = c(
  "chr","source","feature","start","end","dot","strand","dot2","annot"
         ),
         show_col_types = FALSE) |>
  mutate(gene_name = str_replace(annot, '.*gene_name "([^"]*).*', '\\1')) |>
  select(gene_name, strand)

exon_order <-
  exon_count_data |>
  select(-'rna_seq_a', -'rna_seq_b') |>
  separate(
    chr_start_stop_gene_name_enst_id, 
    into = c('chr','start','stop','gene_name'), sep = '-', convert = TRUE) |>
  distinct(chr,start,stop, .keep_all = TRUE) |>
  select(chr, start, stop, gene_name) |>
  left_join(gene_strand, by = "gene_name") |>
  na.omit() |>
  group_by(gene_name) |>
  mutate(exon_number = 1:n()) |>
  mutate(is_last = if_else(strand == "+", 
    (exon_number == max(exon_number)) | (exon_number ==  max(exon_number) - 1), 
    (exon_number == 1) | (exon_number ==  2))) |>
  mutate(is_first = if_else(strand == "+", 
    exon_number == 1, exon_number == max(exon_number)))


data2plot <-
  my_data_wcounts |>
  select(geneSymbol, gene_name, start, stop, counts, assay, my_id, counts_rnaseq) |>
  arrange(my_id) |>
  distinct(my_id, geneSymbol, gene_name, start, stop, counts, counts_rnaseq) |>
  left_join(sum_counts, by = "gene_name") |>
  filter(n_exons > 3) |>
  left_join(exon_order, by = c("gene_name", "start","stop")) |>
  na.omit() |>
  filter(is_last == FALSE) |>
  mutate(mean_counts = sum_counts/n_exons) |>
  mutate(ratio_counts = counts/mean_counts) |>
  filter(counts_rnaseq > 0) |>
  filter(gene_name %in% my_data_wcounts$gene_name) 

data2sample <-
  my_clip_data |>
  left_join(sum_counts, by = "gene_name") |>
  filter(n_exons > 3) |>
  left_join(exon_order |> select(-chr), by = c("gene_name", "start","stop")) |>
  na.omit() |>
  filter(is_last == FALSE) |>
  mutate(mean_counts = sum_counts/n_exons) |>
  mutate(ratio_counts = counts/mean_counts) |>
  filter(gene_name %in% data2plot$gene_name) |>
  left_join(data2plot |> count(gene_name)) |>
  uncount(n, .id = 'myrepeats') |>
  mutate(gene_name = str_c(gene_name,'_',myrepeats)) |>
  group_by(gene_name) 
  

n_row <- data2plot |> nrow()
median(data2plot |> pull(ratio_counts) |> median())
cat("\n")
cat(n_row)
cat("\n")
cat(data2sample |> distinct(gene_name) |> nrow())
cat("\n")
my_medians <- tibble(i = numeric(), median = numeric())
cat("Bootstraping...\n")
nboots <- 1000
pb <- progress_bar$new(total = nboots,
                       format = "bootstraping [:bar] :current/:total (:percent) eta: :eta",)
for (my_i in seq(1,nboots)) {
  my_medians <- my_medians |> rbind(
    tibble(i = my_i, 
           median = data2sample |>
            sample_n(1) |>
            ungroup() |>
            sample_n(n_row) |>
            pull(ratio_counts) |> median()
            # distinct(gene_name) |> nrow()
    )
  )
  
  pb$tick()
  Sys.sleep(1 / 100)
  # if (my_i %% 100 == 0 ){
  #     cat(my_i)
  #     cat("\n")
  # }
}

p1 <- 
  my_medians |>
  ggplot(aes(x = median)) +
  geom_histogram(binwidth = 0.001) +
  geom_vline(xintercept = median(data2plot |> pull(ratio_counts) |> median()), color = 'red') +
  labs(x = 'Median Ratio', y = 'Frequency (n)')

ggsave(
    filename = str_c(my_output,".pdf"), 
    plot = p1,
    device = 'pdf',
    width = 3.5, 
    height = 2)

my_medians |>
  bind_rows(tibble(i = 0, median = median(data2plot |> pull(ratio_counts) |> median()))) |>
  write_csv(str_c(my_output, ".csv"))
