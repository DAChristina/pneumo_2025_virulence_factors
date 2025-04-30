# test
# need to be loaded:
# 1. metadata
# 2. blastn result and blastx result

# I test panvita instead

library(tidyverse)

# load metadata
dat_list_all_data <- dplyr::bind_rows(
  readxl::read_excel("raw_data/serotype 19F-003 and 006.xlsx", sheet = "Serotype 3") %>% 
    select(`ID Pathogen watch`, serotype, `Source of sample`),
  readxl::read_excel("raw_data/serotype 19F-003 and 006.xlsx", sheet = "Serogroup 6") %>% 
    select(`ID Pathogen watch`, serotype, `Source of sample`),
  readxl::read_excel("raw_data/serotype 19F-003 and 006.xlsx", sheet = "19F") %>% 
    select(`ID Pathogen watch`, serotype, `Source of sample`)
) %>% 
  dplyr::mutate(fasta_name = paste0(`ID Pathogen watch`, ".fasta"),
                fasta_name = gsub(" ", "_", fasta_name),
                fasta_name = gsub("_IDN_", "_", fasta_name),
                id_cleaned = gsub(".fasta", "", fasta_name)) %>% 
  dplyr::rename(id_ori = `ID Pathogen watch`,
                source = `Source of sample`) %>% 
  # view() %>%
  glimpse()


all_fasta <- read.table("raw_data/test_available_fasta.txt", header = F) %>% 
  stats::setNames("fasta_name") %>% 
  dplyr::mutate(availability = "available") %>% 
  glimpse()

test_missing_files <- dplyr::left_join(
  dat_list_all_data
  ,
  all_fasta
  ,
  by = "fasta_name"
) %>% 
  # filter(is.na(availability)) %>%
  # view() %>% 
  dplyr::mutate(file_location = paste0("/home/ron/pneumo_2025_virulence_factors/raw_data/compiled_all_fasta/", fasta_name)) %>% 
  dplyr::select(id_cleaned, file_location) %>% 
  glimpse()

# generate qfile
write.table(test_missing_files, "raw_data/list_prokka_3_6_19F_priority.txt",
            row.names = F, col.names = F, sep = "\t", quote = F)

# virulence profile from VFDB ##################################################
vir_combined_results <- dplyr::left_join(
  read.table("outputs/result_blast_virulences/blastn_tabular_setA_nt_VFDB.txt",
             header = F, sep = "\t") %>% 
    stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                      "nt_pident", "nt_length", "nt_mismatch",
                      "nt_gapopen", "nt_qstart", "nt_qend",
                      "nt_sstart", "nt_send",
                      "nt_evalue", "nt_bitscore")) %>% 
    dplyr::left_join(
      read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_nt.txt",
                 header = F, sep = "\t") %>% 
        dplyr::rename(header = V1) %>% 
        dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                      gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                      gene_class = str_extract(header, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    ) %>% 
    dplyr::mutate(
      # temporary_sseqid = case_when(
      #   temporary_sseqid == "VFG005653" ~ "VFG005653(gb|WP_142355754.1)", # different ID for nanA from VFDB (gene) and NCBI (amino acids)
      #   TRUE ~ temporary_sseqid
      # ),
      # nt_sdiff = abs(nt_sstart-nt_send),
      # nt_qdiff = abs(nt_qstart-nt_qend),
      nt_lengthdiff = abs(abs(nt_sstart-nt_send) - abs(nt_qstart-nt_qend))
    )
  ,
  read.table("outputs/result_blast_virulences/blastx_tabular_setA_aa_VFDB.txt",
             header = F, sep = "\t") %>% 
    stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                      "aa_pident", "aa_length", "aa_mismatch",
                      "aa_gapopen", "aa_qstart", "aa_qend",
                      "aa_sstart", "aa_send",
                      "aa_evalue", "aa_bitscore")) %>% 
    dplyr::left_join(
      read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_pro.txt",
                 header = F, sep = "\t") %>% 
        dplyr::rename(header = V1) %>% 
        dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                      gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    )
  ,
  by = c("file_name", "qseqid", "temporary_sseqid")
) %>% 
  dplyr::arrange(dplyr::desc(aa_pident)) %>% 
  dplyr::distinct(file_name, temporary_sseqid, .keep_all = T) %>% 
  dplyr::select(-contains(".y")) %>% 
  dplyr::mutate(
    gene_present = case_when(
      nt_pident >= 70 ~ "present",
      TRUE ~ "absent"
    ),
    aa_lengthdiff = abs(abs(aa_sstart-aa_send) - abs(aa_qstart-aa_qend)),
    # protein_function = case_when(
    #   aa_pident >= 95 & aa_gapopen == 0 ~ "functional",
    #   aa_pident >= 90 & aa_gapopen <= 5 & aa_mismatch <= 30 ~ "variant",
    #   aa_pident <= 85 & aa_mismatch >= 25 ~ "possibly defective",
    #   TRUE ~ "possibly defective"
    # ),
  ) %>% 
  # view() %>% 
  glimpse()

# nanA is not included in VFDB pro list -_-)
test <- vir_combined_results %>% 
  select(file_name,
         nt_pident, nt_length, nt_lengthdiff,
         nt_mismatch, aa_lengthdiff,
         aa_pident, aa_mismatch, aa_gapopen, gene.x, species.x, gene_class,
         gene_present, #protein_function
         ) %>% 
  # view() %>% 
  glimpse()

# test analysis from panvita
# fuzzyjoin requires weird & needy data preparation & arrangement
vir_mtx <- read.csv("outputs/result_panvita/Results_vfdb_25-04-2025_09-51-14/matriz_vfdb.csv", sep = ";") %>% 
  dplyr::select(-X) %>% 
  tidyr::pivot_longer(cols = 2:ncol(.),
                      names_to = "gene",
                      values_to = "aa_percent") %>% 
  dplyr::left_join(
    dat_list_all_data
    ,
    by = c("Strains" = "id_cleaned")
  ) %>% 
  dplyr::mutate(gene = gene %>%
      str_replace("srtC", "srtC-") %>%
        str_replace("-.", "-") %>%
      str_replace("\\.", "/")
  ) %>% 
  glimpse()

# apparently R cannot read huge samples from headers. I filter out headers on linux instead
headers <- read.table("outputs/result_panvita/VFDB_setA_compiled_headers_pro_streptococcus_only.txt",
                      header = F, sep = "\t") %>%
  # omit VFDB original data for Streptococcus genus
  # read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_nt.txt",
  #                     header = F, sep = "\t") %>%
  dplyr::rename(header = V1) %>% 
  dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                gene_class = str_extract(header, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
                species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
  ) %>% 
  # view() %>% 
  glimpse()

headers_non_streptococcal <- read.table("outputs/result_panvita/VFDB_setA_compiled_headers_pro.txt",
                                        header = F, sep = "\t") %>%
  # omit VFDB original data for Streptococcus genus
  # read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_nt.txt",
  #                     header = F, sep = "\t") %>%
  dplyr::rename(header = V1) %>% 
  dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                gene_class = str_extract(header, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
                species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
  ) %>% 
  # view() %>% 
  glimpse()

df_joined <- fuzzyjoin::regex_left_join(
  vir_mtx, headers,
  by = c("gene" = "gene")
) %>% 
  dplyr::mutate(
    gene_class = case_when(
      str_detect(gene.x, "pav|pfb|srt|pce|cbp") ~ "Adherence",
      str_detect(gene.x, "cps") ~ "Immune modulation",
      str_detect(gene.x, "gnd") ~ "Nutritional/Metabolic factor",
      TRUE ~ gene_class
      )
    ) %>% 
  dplyr::mutate(source = factor(source),
                serotype = case_when(
                  serotype == "06E(6B)" ~ "06B",
                  serotype == "003" ~ "03",
                  TRUE ~ serotype
                ),
                serotype = factor(serotype),
                gene.x = case_when(
                  is.na(species) ~ paste0(gene.x, " (non-streptococcal)"),
                  TRUE ~ gene.x
                  ),
                species = case_when(
                  gene.x == "gndA (non-streptococcal)" ~ "Klebsiella pneumoniae subsp. pneumoniae NTUH-K2044",
                  gene.x == "cpsI (non-streptococcal)" ~ "Enterococcus faecalis V583",
                  TRUE ~ species
                  )
                ) %>% 
  # view() %>% 
  glimpse()

filtered_report <- df_joined %>% 
  dplyr::rename(gene = gene.x) %>% 
  dplyr::select(-c("header", "gene.y", "gene_id")) %>% 
  glimpse()

write.csv(df_joined, "outputs/result_panvita/panvita_results_compiled.csv", row.names = F)
write.csv(filtered_report, "report/panvita_results_compiled_report.csv", row.names = F)

test_missing_geneClass <- df_joined %>% 
  dplyr::filter(is.na(gene_class)) %>% 
  dplyr::select(Strains, gene.x, gene_class,
                -c("gene_id")) %>% 
  # view() %>% 
  glimpse()

test_empty_species <- df_joined %>% 
  dplyr::filter(is.na(species)) %>% 
  dplyr::select(Strains, gene.x, gene_class, header, species,
                -c("gene_id")) %>% 
  # view() %>%
  glimpse()
unique(test_empty_species$gene.x)


# test visualisation
# faceted by gene
png("pictures/panvita_aa_percent_1allgenes.png",
    width = 20, height = 20, units = "cm", res = 300)
ggplot(df_joined,
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()

# faceted by gene_class
png("pictures/panvita_aa_percent_2geneClass.png",
    width = 20, height = 15, units = "cm", res = 300)
ggplot(df_joined,
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  facet_wrap(~ gene_class, drop = F) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()


################################################################################
# filtered by source only
png("pictures/panvita_aa_percent_3a_adherence.png",
    width = 20, height = 15, units = "cm", res = 300)
ggplot(df_joined %>% 
         dplyr::filter(gene_class == "Adherence"),
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  labs(title = "Adherence") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()

png("pictures/panvita_aa_percent_3b_exoenzyme.png",
    width = 20, height = 8, units = "cm", res = 300)
ggplot(df_joined %>% 
         dplyr::filter(gene_class == "Exoenzyme"),
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  labs(title = "Exoenzyme") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()

png("pictures/panvita_aa_percent_3c_exotoxin.png",
    width = 10, height = 8, units = "cm", res = 300)
ggplot(df_joined %>% 
         dplyr::filter(gene_class == "Exotoxin"),
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  labs(title = "Exotoxin") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()

png("pictures/panvita_aa_percent_3d_exoenzyme.png",
    width = 20, height = 15, units = "cm", res = 300)
ggplot(df_joined %>% 
         dplyr::filter(gene_class == "Immune modulation"),
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  labs(title = "Immune modulation") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()

png("pictures/panvita_aa_percent_3e_nutritional_metabolicFactor.png",
    width = 20, height = 10, units = "cm", res = 300)
ggplot(df_joined %>% 
         dplyr::filter(gene_class == "Nutritional/Metabolic factor"),
       aes(x = serotype, y = aa_percent, fill = source, colour = source)) +
  geom_violin(trim = TRUE, alpha = 0.4, position = position_dodge(width = 0.8)) +
  # geom_boxplot(alpha = 0.4, aes(fill = source)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1) +
  scale_y_continuous(limits = c(0, 100)) +
  facet_wrap(~ gene.x, drop = F) +
  theme_bw() +
  labs(title = "Nutritional/Metabolic factor") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")
dev.off()
