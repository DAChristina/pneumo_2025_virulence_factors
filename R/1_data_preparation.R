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
  dplyr::mutate(`ID Pathogen watch` = paste0(`ID Pathogen watch`, ".fasta")) %>% 
  # view() %>% 
  glimpse()


all_fasta <- read.table("raw_data/test_available_fasta.txt", header = F) %>% 
  stats::setNames("all_fasta") %>% 
  dplyr::mutate(availability = "available") %>% 
  glimpse()

test_missing_files <- dplyr::left_join(
  dat_list_all_data %>% 
    mutate(`ID Pathogen watch` = gsub(" ", "_", `ID Pathogen watch`),
           `ID Pathogen watch` = gsub("_IDN_", "_", `ID Pathogen watch`)),
  all_fasta
  ,
  by = c(`ID Pathogen watch` = "all_fasta")
) %>% 
  # filter(is.na(availability)) %>%
  # view() %>% 
  dplyr::mutate(file_location = paste0("/home/ron/pneumo_2025_virulence_factors/raw_data/compiled_all_fasta/", `ID Pathogen watch`),
                ID = gsub(".fasta", "", `ID Pathogen watch`)) %>% 
  dplyr::select(ID, file_location) %>% 
  glimpse()

write.table(test_missing_files, "raw_data/list_prokka_3_6_19F_priority.txt",
            row.names = F, col.names = F, sep = "\t", quote = F)




# virulence profile
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
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" and "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    ) %>% 
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
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
    ) %>% 
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
  ,
  by = c("file_name", "qseqid", "temporary_sseqid")
) %>% 
  dplyr::arrange(dplyr::desc(aa_pident)) %>% 
  dplyr::distinct(file_name, temporary_sseqid, .keep_all = T) %>% 
  dplyr::select(-contains(".y")) %>% 
  dplyr::mutate(
    gene_present = case_when(
      nt_pident >= 80 ~ "present",
      TRUE ~ "absent"
    ),
    protein_function = case_when(
      aa_pident >= 95 & aa_gapopen == 0 ~ "functional",
      aa_pident >= 90 & aa_gapopen <= 5 & aa_mismatch <= 30 ~ "variant",
      aa_pident <= 85 & aa_mismatch >= 25 ~ "possibly defective",
      TRUE ~ "possibly defective"
    ),
    bacwgs_check = case_when(
      gene.x %in% c("cbpD", "cps4A", "cps4B", "cps4D", "hysA", "lytA", "lytC", 
                    "nanA", "nanB", "pavA", "pce/cbpE", "pfbA", "ply", "psaA") ~ "detected",
      TRUE ~ "no"
    )
  ) %>% 
  # view() %>% 
  glimpse()

# nanA is not included in VFDB pro list -_-)
test <- vir_combined_results %>% 
  select(nt_pident, nt_mismatch,
         aa_pident, aa_mismatch, aa_gapopen, gene.x, species.x,
         gene_present, protein_function, bacwgs_check) %>% 
  view() %>% 
  glimpse()