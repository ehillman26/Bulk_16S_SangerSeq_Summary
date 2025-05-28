#Make sure all input/output paths exist before running the script, or include checks like dir.exists() and file.exists() where needed.
#There are multiple calls to setwd(). This can be brittle if the script is run on another system or if a path is misspelled. Consider using relative paths or here::here().
#Download BLAST files from NCBI to run BLAST on your local machine. Add the full path to the files to run the programs.

#STEP1 - High quality filtering
# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Biostrings", "seqinr"))
install.packages("stringr")
install.packages("tidyverse")
install.packages("dplyr")

library(stringr)
library(tidyverse)
library(Biostrings)
library(seqinr)

getwd()
# Set working directory (Ensure this path exists!)
setwd("/Users/hillman1/Documents/Michigan_16S_Class_Analysis")

# Initialize storage for all HQ sequences
dir.create("HQ_FASTA", showWarnings = FALSE)

# Initialize storage for all HQ sequences
all_HQ_seqs <- DNAStringSet()

# Process each .seq file in the SEQ directory
seq_files <- list.files("SEQ", pattern = "\\.seq$", full.names = TRUE)

for (file in seq_files) {
  seqs <- readDNAStringSet(file)
  
  for (i in seq_along(seqs)) {
    seq_id <- names(seqs)[i]
    seq <- seqs[[i]]
    seq_length <- width(seqs)[i]
    
    cat(seq_id, seq_length, "\n")  # Optional: useful for debugging
    
    if (seq_length >= 50) {
      # Trim based on sequence length
      if (seq_length >= 900) {
        trimmed_seq <- subseq(seq, start = 51, end = 900)
      } else if (seq_length >= 500) {
        trimmed_seq <- subseq(seq, start = 51)
      } else {
        next
      }
      
      # Clean the sequence ID: keep only text before first space
      clean_id <- strsplit(seq_id, " ")[[1]][1]
      
      # Create DNAStringSet with name
      HQ_seq <- DNAStringSet(trimmed_seq)
      names(HQ_seq) <- clean_id
      
      # Save trimmed HQ sequence to file
      output_filename <- file.path("HQ_FASTA", paste0(clean_id, "_HQ.fasta"))
      writeXStringSet(HQ_seq, filepath = output_filename)
      
      # Add to list for merged file
      all_HQ_seqs <- c(all_HQ_seqs, HQ_seq)
    }
  }
}

# Write the merged FASTA to the parent directory
writeXStringSet(all_HQ_seqs, filepath = "merged_HQ_seqs.fasta")


#STEP 3 - BLAST (ensure BLAST executables are in PATH)
blast_db_fasta <- "NCBI_16S_rRNA_gene_DATABASE.fasta"
query_fasta <- "Merged_HQ_Sequences.fasta"
output_file <- "Bio173_Fall23.tsv"

system2("/Users/hillman1/Documents/ncbi-blast-2.16.0+/bin/makeblastdb", args = c("-in", blast_db_fasta, 
                                      "-dbtype", "nucl", 
                                      "-parse_seqids"))

system2("/Users/hillman1/Documents/ncbi-blast-2.16.0+/bin/blastn", args = c("-db", blast_db_fasta,
                                 "-query", query_fasta,
                                 "-evalue", "0.00001",
                                 "-qcov_hsp_perc", "50",
                                 "-outfmt", "6",
                                 "-max_target_seqs", "1",
                                 "-out", output_file))

#STEP 4 - BLAST Summary
blast_tsv <- "Bio173_Fall23.tsv"
fasta_file <- "NCBI_16S_rRNA_gene_DATABASE.fasta"

fasta_seqs <- readDNAStringSet(fasta_file)

# Extract accession and RNA_id from fasta headers
seq_lib <- tibble(
  full_header = names(fasta_seqs)
) %>%
  mutate(
    # Extract accession: match something like NR_123456.1 or NR_123456.x or similar
    accession = str_extract(full_header, "^[A-Z]{1,3}_\\d+\\.\\d+|^[A-Z]{1,3}_\\d+"),
    
    # Extract RNA_id/description by removing accession and any following whitespace
    RNA_id = str_trim(str_replace(full_header, paste0("^", accession, "\\s*"), ""))
  ) %>%
  select(accession, RNA_id)

blast_results <- read_tsv(
  blast_tsv, 
  col_names = c("qseqid", "accession", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), 
  show_col_types = FALSE)

parse_sample <- function(sample_name) {
  sample_name <- str_replace(sample_name, " Sample_Name=", "")
  sample_name <- str_replace_all(sample_name, "[^[:alnum:]_]+", "_")
  parts <- unlist(str_split(sample_name, "[_-]"))
  tibble(sample = parts[1], isolate_id = ifelse(length(parts) > 1, parts[2], NA))
}

sample_info <- map_dfr(blast_results$qseqid, parse_sample)

blast_summary <- blast_results %>%
  bind_cols(sample_info) %>%
  left_join(seq_lib, by = "accession") %>%
  mutate(Best_hit = coalesce(RNA_id, "no hit")) %>%
  select(qseqid, sample, accession, Best_hit, pident, length, mismatch, gapopen,
         qstart, qend, sstart, send, evalue, bitscore)

write_csv(blast_summary, "Bio173_Fall23_SUMMARY.csv")
cat("Bio173_Fall23_SUMMARY.csv file created in", getwd(), "\n")

#STEP 5 - Add HQ DNA sequences to the BLAST summary
library (dplyr)

# Load merged FASTA again (if needed)
#all_HQ_seqs <- readDNAStringSet(merged_fasta_file)
seq_named_vector <- as.character(all_HQ_seqs)
names(seq_named_vector) <- names(all_HQ_seqs)

# Add DNA sequences to BLAST_summary
if (!"qseqid" %in% names(blast_summary)) {
  stop("BLAST_summary must include a 'qseqid' column matching FASTA headers")
}

blast_summary <- blast_summary %>%
  mutate(query_sequence = seq_named_vector[qseqid])

# Write final CSV output
write.csv(blast_summary, "Bio173_Full_Summary.csv", row.names = FALSE)

cat("Pipeline complete. Final summary written to Bio173_Full_Summary.csv\n")

