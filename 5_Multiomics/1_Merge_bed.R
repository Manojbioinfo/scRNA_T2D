library(dplyr)
library(readr)

## 0) Input and output folders
bed_dir <- "../Data/Mouse_bed/"
out_dir <- "Mouse_bed_merged/"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

bedfiles <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

## 1) Parse info from each filename
# Example: H3K27ac_WT1_GSM5514720_INS1_250kb.bed
parse_bed_name <- function(path) {
  b <- basename(path)
  p <- strsplit(b, "_")[[1]]
  # p = c("H3K27ac","WT1","GSM5514720","INS1","250kb.bed")
  
  mark   <- p[1]               # H3K27ac, H3K27me3, H3K4me1, H3K4me3
  status <- substr(p[2], 1, 2) # "WT" or "KO" from WT1/KO1
  sample <- substr(p[2], 3, 3) # "1" or "2"
  ins    <- p[4]               # "INS1" or "INS2"
  
  data.frame(mark, status, sample, ins, file = path, stringsAsFactors = FALSE)
}

bedinfo <- do.call(rbind, lapply(bedfiles, parse_bed_name))
bedinfo  # check: ins should be INS1/INS2 now

## 2) Groups with exactly 2 samples (1 and 2) per mark/status/ins
groups <- bedinfo %>%
  group_by(mark, status, ins) %>%
  arrange(sample) %>%
  filter(n() == 2) %>%
  ungroup()

groups   # should now have 2 rows per (mark,status,INS1/2)

## 3) Function: read two BEDs, average score, write merged BED in out_dir
merge_two_bed_avg <- function(f1, f2, out_file) {
  bed1 <- readr::read_tsv(f1, col_names = FALSE, show_col_types = FALSE)
  bed2 <- readr::read_tsv(f2, col_names = FALSE, show_col_types = FALSE)
  
  merged <- dplyr::inner_join(
    bed1, bed2,
    by = c("X1", "X2", "X3"),
    suffix = c(".s1", ".s2")
  ) %>%
    dplyr::mutate(
      # coerce numeric explicitly; here I assume X5 is the score column
      score = (as.numeric(X5.s1) + as.numeric(X5.s2)) / 2
    ) %>%
    dplyr::select(X1, X2, X3, score)
  
  readr::write_tsv(merged, out_file, col_names = FALSE)
  out_file
}

## 4) Loop over all groups and create merged files in out_dir
merged_files <- list()
group_list   <- split(groups, list(groups$mark, groups$status, groups$ins))

for (g in group_list) {
  if (nrow(g) != 2) next
  
  mark   <- g$mark[1]
  status <- g$status[1]
  ins    <- g$ins[1]
  
  f1 <- g$file[g$sample == "1"]
  f2 <- g$file[g$sample == "2"]
  
  # e.g. H3K27ac_WT_INS1_250kb_merged.bed
  out_name <- sprintf("%s_%s_%s_250kb_merged.bed", mark, status, ins)
  out_file <- file.path(out_dir, out_name)
  
  merged_files[[length(merged_files) + 1]] <- merge_two_bed_avg(f1, f2, out_file)
}

merged_files
