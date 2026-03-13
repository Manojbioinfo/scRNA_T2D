library(dplyr)
library(GenomicRanges)

# -------------------------------
# Parameters
# -------------------------------
cis_ext   <- 250000     # ±250 kb
trans_ext <- 1000000    # ±1 Mb

# -------------------------------
# Helper function: annotate cis/trans
# -------------------------------
annotate_cis_trans <- function(region_df, gene_coords_df,
                               cis_ext = 250000,
                               trans_ext = 1000000) {
  
  chr <- as.character(region_df$seqnames)
  mid <- region_df$start + floor((region_df$end - region_df$start) / 2)
  
  gene_coords_df |>
    mutate(
      cis_trans = case_when(
        # CIS: same chr, ±cis_ext
        chr == !!chr &
          start <= (mid + cis_ext) &
          end   >= (mid - cis_ext) ~ "cis",
        
        # TRANS: same chr, ±trans_ext but outside cis
        chr == !!chr &
          start <= (mid + trans_ext) &
          end   >= (mid - trans_ext) &
          !(start <= (mid + cis_ext) & end >= (mid - cis_ext)) ~ "trans",
        
        # TRANS: different chromosome
        chr != !!chr ~ "trans",
        
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(cis_trans))
}

# -------------------------------
# INS1
# -------------------------------
ins1_region <- as.data.frame(GRanges(
  seqnames = "chr19",
  ranges = IRanges(
    start = 52264297,
    end   = 52265015
  )
))

ins1_annot <- annotate_cis_trans(ins1_region, gene_coords_df,
                                 cis_ext = cis_ext,
                                 trans_ext = trans_ext)

# Save cis only
ins1_cis <- dplyr::filter(ins1_annot, cis_trans == "cis")
saveRDS(ins1_cis, "INS1_cis_genes.rds")

# Save trans only
ins1_trans <- dplyr::filter(ins1_annot, cis_trans == "trans")
saveRDS(ins1_trans, "INS1_trans_genes.rds")

# -------------------------------
# INS2
# -------------------------------
ins2_region <- as.data.frame(GRanges(
  seqnames = "chr7",
  ranges = IRanges(
    start = 142678656,
    end   = 142679726
  )
))

ins2_annot <- annotate_cis_trans(ins2_region, gene_coords_df,
                                 cis_ext = cis_ext,
                                 trans_ext = trans_ext)

# Save cis only
ins2_cis <- dplyr::filter(ins2_annot, cis_trans == "cis")
saveRDS(ins2_cis, "INS2_cis_genes.rds")

# Save trans only
ins2_trans <- dplyr::filter(ins2_annot, cis_trans == "trans")
saveRDS(ins2_trans, "INS2_trans_genes.rds")
