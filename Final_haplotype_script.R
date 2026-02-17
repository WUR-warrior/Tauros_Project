library(tidyverse)
library(data.table)
library(viridis)

output_dir <- "C:/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/Results/Haplotype Analysis/No LD"
input_file <- file.path(output_dir, "Tauros_Breeds/haplotype_matches_500.csv")
aurochs_file <- file.path(output_dir, "Tauros_Aurochs/SAMPLE_vs_AURG_zscores.csv")

z_df <- fread(input_file)
setnames(z_df, make.names(names(z_df)))
z_df[, c("CHR", "Start", "End") := tstrsplit(Region, ":|-")]
z_df[, `:=`(
  CHR = as.character(CHR),
  Start = as.numeric(Start),
  End = as.numeric(End)
)]

breed_z_threshold <- 1.1
aurochs_z_threshold <- 0.8

sample_cols <- setdiff(names(z_df), c("Region", "CHR", "Start", "End"))
z_long <- melt(z_df, id.vars = c("Region", "CHR", "Start", "End"),
               measure.vars = sample_cols,
               variable.name = "Sample", value.name = "Zscore")
z_long <- z_long[!is.na(Zscore)]
z_long[, Bin5Mb := floor(Start / 5e6) * 5e6]
z_long[, Breed := str_extract(Sample, "^[A-Za-z]+")]

breed_bin_z <- z_long[, .(BreedZ = mean(Zscore)), by = .(CHR, Bin5Mb, Breed)]
breed_bin_z[, BinMid := Bin5Mb + 2.5e6]
breed_bin_z[, CHR := factor(CHR, levels = as.character(1:29), labels = paste0("Chr", 1:29))]

top_breed_bins <- breed_bin_z[BreedZ >= breed_z_threshold]

p1 <- ggplot(top_breed_bins, aes(x = BinMid / 1e6, y = Breed, fill = BreedZ)) +
  geom_tile(width = 5, height = 0.9) +
  facet_wrap(~CHR, ncol = 5) +
  scale_fill_viridis(option = "plasma", name = "Z-score", na.value = "grey90") +
  labs(
    title = "Top Breed Bins (Filtered)",
    subtitle = "Per-bin mean Z-scores for all breeds; only bins with Z > 1.1 shown",
    x = "Genomic Position (Mb)",
    y = "Breed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12, margin = margin(t = 12)),
    axis.text.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  )

ggsave(file.path(output_dir, "top_shared_bins_filtered_simplified.pdf"), p1, width = 22, height = 14)

#CSV
top_breed_bins_wide <- dcast(top_breed_bins, CHR + Bin5Mb ~ Breed, value.var = "BreedZ")

top_breed_summary <- top_breed_bins[, .(
  Breeds = paste(Breed, collapse = ", ")
), by = .(CHR, Bin5Mb)]

collapsed_top_breed_bins <- merge(top_breed_summary, top_breed_bins_wide, by = c("CHR", "Bin5Mb"))

collapsed_top_breed_bins[, `:=`(
  Start = as.character(Bin5Mb),
  End = as.character(Bin5Mb + 5e6),
  Genomic_Location = paste0(CHR, ":", Bin5Mb, "-", Bin5Mb + 5e6)
)]

setcolorder(
  collapsed_top_breed_bins,
  c("CHR", "Genomic_Location", "Breeds", "LM", "MN", "PO", "SA")
)

options(scipen = 999)

fwrite(
  collapsed_top_breed_bins,
  file.path(output_dir, "top_breed_bins_only_collapsed.csv"),
  na = ""
)

p2 <- ggplot(breed_bin_z, aes(x = BinMid / 1e6, y = Breed, fill = BreedZ)) +
  geom_tile(width = 5, height = 0.9) +
  facet_wrap(~CHR, ncol = 5) +
  scale_fill_viridis(option = "plasma", name = "Z-score", na.value = "grey90") +
  labs(
    title = "All Breed Bins (Unfiltered)",
    subtitle = "Per-bin mean Z-scores for all breeds; no threshold filtering",
    x = "Genomic Position (Mb)",
    y = "Breed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12, margin = margin(t = 12)),
    axis.text.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  )

ggsave(file.path(output_dir, "all_breed_bins_unfiltered_simplified.pdf"), p2, width = 22, height = 14)

p3 <- ggplot(breed_bin_z, aes(x = BinMid / 1e6, y = BreedZ, color = Breed)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~CHR, scales = "free_x", ncol = 5) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Breed Z-scores per Genomic Bin",
    subtitle = "Per-bin mean Z-scores across all individuals per breed; reflects similarity to Tauros",
    x = "Genomic Position (Mb)",
    y = "Z-score"
  ) +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "breed_zscore_lineplot.pdf"), p3, width = 18, height = 10)

aurochs <- fread(aurochs_file)
aurochs[, c("CHR", "Start", "End") := tstrsplit(gsub(":", "-", Region), "-")]
aurochs[, `:=`(
  CHR = as.character(CHR),
  Start = as.numeric(Start),
  End = as.numeric(End),
  Zscore = as.numeric(AURG)
)]
aurochs <- aurochs[!is.na(Zscore)]
aurochs[, Bin5Mb := floor(Start / 5e6) * 5e6]
aurochs_summary <- aurochs[, .(MeanZ = mean(Zscore)), by = .(CHR, Bin5Mb)]
aurochs_summary[, BinMid := Bin5Mb + 2.5e6]
aurochs_summary[, CHR := factor(CHR, levels = as.character(1:29), labels = paste0("Chr", 1:29))]

p4 <- ggplot(aurochs_summary, aes(x = BinMid / 1e6, y = MeanZ)) +
  geom_line(color = "darkred") +
  facet_wrap(~CHR, scales = "free_x", ncol = 5) +
  labs(
    title = "Aurochs Z-scores per Genomic Bin",
    subtitle = "Per-bin mean Z-scores reflecting similarity of Aurochs to Tauros",
    x = "Genomic Position (Mb)",
    y = "Z-score"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(output_dir, "aurochs_zscore_lineplot.pdf"), p4, width = 18, height = 10)

shared_bins <- merge(
  breed_bin_z[BreedZ >= breed_z_threshold],
  aurochs_summary[MeanZ >= aurochs_z_threshold],
  by = c("CHR", "Bin5Mb")
)

p5 <- ggplot(shared_bins, aes(x = (Bin5Mb + 2.5e6) / 1e6, y = Breed, fill = BreedZ)) +
  geom_tile(width = 5, height = 0.9) +
  facet_wrap(~CHR, scales = "fixed", ncol = 5, strip.position = "top") +
  scale_fill_viridis(option = "magma", name = "Breed Z-score") +
  labs(
    title = "Shared High Z-score Bins (Breed Z ≥ 1.1 and Aurochs Z ≥ 0.8)",
    subtitle = "Per-bin mean Z-scores shared between Aurochs and Breeds",
    x = "Genomic Position (Mb)",
    y = "Breed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.placement = "outside",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.x = element_text(margin = margin(t = 6)),
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.y = element_text(margin = margin(r = 8))
  )


ggsave(file.path(output_dir, "shared_high_bins_tile_plot.pdf"), p5, width = 18, height = 10)


breed_z_wide <- dcast(shared_bins, CHR + Bin5Mb ~ Breed, value.var = "BreedZ")

aurochs_z <- shared_bins[, .(Aurochs_Z = round(unique(MeanZ), 6)), by = .(CHR, Bin5Mb)]

shared_bins_export <- merge(aurochs_z, breed_z_wide, by = c("CHR", "Bin5Mb"))

shared_bins_export[, Genomic_Location := paste0("Chr", CHR, ":", Bin5Mb, "-", Bin5Mb + 5e6)]

shared_bins_export[, CHR := paste0("Chr", CHR)]

setcolorder(shared_bins_export, c("CHR", "Genomic_Location", "Aurochs_Z"))

fwrite(shared_bins_export, file.path(output_dir, "shared_high_bins_plot5.csv"), quote = FALSE, sep = ",", na = "")

library(magick)

plot_pdf_files <- c(
  "top_breed_bins_only_tile_plot.pdf",
  "all_breed_bins_unfiltered_tile_plot.pdf",
  "breed_zscore_lineplot.pdf",
  "aurochs_zscore_lineplot.pdf",
  "shared_high_bins_tile_plot.pdf"
)

for (pdf_name in plot_pdf_files) {
  pdf_path <- file.path(output_dir, pdf_name)
  png_path <- sub("\\.pdf$", ".png", pdf_path)
  if (file.exists(pdf_path)) {
    message("Converting to PNG: ", basename(pdf_path))
    img <- image_read_pdf(pdf_path, density = 300)
    image_write(img, path = png_path, format = "png")
  } else {
    warning("File not found: ", pdf_path)
  }
}

library(dplyr)
library(readr)

df <- read_csv("C:/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/Results/Haplotype Analysis/No LD/top_breed_bins_only_collapsed.csv")

options(scipen = 999)

df_filtered <- df %>%
  filter(if_any(c(LM, MN, PO, SA), ~ . >= 1.1)) %>%
  mutate(
    Start = as.character(Start),
    End = as.character(End),
    Genomic_Location = paste0(CHR, ":", Start, "-", End)
  ) %>%
  select(CHR, Genomic_Location, LM, MN, PO, SA)

write_csv(df_filtered, "C:/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/Results/Haplotype Analysis/No LD/top_breed_bins_simplified.csv", na="")


