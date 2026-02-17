required_packages <- c("ggplot2", "readr", "dplyr", "stringr","tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(ggplot2)
library(readr)

data_dir <- "/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/BAM_Stats"

file_pattern <- "*edited_coverage_stats.txt$"

coverage_files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)

if(length(coverage_files) == 0) {
  stop("No matching files found. Check your file pattern and directory.")
}

cat("Found", length(coverage_files), "files to process:\n")
for(file in coverage_files) {
  cat(" -", basename(file), "\n")
}
cat("\n")

output_dir <- file.path(data_dir, "coverage_plots")
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

all_samples_data <- NULL

for(file_path in coverage_files) {
  sample_name <- str_replace(basename(file_path), "_coverage_stats\\.txt$", "")
  if(sample_name == basename(file_path)) {
    sample_name <- str_replace(basename(file_path), "\\.txt$", "")
  }
  
  cat("Processing", sample_name, "...\n")
  
  tryCatch({
    coverage_data <- read_tsv(file_path, show_col_types = FALSE)
    
    required_columns <- c("#rname", "meandepth", "%_covered")
    missing_columns <- required_columns[!required_columns %in% names(coverage_data)]
    
    if(length(missing_columns) > 0) {
      if("%_covered" %in% missing_columns) {
        pct_col <- names(coverage_data)[grep("coverage|percent|covered", names(coverage_data), ignore.case=TRUE)]
        if(length(pct_col) > 0) {
          required_columns[required_columns == "%_covered"] <- pct_col[1]
        }
      }
      if("#rname" %in% missing_columns) {
        chr_col <- names(coverage_data)[grep("chr|region|contig|scaffold|rname", names(coverage_data), ignore.case=TRUE)]
        if(length(chr_col) > 0) {
          required_columns[required_columns == "#rname"] <- chr_col[1]
        }
      }
      if("meandepth" %in% missing_columns) {
        depth_col <- names(coverage_data)[grep("depth|coverage|mean", names(coverage_data), ignore.case=TRUE)]
        if(length(depth_col) > 0) {
          required_columns[required_columns == "meandepth"] <- depth_col[1]
        }
      }
    }
    
    if(!(required_columns[1] %in% names(coverage_data) && 
         required_columns[2] %in% names(coverage_data))) {
      cat("  SKIPPING: Cannot find required columns for analysis.\n")
      next
    }
    
    coverage_data$sample <- sample_name
    if(is.null(all_samples_data)) {
      all_samples_data <- coverage_data
    } else {
      all_samples_data <- bind_rows(all_samples_data, coverage_data)
    }
    
    summary_data <- if("length" %in% names(coverage_data)) {
      coverage_data %>%
        summarize(
          sample = sample_name,
          total_bases = sum(length, na.rm = TRUE),
          mean_coverage = weighted.mean(!!sym(required_columns[2]), length, na.rm = TRUE),
          median_coverage = median(!!sym(required_columns[2]), na.rm = TRUE),
          min_coverage = min(!!sym(required_columns[2]), na.rm = TRUE),
          max_coverage = max(!!sym(required_columns[2]), na.rm = TRUE)
        )
    } else {
      coverage_data %>%
        summarize(
          sample = sample_name,
          total_regions = n(),
          mean_coverage = mean(!!sym(required_columns[2]), na.rm = TRUE),
          median_coverage = median(!!sym(required_columns[2]), na.rm = TRUE),
          min_coverage = min(!!sym(required_columns[2]), na.rm = TRUE),
          max_coverage = max(!!sym(required_columns[2]), na.rm = TRUE)
        )
    }
    
    if(!exists("all_summaries")) {
      all_summaries <- summary_data
    } else {
      all_summaries <- bind_rows(all_summaries, summary_data)
    }
    
    cat("  Completed processing", sample_name, "\n")
  }, error = function(e) {
    cat("  ERROR processing file:", e$message, "\n")
  })
}

if(exists("all_summaries") && nrow(all_summaries) > 0) {
  all_summaries <- all_summaries %>% arrange(desc(mean_coverage))
  write.csv(all_summaries, file.path(output_dir, "all_samples_summary.csv"), row.names = FALSE)
  
  p_mean <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage), y = mean_coverage)) +
    geom_col(fill = "#4B9CD3", width = 0.7) +
    theme_classic(base_size = 12) +
    labs(title = "Mean Coverage Depth", x = "Sample", y = "Mean Depth") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "all_samples_mean_coverage.png"), p_mean, width = 9, height = 5, dpi = 300)
  
  p_mean_median <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage))) +
    geom_col(aes(y = mean_coverage, fill = "Mean"), width = 0.4, position = position_dodge(width = 0.5)) +
    geom_col(aes(y = median_coverage, fill = "Median"), width = 0.4, position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = c("Mean" = "#4B9CD3", "Median" = "#E87461")) +
    theme_classic(base_size = 12) +
    labs(title = "Mean vs Median Coverage", x = "Sample", y = "Coverage", fill = "Metric") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "mean_vs_median_coverage.png"), p_mean_median, width = 9, height = 5, dpi = 300)
  
  p_range <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage))) +
    geom_linerange(aes(ymin = min_coverage, ymax = max_coverage), color = "grey50", size = 1) +
    geom_point(aes(y = mean_coverage), color = "#D7263D", size = 2) +
    theme_classic(base_size = 12) +
    labs(title = "Coverage Range per Sample", x = "Sample", y = "Coverage Depth") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "coverage_range.png"), p_range, width = 9, height = 5, dpi = 300)
  
  if(!is.null(all_samples_data)) {
    if(length(unique(all_samples_data$sample)) > 10) {
      top_samples <- all_summaries$sample[1:10]
      all_samples_data <- filter(all_samples_data, sample %in% top_samples)
    }
    
    chrom_col <- required_columns[1]
    all_samples_data[[chrom_col]] <- factor(all_samples_data[[chrom_col]], 
                                            levels = unique(all_samples_data[[chrom_col]]))
    
    p_heatmap <- ggplot(all_samples_data, 
                        aes(x = !!sym(chrom_col), y = sample, fill = !!sym(required_columns[2]))) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "white", high = "#3E8E7E", name = "Depth") +
      theme_minimal(base_size = 11) +
      labs(title = "Coverage Depth by Region", x = "Chromosome / Contig", y = "Sample") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            plot.title = element_text(hjust = 0.5),
            legend.position = "right")
    ggsave(file.path(output_dir, "chromosome_coverage_heatmap.png"), p_heatmap, width = 12, height = 6, dpi = 300)
  }
  
  cat("\nAll analysis complete! Results saved to:", output_dir, "\n")
} else {
  cat("\nNo files were processed successfully. Check file format and column names.\n")
}
