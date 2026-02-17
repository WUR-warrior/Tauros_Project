# Check and install required packages
required_packages <- c("ggplot2", "readr", "dplyr", "stringr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

# Define the directory where your files are located
data_dir <- "/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/BAM_Stats"  # Current directory, change if needed
# setwd("/path/to/your/data")  # Uncomment and adjust this if needed

# Define the naming pattern for your files
# For example, if your files are named "sample1_coverage.txt", "sample2_coverage.txt", etc.
file_pattern <- "*_coverage_stats.txt$"  # Adjust this pattern to match your files

# Get list of all matching files
coverage_files <- list.files(path = data_dir, pattern = file_pattern, full.names = TRUE)

if(length(coverage_files) == 0) {
  stop("No matching files found. Check your file pattern and directory.")
}

# Print the files that will be processed
cat("Found", length(coverage_files), "files to process:\n")
for(file in coverage_files) {
  cat(" -", basename(file), "\n")
}
cat("\n")

# Create output directory for plots if it doesn't exist
output_dir <- file.path(data_dir, "coverage_plots")
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Container for all chromosome-level data
all_samples_data <- NULL

# Process each file
for(file_path in coverage_files) {
  # Extract sample name from filename
  sample_name <- str_replace(basename(file_path), "_coverage_stats\\.txt$", "")
  if(sample_name == basename(file_path)) {  # If replacement didn't work
    sample_name <- str_replace(basename(file_path), "\\.txt$", "")  # Just remove .txt
  }
  
  cat("Processing", sample_name, "...\n")
  
  # Read the data
  tryCatch({
    coverage_data <- read_tsv(file_path, show_col_types = FALSE)
    
    # Check if expected columns exist
    required_columns <- c("#rname", "meandepth", "%_covered")
    missing_columns <- required_columns[!required_columns %in% names(coverage_data)]
    
    if(length(missing_columns) > 0) {
      cat("  WARNING: Missing required columns:", paste(missing_columns, collapse=", "), "\n")
      cat("  Available columns are:", paste(names(coverage_data), collapse=", "), "\n")
      
      # Try to identify alternative column names
      if("%_covered" %in% missing_columns) {
        pct_col <- names(coverage_data)[grep("coverage|percent|covered", names(coverage_data), ignore.case=TRUE)]
        if(length(pct_col) > 0) {
          cat("  Found potential alternative for %_covered:", pct_col[1], "\n")
          required_columns[required_columns == "%_covered"] <- pct_col[1]
        }
      }
      
      if("#rname" %in% missing_columns) {
        chr_col <- names(coverage_data)[grep("chr|region|contig|scaffold|rname", names(coverage_data), ignore.case=TRUE)]
        if(length(chr_col) > 0) {
          cat("  Found potential alternative for #rname:", chr_col[1], "\n")
          required_columns[required_columns == "#rname"] <- chr_col[1]
        }
      }
      
      if("meandepth" %in% missing_columns) {
        depth_col <- names(coverage_data)[grep("depth|coverage|mean", names(coverage_data), ignore.case=TRUE)]
        if(length(depth_col) > 0) {
          cat("  Found potential alternative for meandepth:", depth_col[1], "\n")
          required_columns[required_columns == "meandepth"] <- depth_col[1]
        }
      }
    }
    
    # Check if we have the data needed for analysis
    if(!(required_columns[1] %in% names(coverage_data) && 
         required_columns[2] %in% names(coverage_data))) {
      cat("  SKIPPING: Cannot find required columns for analysis.\n")
      next
    }
    
    # Add sample column to the data
    coverage_data$sample <- sample_name
    
    # Collect all chromosomal data for combined visualization
    if(is.null(all_samples_data)) {
      all_samples_data <- coverage_data
    } else {
      all_samples_data <- bind_rows(all_samples_data, coverage_data)
    }
    
    # Create a summary of this sample
    tryCatch({
      # Calculate stats with length weighting if available
      if("length" %in% names(coverage_data)) {
        summary_data <- coverage_data %>%
          summarize(
            sample = sample_name,
            total_bases = sum(length, na.rm = TRUE),
            mean_coverage = weighted.mean(!!sym(required_columns[2]), length, na.rm = TRUE),
            median_coverage = median(!!sym(required_columns[2]), na.rm = TRUE),
            min_coverage = min(!!sym(required_columns[2]), na.rm = TRUE),
            max_coverage = max(!!sym(required_columns[2]), na.rm = TRUE)
          )
      } else {
        summary_data <- coverage_data %>%
          summarize(
            sample = sample_name,
            total_regions = n(),
            mean_coverage = mean(!!sym(required_columns[2]), na.rm = TRUE),
            median_coverage = median(!!sym(required_columns[2]), na.rm = TRUE),
            min_coverage = min(!!sym(required_columns[2]), na.rm = TRUE),
            max_coverage = max(!!sym(required_columns[2]), na.rm = TRUE)
          )
      }
      
      # Store for combined summary
      if(!exists("all_summaries")) {
        all_summaries <- summary_data
      } else {
        all_summaries <- bind_rows(all_summaries, summary_data)
      }
      
      cat("  Completed processing", sample_name, "\n")
    }, error = function(e) {
      cat("  ERROR creating summary:", e$message, "\n")
    })
    
  }, error = function(e) {
    cat("  ERROR reading file:", e$message, "\n")
  })
}

# Proceed with combined analysis if data was collected
if(exists("all_summaries") && nrow(all_summaries) > 0) {
  # Sort samples by mean coverage (for better visualization)
  all_summaries <- all_summaries %>% 
    arrange(desc(mean_coverage))
  
  # Save combined summary table
  write.csv(all_summaries, file.path(output_dir, "all_samples_summary.csv"), row.names = FALSE)
  cat("Saved combined statistics to:", file.path(output_dir, "all_samples_summary.csv"), "\n")
  
  # Create publication-quality comparison plots
  
  # 1. Mean Coverage Comparison
  p_mean <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage), y = mean_coverage)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    geom_text(aes(label = round(mean_coverage, 1)), 
              position = position_stack(vjust = 0.94), 
              color = "white", size = 3.5) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "Mean Coverage Depth Comparison",
      x = "Sample",
      y = "Mean Coverage Depth"
    )
  
  # Save plot
  ggsave(file.path(output_dir, "all_samples_mean_coverage.png"), p_mean, 
         width = 10, height = 6, dpi = 300)
  cat("Saved mean coverage comparison plot\n")
  
  # 2. Mean vs Median Coverage
  p_mean_median <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage))) +
    geom_bar(aes(y = mean_coverage, fill = "Mean"), stat = "identity", width = 0.6, position = position_dodge(width = 0.7)) +
    geom_bar(aes(y = median_coverage, fill = "Median"), stat = "identity", width = 0.6, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = c("Mean" = "steelblue", "Median" = "darkgreen"), name = "Metric") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Mean vs Median Coverage Depth",
      x = "Sample",
      y = "Coverage Depth"
    )
  
  # Save plot
  ggsave(file.path(output_dir, "mean_vs_median_coverage.png"), p_mean_median, 
         width = 10, height = 6, dpi = 300)
  cat("Saved mean vs median comparison plot\n")
  
  # 3. Coverage Range (Min to Max)
  p_range <- ggplot(all_summaries, aes(x = reorder(sample, mean_coverage))) +
    geom_linerange(aes(ymin = min_coverage, ymax = max_coverage), 
                   color = "darkgray", size = 1.2) +
    geom_point(aes(y = mean_coverage), color = "red", size = 3) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Coverage Range by Sample",
      x = "Sample",
      y = "Coverage Depth (Min to Max, Red dot = Mean)"
    )
  
  # Save plot
  ggsave(file.path(output_dir, "coverage_range.png"), p_range, 
         width = 10, height = 6, dpi = 300)
  cat("Saved coverage range plot\n")
  
  # Process chromosome-level data if available
  if(!is.null(all_samples_data)) {
    # Only get top samples by coverage for better visualization
    if(length(unique(all_samples_data$sample)) > 10) {
      top_samples <- all_summaries$sample[1:10]  # Top 10 samples
      all_samples_data <- all_samples_data %>%
        filter(sample %in% top_samples)
    }
    
    # Ensure chromosome/region names are consistently ordered
    chrom_col <- required_columns[1]
    all_samples_data[[chrom_col]] <- factor(all_samples_data[[chrom_col]], 
                                            levels = unique(all_samples_data[[chrom_col]]))
    
    # 4. Chromosome-level coverage heatmap
    p_heatmap <- ggplot(all_samples_data, 
                        aes(x = !!sym(chrom_col), y = sample, fill = !!sym(required_columns[2]))) +
      geom_tile() +
      scale_fill_viridis_c(name = "Coverage Depth") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      ) +
      labs(
        title = "Coverage Depth by Chromosome/Region",
        x = "Chromosome/Region",
        y = "Sample"
      )
    
    # Save plot
    ggsave(file.path(output_dir, "chromosome_coverage_heatmap.png"), p_heatmap, 
           width = 12, height = 8, dpi = 300)
    cat("Saved chromosome-level heatmap\n")
  }
  
  cat("\nAll analysis complete! Results saved to:", output_dir, "\n")
} else {
  cat("\nNo files were processed successfully. Check file format and column names.\n")
}

# Define the pattern for depth distribution files
depth_pattern <- "*depth_distribution.txt$"  # Adjust as needed

# List all matching depth distribution files
depth_files <- list.files(path = data_dir, pattern = depth_pattern, full.names = TRUE)

# Process each file
for (file_path in depth_files) {
  sample_name <- str_replace(basename(file_path), "_depth_distribution\\.txt$", "")
  
  cat("Processing:", sample_name, "\n")
  
  tryCatch({
    # Read depth distribution file
    depth_dist <- read.table(file_path, col.names = c("Count", "Depth"))
    
    # Optional: Filter out extreme values to focus on main distribution
    depth_max <- quantile(depth_dist$Depth, 0.99)  # Adjust if needed
    depth_dist_trim <- depth_dist[depth_dist$Depth <= depth_max, ]
    
    # Plot
    p <- ggplot(depth_dist_trim, aes(x = Depth, y = Count)) +
      geom_col(fill = "steelblue", width = 1) +
      scale_y_log10(labels = scales::comma) +  # Log Y for better visibility of large ranges
      theme_bw(base_size = 14) +
      labs(
        title = paste("Read Depth Distribution:", sample_name),
        x = "Read Depth (coverage)",
        y = "Number of Genomic Positions (log scale)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))
      )
    
    # Save plot
    ggsave(filename = file.path(output_dir, paste0(sample_name, "_depth_distribution.png")),
           plot = p, width = 10, height = 6)
    
    cat("✔️  Saved:", sample_name, "\n")
  }, error = function(e) {
    cat("❌ Error in", sample_name, ":", e$message, "\n")
  })
}
