library(ggplot2)
library(data.table)
library(dplyr)

Q_file <- "/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/PLINK_pruned_final_relaxed_no_ancestral.6.Q"
fam_file <- "/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/PLINK_pruned_final_relaxed_no_ancestral.fam"
K <- 6

Q <- fread(Q_file, header = FALSE)
fam <- fread(fam_file)
stopifnot(nrow(Q) == nrow(fam))  # Safety check

colnames(Q)[1:K] <- paste0("Cluster", 1:K)

Q$IID <- fam$V2
Q$population <- fam$V1

Q_long <- melt(Q, id.vars = c("IID", "population"), 
               variable.name = "Ancestry", value.name = "Proportion")

Q_long <- Q_long %>%
  group_by(IID) %>%
  mutate(ordering = Proportion[Ancestry == "Cluster1"]) %>%
  ungroup()

ordered_ids <- Q_long %>%
  distinct(IID, population, ordering) %>%
  arrange(population, -ordering) %>%
  pull(IID)

Q_long$IID <- factor(Q_long$IID, levels = ordered_ids)

ggplot(Q_long, aes(x = IID, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ population, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "top"
  ) +
  labs(x = "Individuals", y = "Ancestry Proportion", fill = paste0("Ancestry (K=", K, ")"))

