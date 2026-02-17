pca <- read.table("/Users/urban_jkr2clv/OneDrive/Documents/Semester_13_Wag/Results/Plink and Admixture/Plink_with_LD_trimming/PLINK_PCA_relaxed.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca)-2)))

pca$breed <- gsub("([A-Z]+)([0-9]+)?", "\\1", pca$IID)

pca$breed[pca$breed == "SAMPLE"] <- "Tauros"

pca$breed[pca$breed == "AURG"] <- "Aurochs"

unique_breeds <- unique(pca$breed)

breed_colors <- c("MN" = "red", "MA" = "blue", "PA" = "green", 
                  "SA" = "purple", "PO" = "orange", "LM" = "brown",
                  "Aurochs" = "pink", "Tauros" = "black")

colors <- breed_colors[unique_breeds]

pca$color <- colors[pca$breed]

plot(pca$PC1, pca$PC2,
     xlab="PC1",
     ylab="PC2",
     main="PCA Plot Colored by Breed",
     pch=19, 
     col=pca$color,
     cex=1.2)

pca$display_label <- pca$IID
pca$display_label[pca$display_label == "SAMPLE"] <- "Tauros"
pca$display_label[pca$display_label == "AURG"] <- "Aurochs"

text(pca$PC1, pca$PC2, labels=pca$display_label, pos=4, cex=0.6, col="black")

legend("bottomright", 
       legend=unique_breeds, 
       col=colors[unique_breeds], 
       pch=19,  
       title="Breeds",
       cex=0.6,
       pt.cex=0.8,
       ncol=2, 
       bg="white",
       box.col="black",
       box.lwd=1)

print("Breed assignments:")
print(table(pca$breed))
