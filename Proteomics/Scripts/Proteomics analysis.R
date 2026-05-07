# 1. Install packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "limma",
  "clusterProfiler",
  "org.Mm.eg.db",
  "WGCNA"
)

cran_pkgs <- c(
  "openxlsx",
  "readxl",
  "writexl",
  "tibble",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "broom",
  "car",
  "ggplot2",
  "ggrepel",
  "scales",
  "igraph",
  "ggraph",
  "tidygraph",
  "UpSetR",
  "doParallel"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = TRUE)
install.packages(cran_pkgs, dependencies = TRUE)

# 2. Load packages ----
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(WGCNA)

library(openxlsx)
library(readxl)
library(writexl)

library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)
library(car)

library(ggplot2)
library(ggrepel)
library(scales)

library(igraph)
library(ggraph)
library(tidygraph)
library(UpSetR)

library(doParallel)
# 3. Cleaning up excel file and making data frame into a matrix----
proteomics <- read_excel("Proteomics_v2.xlsx")
# Clean and prepare the data
proteomics <- proteomics[, -(1:3)]                    
proteomics <- proteomics[, -c(17)]  

# Convert tibble to base R data frame
proteomics <- as.data.frame(proteomics)  

# Set row names using the 'Proteins' column
rownames(proteomics) <- toupper(make.unique(as.character(proteomics$Proteins)))

# Remove the 'protein' column
proteomics$Proteins <- NULL

# Convert to matrix
proteomics_matrix <- as.matrix(proteomics)

# Convert all 0 values to NA
proteomics_matrix[proteomics_matrix == 0] <- NA

# 4. Histogram of raw values----
num_cols <- ncol(proteomics_matrix)
layout_dims <- ceiling(sqrt(num_cols))
par(mfrow = c(layout_dims, ceiling(num_cols / layout_dims)))

for (i in 1:num_cols) {
  hist(proteomics_matrix[, i],
       main = paste(colnames(proteomics_matrix)[i]),
       xlab = "Protein Values",
       ylab = "Frequency",
       col = "lightblue",
       breaks = 30)  # You can adjust number of bins here
}
par(mfrow = c(1,1))

# 5. Transforming data ----

# Log2 transformation
log2_mat <- log2(proteomics_matrix)

# Median normalization

# Median normalization on the rest of the data
sample_median <- apply(log2_mat, 2, median, na.rm = TRUE)
global_median <- median(as.vector(log2_mat), na.rm = TRUE)
median_normalised_data <- sweep(log2_mat, 2, sample_median, FUN = "-") + global_median

# Histogram of log2 and median normalised data

num_cols <- ncol(median_normalised_data)
layout_dims <- ceiling(sqrt(num_cols))  # e.g., 5 if there are 25 columns

par(mfrow = c(layout_dims, ceiling(num_cols / layout_dims)))

for (i in 1:num_cols) {
  hist(median_normalised_data[, i],
       main = paste(colnames(median_normalised_data)[i]),
       xlab = "Log2 Protein Values",
       col = "lightblue", breaks = 30)
}

par(mfrow = c(1,1))  # Reset layout

# Boxplot of normalized data
boxplot(median_normalised_data,
        las = 2,
        col = "lightblue",
        main = "Median normalisation of Log2 intensities",
        ylab = "Log2 Intensity")
abline(h = median(as.vector(median_normalised_data), na.rm = TRUE), col = "red", lty = 2)

# 6. PCA ----
# Filter out proteins with NA and zero variance
noNA <- median_normalised_data[complete.cases(median_normalised_data), ]
var_proteins <- apply(noNA, 1, var, na.rm = TRUE)
var_samples <- apply(noNA, 2, var, na.rm = TRUE)
filtered_data <- noNA[var_proteins > 0, var_samples > 0]

# Extract genotype by removing digits at the end
sample_names <- colnames(filtered_data)
genotype <- gsub("[0-9]+$", "", sample_names)
genotype <- factor(genotype)

pca_res <- prcomp(t(filtered_data), center = TRUE, scale. = TRUE)

genotypes <- as.character(genotype)
col_map <- c("NLFTau" = "#0F80FF", "NLF" = "#0F80FF", "WT" = "black", "Tau" = "black")
bg_map  <- c("NLFTau" = "#0F80FF", "NLF" = "white", "WT" = "white", "Tau" = "black")
pch_map <- c("NLFTau" = 21, "NLF" = 21, "WT" = 21, "Tau" = 21)

point_col <- col_map[genotypes]
point_bg  <- bg_map[genotypes]
point_pch <- pch_map[genotypes]

png("PCA_plot.png", width = 8000, height = 8000, res = 1200)

# Plot without labels or legend
plot(pca_res$x[, 1], pca_res$x[, 2],
     xlab = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
     main = "PCA",
     pch = point_pch, col = point_col, bg = point_bg, cex = 1.5)

dev.off()

# 7. Assessing data exclusion and potential imputation ----
# Function to check if a protein has ≥2 NAs in any group
has_two_or_more_NA_in_any_group <- function(row, group_factor) {
  for (grp in unique(group_factor)) {
    samples_in_group <- which(group_factor == grp)
    na_count <- sum(is.na(row[samples_in_group]))
    if (na_count >= 2) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Identify proteins to exclude
proteins_to_exclude <- rownames(median_normalised_data)[
  apply(median_normalised_data, 1, has_two_or_more_NA_in_any_group, group_factor = genotype)
]

cat("Number of proteins excluded due to >=2 NAs in any group:", length(proteins_to_exclude), "\n")

# Filter the data to keep only proteins passing the filter
median_normalised_data_filtered <- median_normalised_data[!rownames(median_normalised_data) %in% proteins_to_exclude, ]

# Create workbook and add sheets
wb <- createWorkbook()

addWorksheet(wb, "Excluded_Proteins")
writeData(wb, "Excluded_Proteins", data.frame(Protein = proteins_to_exclude))

addWorksheet(wb, "Filtered_Data")
writeData(wb, "Filtered_Data", median_normalised_data_filtered, rowNames = TRUE)

# Save workbook
saveWorkbook(wb, "Proteins_Filtered_And_Excluded.xlsx", overwrite = TRUE)
cat("Workbook saved as 'Proteins_Filtered_And_Excluded.xlsx'\n")

#Saving data file as excel file 
# Convert row names to a new column called "Row"
df1 <- data.frame(Row = rownames(median_normalised_data_filtered),
                 as.data.frame(median_normalised_data_filtered))

# Save with row names included
write_xlsx(df1, "median_normalised_data_filtered.xlsx")


# 8. Limma ----
adjP_cutoff        <- 0.1
logFC_label_cutoff <- 1

sample_names <- colnames(median_normalised_data_filtered)

genotype <- gsub("[0-9]+$", "", sample_names)
genotype <- gsub("/", "_", genotype)
genotype <- factor(genotype)
names(genotype) <- sample_names

design <- model.matrix(~ 0 + genotype)
colnames(design) <- levels(genotype)

fit <- lmFit(median_normalised_data_filtered, design)

contrasts_list <- list(
  NLFTau_vs_Tau = "NLFTau - Tau",
  NLFTau_vs_NLF = "NLFTau - NLF",
  NLF_vs_WT     = "NLF - WT",
  Tau_vs_WT     = "Tau - WT",
  NLFTau_vs_WT  = "NLFTau - WT"
)

dir.create("Volcano_BH", showWarnings = FALSE)

wb <- createWorkbook()
results_list <- list()

for (contrast_name in names(contrasts_list)) {
  
  contrast <- makeContrasts(
    contrasts = contrasts_list[[contrast_name]],
    levels    = design
  )
  
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  results <- topTable(
    fit2,
    coef          = 1,
    number        = Inf,
    adjust.method = "BH",
    sort.by       = "none"
  )
  
  results$Protein <- rownames(results)
  
  results$Status_BH <- "Not Significant"
  results$Status_BH[results$adj.P.Val < adjP_cutoff & results$logFC > 0] <- "Upregulated"
  results$Status_BH[results$adj.P.Val < adjP_cutoff & results$logFC < 0] <- "Downregulated"
  
  results$Status_BH <- factor(
    results$Status_BH,
    levels = c("Downregulated", "Upregulated", "Not Significant")
  )
  
  results_list[[contrast_name]] <- results
  
  addWorksheet(wb, contrast_name)
  writeData(wb, contrast_name, results)
}

plot_volcano_BH <- function(results, contrast_name) {
  
  y_all <- -log10(results$adj.P.Val)
  
  label_data <- if (contrast_name %in% c("NLF_vs_WT", "NLFTau_vs_Tau")) {
    subset(results, adj.P.Val < adjP_cutoff)
  } else {
    subset(results, adj.P.Val < adjP_cutoff & abs(logFC) > logFC_label_cutoff)
  }
  
  label_data$y_lab <- -log10(label_data$adj.P.Val)
  
  p <- ggplot(results, aes(x = logFC, y = y_all)) +
    geom_point(
      aes(color = Status_BH, fill = Status_BH, size = Status_BH),
      shape = 21,
      alpha = 0.6,
      stroke = 0.3
    ) +
    scale_color_manual(values = c(
      "Not Significant" = "gray",
      "Upregulated"     = "black",
      "Downregulated"   = "black"
    )) +
    scale_fill_manual(values = c(
      "Not Significant" = "gray",
      "Upregulated"     = "#FF9999",
      "Downregulated"   = "#A7C7E7"
    )) +
    scale_size_manual(values = c(
      "Not Significant" = 1,
      "Upregulated"     = 3,
      "Downregulated"   = 3
    )) +
    labs(
      x = expression(Log[2] * "(fold change)"),
      y = "-Log10 (adjusted P-value)"
    ) +
    geom_hline(
      yintercept = -log10(adjP_cutoff),
      linetype   = "dashed",
      color      = "red"
    ) +
    scale_y_continuous(trans = "pseudo_log") +
    theme_gray(base_size = 16) +
    theme(
      legend.position  = "none",
      axis.text        = element_text(size = 18),
      axis.title       = element_text(size = 18),
      panel.background = element_rect(fill = "#F6F6F6"),
      axis.line        = element_line(color = "black")
    ) +
    geom_text_repel(
      data          = label_data,
      aes(x = logFC, y = y_lab, label = Protein),
      size          = 5,
      box.padding   = 0.6,
      point.padding = 0.6,
      segment.color = "grey50"
    )
  
  ggsave(
    filename = file.path("Volcano_BH", paste0(contrast_name, "_BH.pdf")),
    plot     = p,
    width    = 8,
    height   = 11
  )
}

for (contrast_name in names(results_list)) {
  plot_volcano_BH(results_list[[contrast_name]], contrast_name)
}

saveWorkbook(
  wb,
  "Limma_Results_BH.xlsx",
  overwrite = TRUE
)

# 9. GO analysis ----
run_go_enrichment <- function(results, contrast_name, adjP_cutoff = 0.1, id_type = "SYMBOL") {
  
  # Split into upregulated and downregulated significant proteins
  up_proteins <- results$Protein[results$adj.P.Val < adjP_cutoff & results$logFC > 0]
  down_proteins <- results$Protein[results$adj.P.Val < adjP_cutoff & results$logFC < 0]
  
  up_proteins <- unique(na.omit(as.character(up_proteins)))
  down_proteins <- unique(na.omit(as.character(down_proteins)))
  
  # Convert human-style uppercase symbols to mouse-style symbols
  convert_to_mouse_style <- function(proteins) {
    paste0(
      toupper(substr(proteins, 1, 1)),
      tolower(substring(proteins, 2))
    )
  }
  
  up_proteins <- convert_to_mouse_style(up_proteins)
  down_proteins <- convert_to_mouse_style(down_proteins)
  
  # Map protein symbols -> ENTREZ IDs
  up_df <- bitr(up_proteins, fromType = id_type, toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  down_df <- bitr(down_proteins, fromType = id_type, toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  up_df <- unique(up_df)
  down_df <- unique(down_df)
  
  ontologies <- c("BP", "MF", "CC")
  wb_go <- createWorkbook()
  
  for (ont in ontologies) {
    
    #### Upregulated proteins ####
    if (nrow(up_df) > 0) {
      ego_up <- enrichGO(
        gene          = up_df$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize     = 3,
        readable      = TRUE
      )
      
      ego_up_df <- as.data.frame(ego_up)
      if (nrow(ego_up_df) > 0) {
        ego_up_df <- ego_up_df[ego_up_df$Count >= 3, ]
      }
    } else {
      ego_up_df <- data.frame()
    }
    
    addWorksheet(wb_go, paste0(ont, "_Up_Proteins"))
    writeData(wb_go, sheet = paste0(ont, "_Up_Proteins"), ego_up_df)
    
    
    #### Downregulated proteins ####
    if (nrow(down_df) > 0) {
      ego_down <- enrichGO(
        gene          = down_df$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize     = 3,
        readable      = TRUE
      )
      
      ego_down_df <- as.data.frame(ego_down)
      if (nrow(ego_down_df) > 0) {
        ego_down_df <- ego_down_df[ego_down_df$Count >= 3, ]
      }
    } else {
      ego_down_df <- data.frame()
    }
    
    addWorksheet(wb_go, paste0(ont, "_Down_Proteins"))
    writeData(wb_go, sheet = paste0(ont, "_Down_Proteins"), ego_down_df)
  }
  
  saveWorkbook(wb_go, file = paste0("GO_Enrichment_", contrast_name, ".xlsx"), overwrite = TRUE)
  message("GO enrichment completed and saved for ", contrast_name)
}

contrast_names_to_run <- c("NLFTau_vs_NLF", "Tau_vs_WT", "NLFTau_vs_WT")

for (contrast_name in contrast_names_to_run) {
  run_go_enrichment(results = results_list[[contrast_name]], contrast_name = contrast_name)
}

# 10. Plotting GO terms for NLFTau vs NLF ----
oob_keep <- scales::oob_keep

# Read data
all_data <- read_excel("GO_Enrichment_NLFTau_vs_NLF_updated.xlsx")

# Rename column so filter works
all_data <- all_data |> rename(ONTOLOGY = Ontology)

# Convert to numeric
all_data$GeneRatioNumeric <- sapply(all_data$GeneRatio, function(x) eval(parse(text = x)))
all_data$Count <- as.numeric(all_data$Count)

# Constants
max_terms   <- 9
fill_limits <- c(0, 0.13)
size_limits <- c(3, 6)
contrast    <- "NLFTau_vs_NLF"

# Choose manually
ont <- "BP"
direction <- "DOWN"

# Filter and process data
plot_data <- all_data %>%
  filter(ONTOLOGY == ont, Direction == direction) %>%
  arrange(p.adjust) %>%
  slice_head(n = max_terms)

# Pad if fewer than max_terms
pad_n <- max_terms - nrow(plot_data)
if (pad_n > 0 && nrow(plot_data) > 0) {
  pad_rows <- plot_data[1, , drop = FALSE]
  pad_rows[,] <- NA
  pad_rows$Description <- ""
  pad_rows <- pad_rows[rep(1, pad_n), ]
  plot_data <- rbind(pad_rows, plot_data)
}

# If completely empty, stop
if (nrow(plot_data) == 0) {
  stop(paste("No data found for", ont, direction))
}

# Slot and label setup
plot_data$slot <- factor(
  paste0("row_", sprintf("%02d", 1:nrow(plot_data))),
  levels = rev(paste0("row_", sprintf("%02d", 1:nrow(plot_data))))
)

y_labels <- plot_data$Description
names(y_labels) <- plot_data$slot

p <- ggplot(
  plot_data,
  aes(
    x    = -log10(p.adjust),
    y    = slot,
    size = Count,
    fill = GeneRatioNumeric
  )
) +
  geom_segment(
    data = plot_data[!is.na(plot_data$p.adjust), ],
    aes(
      x    = -log10(0.1),
      xend = -log10(p.adjust),
      y    = slot,
      yend = slot
    ),
    color     = "#3b3b3b",
    linewidth = 0.7
  ) +
  geom_point(
    data   = plot_data[!is.na(plot_data$p.adjust), ],
    shape  = 21,
    color  = "#3b3b3b",
    stroke = 1
  ) +
  scale_x_continuous(
    limits = c(-log10(0.1), -log10(1e-2)),
    breaks = -log10(c(0.1, 0.01)),
    labels = c("1e-1", "1e-2"),
    expand = c(0.02, 0.02),
    oob    = oob_keep
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_gradientn(
    colors = c("#FFFFCC", "#FFCC66", "#FF9999", "#CC66CC", "#660099"),
    limits = fill_limits,
    breaks = seq(0, 0.15, 0.05),
    name   = "Protein Ratio"
  ) +
  scale_size_continuous(
    limits = size_limits,
    range  = c(3, 6),
    breaks = seq(3, 6, 1),
    name   = "Protein Count"
  ) +
  scale_y_discrete(
    labels = y_labels,
    expand = expansion(mult = c(0.22, 1.0))
  ) +
  theme_minimal() +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y   = element_text(size = 12, color = "black", face = "plain"),
    axis.text.x   = element_text(size = 10, color = "black", face = "plain"),
    axis.title.y  = element_blank(),
    axis.title.x  = element_text(size = 10, color = "black", face = "plain",
                                 margin = margin(t = 10)),
    
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 10),
    
    aspect.ratio  = 3
  )+
  labs(
    x     = "Adjusted P-value (FDR)",
    y     = paste("GO", ont),
    size  = "Protein Count",
    fill  = "Protein Ratio",
    title = paste("Proteomics enrichment:", contrast, "-", ont, "-", direction)
  )
filename <- paste0(contrast, "_", ont, "_", direction, ".pdf")

ggsave(
  filename = filename,
  plot     = p,
  width    = 20,
  height   = 7
)
# 11. Plotting GO terms ----
oob_keep <- scales::oob_keep

# Read data
all_data <- read_excel("GO_Enrichment_Tau_vs_WT_updated.xlsx")

# Rename column so filter works
all_data <- all_data |> rename(ONTOLOGY = Ontology)

# Convert to numeric
all_data$GeneRatioNumeric <- sapply(all_data$GeneRatio, function(x) eval(parse(text = x)))
all_data$Count <- as.numeric(all_data$Count)

# Constants
max_terms   <- 4
fill_limits <- c(0, 0.6)
size_limits <- c(3, 6)
contrast    <- "Tau_vs_WT"

# Choose manually
ont <- "BP"
direction <- "UP" 

# Filter and process data
plot_data <- all_data %>%
  filter(ONTOLOGY == ont, Direction == direction) %>%
  arrange(p.adjust) %>%
  slice_head(n = max_terms)

# Pad if fewer than max_terms
pad_n <- max_terms - nrow(plot_data)
if (pad_n > 0 && nrow(plot_data) > 0) {
  pad_rows <- plot_data[1, , drop = FALSE]
  pad_rows[,] <- NA
  pad_rows$Description <- ""
  pad_rows <- pad_rows[rep(1, pad_n), ]
  plot_data <- rbind(pad_rows, plot_data)
}

# If completely empty, stop
if (nrow(plot_data) == 0) {
  stop(paste("No data found for", ont, direction))
}

# Slot and label setup
plot_data$slot <- factor(
  paste0("row_", sprintf("%02d", 1:nrow(plot_data))),
  levels = rev(paste0("row_", sprintf("%02d", 1:nrow(plot_data))))
)

y_labels <- plot_data$Description
names(y_labels) <- plot_data$slot

p <- ggplot(
  plot_data,
  aes(
    x    = -log10(p.adjust),
    y    = slot,
    size = Count,
    fill = GeneRatioNumeric
  )
) +
  geom_segment(
    data = plot_data[!is.na(plot_data$p.adjust), ],
    aes(
      x    = -log10(0.1),
      xend = -log10(p.adjust),
      y    = slot,
      yend = slot
    ),
    color     = "#3b3b3b",
    linewidth = 0.7
  ) +
  geom_point(
    data   = plot_data[!is.na(plot_data$p.adjust), ],
    shape  = 21,
    color  = "#3b3b3b",
    stroke = 1
  ) +
  scale_x_continuous(
    limits = c(-log10(0.1), -log10(1e-4)),
    breaks = -log10(c(0.1, 0.01, 0.001, 0.0001)),
    labels = c("1e-1", "1e-2", "1e-3", "1e-4"),
    expand = c(0.02, 0.02),
    oob    = oob_keep
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_gradientn(
    colors = c("#FFFFCC", "#FFCC66", "#FF9999", "#CC66CC", "#660099"),
    limits = fill_limits,
    breaks = seq(0, 0.6, 0.2),
    name   = "Protein Ratio"
  ) +
  scale_size_continuous(
    limits = size_limits,
    range  = c(3, 6),
    breaks = seq(3, 6, 1),
    name   = "Protein Count"
  ) +
  scale_y_discrete(
    labels = y_labels,
    expand = expansion(mult = c(0.22, 4.5))
  ) +
  theme_minimal() +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y   = element_text(size = 12, color = "black", face = "plain"),
    axis.text.x   = element_text(size = 10, color = "black", face = "plain"),
    axis.title.y  = element_blank(),
    axis.title.x  = element_text(size = 10, color = "black", face = "plain",
                                 margin = margin(t = 10)),
    
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 10),
    
    aspect.ratio  = 3
  )+
  labs(
    x     = "Adjusted P-value (FDR)",
    y     = paste("GO", ont),
    size  = "Protein Count",
    fill  = "Protein Ratio",
    title = paste("Proteomics enrichment:", contrast, "-", ont, "-", direction)
  )
filename <- paste0(contrast, "_", ont, "_", direction, ".pdf")

ggsave(
  filename = filename,
  plot     = p,
  width    = 20,
  height   = 7
)
# 12. UpSet plot for shared DEPs ----
dep_sets <- purrr::imap(
  results_list,
  ~ .x %>%
    dplyr::filter(!is.na(Protein), Protein != "") %>%
    dplyr::filter(adj.P.Val < adjP_cutoff) %>%
    dplyr::pull(Protein) %>%
    unique()
)

# Remove unwanted contrast
dep_sets <- dep_sets[names(dep_sets) != "NLFTau_vs_WT"]

# Remove empty sets
dep_sets <- dep_sets[lengths(dep_sets) > 0]

# Convert list to UpSetR input
upset_df_prot <- UpSetR::fromList(dep_sets)

# Order sets by number of DEPs, largest at bottom like your example
set_order <- names(sort(colSums(upset_df_prot), decreasing = TRUE))

pdf("UpSet_shared_DEPs_BH_UpSetR.pdf", width = 5, height = 5)

UpSetR::upset(
  upset_df_prot,
  sets = set_order,
  keep.order = TRUE,
  
  order.by = "freq",
  
  mainbar.y.label = "Intersection size",
  sets.x.label = "DEPs per contrast",
  
  main.bar.color = "gray20",
  sets.bar.color = "gray20",
  matrix.color = "gray20",
  shade.color = "gray95",
  
  point.size = 3,
  line.size = 1,
  
  mb.ratio = c(0.65, 0.35),
  text.scale = c(1.8, 1.8, 1.5, 1.5, 1.8, 1.5),
  number.angles = 0
)

dev.off() 
# 13. WGCNA ----
cl <- makeCluster(8)
registerDoParallel(cl)
enableWGCNAThreads(nThreads = 8)

options(stringsAsFactors = FALSE)

# Transpose: rows = samples, columns = proteins
datExpr <- t(median_normalised_data_filtered)

# Check dimensions
dim(datExpr)  

#Remove any NA proteins because can't be handled 
datExpr <- datExpr[, colSums(is.na(datExpr)) == 0]
# Check that no NAs remain
sum(is.na(datExpr))  # Should be 0

#Deciding pick soft threshold value 
powers <- c(1:20)

# Run scale-free topology analysis
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale-Free Topology Model Fit (signed R²)", 
     type = "n", main = "Scale Independence")

text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red")

abline(h = 0.85, col = "blue", lty = 2)  # Threshold for good scale-free fit, R^2 needs to be above 0.85 minimum number above this is the chosen soft-thresholding power number 

# --- Build modules (your prefs) ---
net <- blockwiseModules(
  datExpr,
  power              = 6,
  corType            = "bicor",
  networkType        = "signed",
  TOMType            = "signed",
  TOMDenom           = "mean",
  deepSplit          = 2,
  minModuleSize      = 30,
  pamRespectsDendro  = TRUE,
  mergeCutHeight     = 0.25,
  reassignThreshold  = 0.05,
  minKMEtoStay       = 0.30,
  numericLabels      = TRUE,
  saveTOMs           = FALSE,
  verbose            = 5
)

moduleColors_init <- labels2colors(net$colors)
names(moduleColors_init) <- colnames(datExpr)     # <-- give it names immediately

colors_clean <- moduleColors_init                 # cleaning works on a *named* vector

# When building ME list for cleaning, ignore grey to be safe:
MEs_init <- orderMEs(moduleEigengenes(datExpr, colors = moduleColors_init)$eigengenes)
kME_init <- signedKME(datExpr, MEs_init, outputColumnName = "kME")
mods <- sub("^ME", "", colnames(MEs_init))
mods <- setdiff(mods, "grey")                     # don't try to "clean" grey

for (mod in mods) {
  kcol <- paste0("kME", mod)
  in_mod <- which(colors_clean == mod)
  if (!kcol %in% colnames(kME_init) || length(in_mod) == 0) next
  drop_idx <- in_mod[ is.na(kME_init[in_mod, kcol]) | (kME_init[in_mod, kcol] < 0.30) ]
  if (length(drop_idx)) colors_clean[drop_idx] <- "grey"
}

# --- NEW: Handle tiny modules (<30) by reassignment ---
tab_clean <- table(colors_clean)
tiny_mods <- names(tab_clean)[tab_clean < 30 & names(tab_clean) != "grey"]

if (length(tiny_mods)) {
  for (tm in tiny_mods) {
    tm_idx <- which(colors_clean == tm)
    kME_tm <- kME_init[tm_idx, , drop = FALSE]
    best_ix  <- apply(kME_tm, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
    best_val <- apply(kME_tm, 1, max, na.rm = TRUE)
    
    for (i in seq_along(tm_idx)) {
      if (!is.na(best_ix[i]) && is.finite(best_val[i]) && best_val[i] >= 0.30) {
        best_mod <- sub("^kME", "", colnames(kME_tm)[best_ix[i]])
        colors_clean[tm_idx[i]] <- best_mod
      } else {
        colors_clean[tm_idx[i]] <- "grey"
      }
    }
  }
}

# --- Final eigengenes & kME for downstream ---
MEs_final <- orderMEs(moduleEigengenes(datExpr, colors = colors_clean)$eigengenes)
kME_final <- signedKME(datExpr, MEs_final, outputColumnName = "kME")

# --- Quick sanity checks ---
cat("Modules (pre-clean):", length(unique(moduleColors_init[moduleColors_init!="grey"])), "\n")
cat("Modules (final):    ", length(unique(colors_clean[colors_clean!="grey"])), "\n")

# Save dendrogram to PNG
png("Final_Module_Dendrogram.png", width = 2000, height = 1500, res = 300)
plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = labels2colors(net$colors[net$blockGenes[[1]]]),
  groupLabels = "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

# Saving ME data
me_df <- MEs_final %>%
  as.data.frame() %>%
  rownames_to_column("Sample")   # Sample IDs like WT1, WT2, ...

# Write directly to Excel
wb <- createWorkbook()
addWorksheet(wb, "ME_values")
writeData(wb, "ME_values", me_df)
saveWorkbook(wb, "ME_values_by_sample.xlsx", overwrite = TRUE)

message("Wrote: ME_values_by_sample.xlsx")

# 14. WGCNA GO analysis ----
# Define modules of interest
modules_of_interest <- c("red", "brown", "blue", "pink")

# Get gene-module assignments
gene_modules <- data.frame(
  gene = colnames(datExpr),
  module = labels2colors(net$colors)
)

# Extract genes per module
selected_genes <- lapply(modules_of_interest, function(mod) {
  gene_modules$gene[gene_modules$module == mod]
})
names(selected_genes) <- modules_of_interest

# Run GO enrichment and store in all_ego

all_ego <- list()
wb <- createWorkbook()

for (mod in names(selected_genes)) {
  message("Running GO for module: ", mod)
  genes <- stringr::str_to_title(selected_genes[[mod]])
  valid_genes <- genes[genes %in% keys(org.Mm.eg.db, keytype = "SYMBOL")]
  if (length(valid_genes) < 5) next
  
  gene_df <- bitr(valid_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  entrez_ids <- gene_df$ENTREZID
  
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrichGO(
      gene          = entrez_ids,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      minGSSize     = 5,
      readable      = TRUE
    )
    
    all_ego[[mod]][[ont]] <- ego
    
    ego_df <- as.data.frame(ego)
    ego_df_filtered <- ego_df %>% filter(p.adjust < 0.05 & Count >= 5)
    
    addWorksheet(wb, paste0(mod, "_", ont))
    writeData(wb, paste0(mod, "_", ont), ego_df_filtered)
  }
}

# Save Excel workbook
saveWorkbook(wb, "GO_Enrichment_SelectedModules.xlsx", overwrite = TRUE)

# Getting list of proteins within modules 
modules_of_interest <- c("red", "brown", "blue", "pink")

# Get gene-module assignments
gene_modules <- data.frame(
  gene = colnames(datExpr),
  module = labels2colors(net$colors)
)

# Extract genes per module
selected_genes <- lapply(modules_of_interest, function(mod) {
  gene_modules$gene[gene_modules$module == mod]
})
names(selected_genes) <- modules_of_interest

# === Save to Excel ===
wb <- createWorkbook()

for (mod in names(selected_genes)) {
  addWorksheet(wb, mod)
  writeData(wb, mod, data.frame(Protein = selected_genes[[mod]]))
}

saveWorkbook(wb, file = "Modules_protein_lists.xlsx", overwrite = TRUE)

# 15. PPI Network ----
# Parameters
modules_of_interest <- c("brown", "red", "blue", "pink")  # Your selected modules
threshold <- 0.7  # Correlation cutoff for edges
module_colors <- labels2colors(net$colors)  # Convert numeric module labels to color names

gene_modules <- data.frame(
  Protein = colnames(datExpr),
  ModuleColor = module_colors,
  stringsAsFactors = FALSE
)

# Ensure column names in datExpr are all uppercase to match later comparisons
colnames(datExpr) <- toupper(colnames(datExpr))

# Recalculate eigengenes using color labels
MEs_colored <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
# Loop through each module
for (mod in modules_of_interest) {
  
  # 1. Get gene names in this module
  mod_genes <- gene_modules %>%
    filter(ModuleColor == mod) %>%
    pull(Protein) %>%
    toupper()  # Ensure case matches datExpr colnames
  
  mod_expr <- datExpr[, colnames(datExpr) %in% mod_genes, drop = FALSE]
  
  # 2. Compute correlation matrix
  cor_mat <- cor(mod_expr, method = "pearson", use = "pairwise.complete.obs")
  diag(cor_mat) <- 0  # Remove self-correlations
  
  # 3. Apply correlation threshold
  cor_mat[abs(cor_mat) < threshold] <- 0
  
  # 4. Convert to edge list
  edge_df <- as.data.frame(as.table(cor_mat))
  colnames(edge_df) <- c("Source", "Target", "Weight")
  edge_df <- edge_df %>%
    filter(Weight != 0) %>%
    filter(as.character(Source) < as.character(Target))  # Avoid duplicate undirected edges
  
  # 5.5. Compute degree from edge list
  degree_df <- edge_df %>%
    select(Source, Target) %>%
    pivot_longer(cols = everything(), values_to = "Protein") %>%
    count(Protein, name = "Degree")
  
  # 6. Calculate kME: correlation with module eigengene
  module_eigengene <- MEs_colored[, paste0("ME", mod)]
  kME_values <- apply(mod_expr, 2, function(gene_expr) cor(gene_expr, module_eigengene, use = "pairwise.complete.obs"))
  
  # 7. Create node table with kME and degree
  node_df <- data.frame(
    name = names(kME_values),
    Label = names(kME_values),
    kME = kME_values,
    Module = mod,
    stringsAsFactors = FALSE
  ) %>%
    left_join(degree_df, by = c("name" = "Protein")) %>%
    mutate(Degree = ifelse(is.na(Degree), 0, Degree))  # 0 if not connected to others
  
  # Also update edge list to only include remaining nodes
  edge_df <- edge_df %>%
    filter(Source %in% node_df$name & Target %in% node_df$name)
  
  # 8. Save edge and node tables for Cytoscape
  write.table(node_df, paste0("Nodes_", mod, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(edge_df, paste0("Edges_", mod, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# 16. Plotting WGCNA GO terms ----
oob_keep <- scales::oob_keep

# Read data
all_data <- readxl::read_excel("WGCNA_GO_terms.xlsx")

# Convert to numeric
all_data$GeneRatioNumeric <- sapply(all_data$GeneRatio, function(x) eval(parse(text = x)))
all_data$Count <- as.numeric(all_data$Count)

# Constants
max_terms   <- 20
fill_limits <- c(0, 0.10)
size_limits <- c(0, 20)

# Choose manually
module <- "Red"
# Filter and process data
plot_data <- all_data %>%
  dplyr::filter(Module == module) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = max_terms)

# Pad if fewer than max_terms
pad_n <- max_terms - nrow(plot_data)

if (pad_n > 0 && nrow(plot_data) > 0) {
  pad_rows <- plot_data[1, , drop = FALSE]
  pad_rows[,] <- NA
  pad_rows$Description <- ""
  pad_rows <- pad_rows[rep(1, pad_n), ]
  plot_data <- rbind(pad_rows, plot_data)
}

# If completely empty, stop
if (nrow(plot_data) == 0) {
  stop(paste("No GO terms found for module:", module))
}

# Slot and label setup
plot_data$slot <- factor(
  paste0("row_", sprintf("%02d", 1:nrow(plot_data))),
  levels = rev(paste0("row_", sprintf("%02d", 1:nrow(plot_data))))
)

y_labels <- plot_data$Description
names(y_labels) <- plot_data$slot

# Plot
p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x    = -log10(p.adjust),
    y    = slot,
    size = Count,
    fill = GeneRatioNumeric
  )
) +
  ggplot2::geom_segment(
    data = plot_data[!is.na(plot_data$p.adjust), ],
    ggplot2::aes(
      x    = -log10(0.1),
      xend = -log10(p.adjust),
      y    = slot,
      yend = slot
    ),
    color     = "#3b3b3b",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    data   = plot_data[!is.na(plot_data$p.adjust), ],
    shape  = 21,
    color  = "#3b3b3b",
    stroke = 1
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-log10(0.1), -log10(1e-7)),
    breaks = -log10(c(0.1, 1e-3, 1e-5, 1e-7)),
    labels = c("1e-1", "1e-3", "1e-5", "1e-7"),
    expand = c(0.02, 0.02),
    oob    = oob_keep
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::scale_fill_gradientn(
    colors = c("#FFFFCC", "#FFCC66", "#FF9999", "#CC66CC", "#660099"),
    limits = fill_limits,
    breaks = seq(0, 0.10, 0.05),
    name   = "Protein Ratio"
  ) +
  ggplot2::scale_size_continuous(
    limits = size_limits,
    range  = c(1.5, 8),
    breaks = seq(5, 20, 5),
    name   = "Protein Count"
  ) +
  ggplot2::scale_y_discrete(
    labels = y_labels,
    expand = ggplot2::expansion(mult = c(0.22, 1.0))
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.border  = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y   = ggplot2::element_text(size = 12, color = "black", face = "plain"),
    axis.text.x   = ggplot2::element_text(size = 10, color = "black", face = "plain"),
    axis.title.y  = ggplot2::element_blank(),
    axis.title.x  = ggplot2::element_text(
      size = 10,
      color = "black",
      face = "plain",
      margin = ggplot2::margin(t = 10)
    ),
    
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.title     = ggplot2::element_text(size = 12),
    legend.text      = ggplot2::element_text(size = 10),
    
    aspect.ratio = 3
  ) +
  ggplot2::labs(
    x     = "Adjusted P-value (FDR)",
    y     = "GO BP",
    size  = "Protein Count",
    fill  = "Protein Ratio",
    title = paste("WGCNA GO BP enrichment:", module)
  )

# Save
filename <- paste0("WGCNA_GO_BP_", module, ".pdf")

ggplot2::ggsave(
  filename = filename,
  plot     = p,
  width    = 20,
  height   = 7
)