# 1. Install Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Bioconductor packages
BiocManager::install(c(
  "limma", "edgeR", "IHW", "org.Mm.eg.db", "biomaRt",
  "clusterProfiler", "ReactomePA", "GOSemSim",
  "impute", "WGCNA", "preprocessCore"
), update = TRUE, ask = FALSE)

## CRAN packages
install.packages(c(
  "tidyverse", "openxlsx", "readxl", "writexl",
  "ggrepel", "scales", "doParallel", "broom", "ggbreak"
))
install.packages(c("ComplexUpset", "dplyr", "ggplot2"))
options(stringsAsFactors = FALSE)

# 2. Load libraries ----

library(tidyverse)   
library(openxlsx)
library(readxl)
library(writexl)

library(edgeR)
library(limma)
library(biomaRt)

library(ggrepel)
library(scales)
library(broom)
library(dplyr)

# 3. Load data ----
load("Gene_Quantification_Output_2021.RDA")
gene_expression_data <- txi_lengthScaled
sample_names <- colnames(gene_expression_data$counts) #extract Sample_ID names
sample_names_df <- data.frame(Sample_ID = sample_names)
write.xlsx(sample_names_df, file = "Sample_IDs_from_counts.xlsx", rowNames = FALSE)
sample_list <- read_excel("24_month_to_keep.xlsx") 
samples_to_keep <- sample_list$Sample_ID #Read in 24-month-old information and only include those IDs
all(samples_to_keep %in% colnames(gene_expression_data$counts))
gene_expression_subset <- list()

gene_expression_subset$counts <-
  gene_expression_data$counts[, samples_to_keep, drop = FALSE]

gene_expression_subset$abundance <-
  gene_expression_data$abundance[, samples_to_keep, drop = FALSE]

gene_expression_subset$length <-
  gene_expression_data$length[, samples_to_keep, drop = FALSE]
gene_expression_subset$infReps <-
  gene_expression_data$infReps[samples_to_keep]
gene_expression_subset$countsFromAbundance <-
  gene_expression_data$countsFromAbundance

# 4. Filter genes ----
counts_24m <- gene_expression_subset$counts
dge_24m <- DGEList(counts = counts_24m)
group_24m <- factor(sample_list$Genotype)
design_filter <- model.matrix(~ 0 + group_24m)
colnames(design_filter) <- levels(group_24m)
keep_genes <- filterByExpr(dge_24m, design = design_filter)
table(keep_genes)
dge_24m <- dge_24m[keep_genes, , keep.lib.sizes = FALSE]

ensembl_ids   <- rownames(dge_24m)
any(duplicated(ensembl_ids))

ensembl_ids <- rownames(dge_24m)

options(timeout = 300)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ann <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = mart
)

ann_pc <- ann[ann$gene_biotype == "protein_coding", ]

ensembl_ids <- rownames(dge_24m)
pc_ids <- intersect(ensembl_ids, ann_pc$ensembl_gene_id)

dge_24m_pc <- dge_24m[pc_ids, , keep.lib.sizes = FALSE]

ann_pc_ordered <- ann_pc[match(rownames(dge_24m_pc), ann_pc$ensembl_gene_id), ]

dge_24m_pc$genes <- data.frame(
  ensembl_id  = ann_pc_ordered$ensembl_gene_id,
  mgi_symbol  = ann_pc_ordered$mgi_symbol,
  biotype     = ann_pc_ordered$gene_biotype,
  stringsAsFactors = FALSE
)

nrow(dge_24m_pc) # this is to check how many protein-coding genes they are

# 5. Normalisation ----
dge_24m_pc <- calcNormFactors(dge_24m_pc, method = "TMM")
dge_24m_pc$samples[, c("lib.size", "norm.factors")]

hk_symbols <- c("Actb", "Actg1")
hk_rows <- which(dge_24m_pc$genes$mgi_symbol %in% hk_symbols)

dge_24m_pc$genes[hk_rows, ]

hk_cpm <- cpm(dge_24m_pc, log = FALSE, normalized.lib.sizes = TRUE)[hk_rows, , drop = FALSE] # getting TMM CPM for HK genes

hk_ref <- colMeans(hk_cpm)
hk_scale <- hk_ref / exp(mean(log(hk_ref))) # HK normalisation 

dge_24m_hk <- dge_24m_pc
dge_24m_hk$samples$norm.factors <- dge_24m_pc$samples$norm.factors * hk_scale

dge_24m_hk$samples$norm.factors <- dge_24m_hk$samples$norm.factors /
  exp(mean(log(dge_24m_hk$samples$norm.factors)))

dge_24m_hk$samples[, c("lib.size", "norm.factors")]

hk_cpm_after <- cpm(dge_24m_hk, log = FALSE, normalized.lib.sizes = TRUE)[hk_rows, , drop = FALSE]

# 6. Mitochondrial limma analysis ----

# Read mitochondrial gene list from Excel
mito_df <- read_excel("Mitochondria.xlsx")

mito_symbols <- mito_df$Symbol
mito_symbols <- unique(na.omit(trimws(mito_symbols)))

length(mito_symbols)
head(mito_symbols)

# Match mitochondrial genes to protein-coding gene list 
gene_symbols_all <- dge_24m_hk$genes$mgi_symbol

matched_mito <- gene_symbols_all %in% mito_symbols

table(matched_mito)

# subset DGEList to matched mitochondrial genes
dge_24m_mito <- dge_24m_hk[matched_mito, , keep.lib.sizes = TRUE]

#  Limma-voom differential expression analysis

adjP_cutoff  <- 0.05
logFC_cutoff <- 0

genotype_24m <- factor(
  sample_list$Genotype,
  levels = c("WT", "NLF", "hTau", "NLFhTau")
)

design_24m_hk <- model.matrix(~ 0 + genotype_24m)
colnames(design_24m_hk) <- levels(genotype_24m)

contrast_matrix_24m_hk <- limma::makeContrasts(
  hTau_vs_WT      = hTau - WT,
  NLF_vs_WT       = NLF - WT,
  NLFhTau_vs_NLF  = NLFhTau - NLF,
  NLFhTau_vs_hTau = NLFhTau - hTau,
  NLFhTau_vs_WT   = NLFhTau - WT,
  levels = design_24m_hk
)

v_24m_mito_hk <- limma::voom(
  dge_24m_mito,
  design_24m_hk,
  plot = TRUE
)

fit_24m_mito_hk <- limma::lmFit(v_24m_mito_hk, design_24m_hk)

fit2_24m_mito_hk <- fit_24m_mito_hk %>%
  limma::contrasts.fit(contrast_matrix_24m_hk) %>%
  limma::eBayes()

# Extract and export all contrast results
get_limma_results <- function(fit, contrast_name, adjP_cutoff = 0.05, logFC_cutoff = 0) {
  
  results <- limma::topTable(
    fit,
    coef          = contrast_name,
    number        = Inf,
    adjust.method = "BH",
    sort.by       = "none"
  )
  
  results <- results %>%
    tibble::rownames_to_column("Ensembl_ID") %>%
    dplyr::select(-dplyr::any_of(c("ensembl_id", "biotype"))) %>%
    dplyr::mutate(
      Status_BH = dplyr::case_when(
        adj.P.Val < adjP_cutoff & logFC >  logFC_cutoff ~ "Upregulated",
        adj.P.Val < adjP_cutoff & logFC < -logFC_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Status_BH = factor(
        Status_BH,
        levels = c("Downregulated", "Upregulated", "Not Significant")
      )
    ) %>%
    dplyr::relocate(dplyr::any_of(c("Ensembl_ID", "mgi_symbol")))
  
  return(results)
}

contrast_names <- colnames(contrast_matrix_24m_hk)

results_list_mito <- purrr::map(
  contrast_names,
  ~ get_limma_results(
    fit           = fit2_24m_mito_hk,
    contrast_name = .x,
    adjP_cutoff   = adjP_cutoff,
    logFC_cutoff  = logFC_cutoff
  )
)

names(results_list_mito) <- contrast_names

# Export to Excel

wb_mito <- openxlsx::createWorkbook()

purrr::iwalk(results_list_mito, function(results, contrast_name) {
  openxlsx::addWorksheet(wb_mito, contrast_name)
  openxlsx::writeData(wb_mito, contrast_name, results)
})

openxlsx::saveWorkbook(
  wb_mito,
  "RNAseq_DEG_24m_mito_all_contrasts_HKscaled.xlsx",
  overwrite = TRUE
)

# Volcano plot for selected contrast

volcano_contrast <- "NLFhTau_vs_NLF"

genes_to_label <- c(
  "Cox8b",
  "mt-Atp8",
  "mt-Nd3",
  "Slc25a48"
)

label_logFC_cutoff <- -0.5

results_volcano <- results_list_mito[[volcano_contrast]] %>%
  dplyr::mutate(
    adj.P.Val.plot = pmax(adj.P.Val, .Machine$double.xmin),
    negLog10AdjP   = -log10(adj.P.Val.plot),
    Significance   = Status_BH,
    Label = dplyr::case_when(
      adj.P.Val < adjP_cutoff & logFC < label_logFC_cutoff ~ mgi_symbol,
      mgi_symbol %in% genes_to_label ~ mgi_symbol,
      TRUE ~ NA_character_
    )
  )
p_volcano <- ggplot2::ggplot(
  results_volcano,
  ggplot2::aes(x = logFC, y = negLog10AdjP)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = Significance),
    alpha = 0.8,
    size = 2
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Downregulated" = "blue",
      "Upregulated" = "red",
      "Not Significant" = "grey75"
    )
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(adjP_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50",
    alpha = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(results_volcano, !is.na(Label)),
    ggplot2::aes(label = Label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    seed = 123
  ) +
  ggplot2::coord_cartesian(xlim = c(-2.5, 0.5)) +
  ggplot2::labs(
    title = volcano_contrast,
    x = expression(Log[2]*"(fold change)"),
    y = expression(-Log[10]*"(adjusted P-value)")
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
    axis.title = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank()
  )

print(p_volcano)

ggplot2::ggsave(
  filename = "Volcano_NLFhTau_vs_NLF_mito.pdf",
  plot = p_volcano,
  width = 6,
  height = 5,
  bg = "white"
)

# 7. Export significant limma genes + expression values for Prism  MITO ----
export_prism_data <- function(
    fit,
    voom_object,
    sample_metadata,
    contrast_name,
    output_file,
    adjP_cutoff = 0.05,
    genotype_levels = c("WT", "NLF", "hTau", "NLFhTau")
) {
  
  # Extract limma results
  results <- limma::topTable(
    fit,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  ) %>%
    tibble::rownames_to_column("Ensembl_ID")
  
  # Keep significant genes only
  sig_genes <- results %>%
    dplyr::filter(adj.P.Val < adjP_cutoff)
  
  # Extract voom expression matrix
  expr_mat <- voom_object$E
  
  # Keep genes present in expression matrix
  sig_gene_ids_present <- intersect(sig_genes$Ensembl_ID, rownames(expr_mat))
  
  expr_sig <- expr_mat[sig_gene_ids_present, , drop = FALSE]
  
  # Match limma results to expression matrix row order
  sig_genes <- sig_genes %>%
    dplyr::filter(Ensembl_ID %in% sig_gene_ids_present) %>%
    dplyr::slice(match(rownames(expr_sig), Ensembl_ID))
  
  stopifnot(identical(sig_genes$Ensembl_ID, rownames(expr_sig)))
  
  # Build sample metadata
  sample_info <- data.frame(
    Sample_ID = colnames(expr_sig),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      sample_metadata %>%
        dplyr::select(Sample_ID, Genotype),
      by = "Sample_ID"
    ) %>%
    dplyr::mutate(
      Genotype = factor(Genotype, levels = genotype_levels)
    )
  
  if (any(is.na(sample_info$Genotype))) {
    stop("Some expression columns could not be matched to sample metadata.")
  }
  
  # Reorder samples by genotype
  sample_info <- sample_info %>%
    dplyr::arrange(Genotype, Sample_ID) %>%
    dplyr::mutate(
      Prism_Colname = paste(Genotype, Sample_ID, sep = "_")
    )
  
  expr_sig <- expr_sig[, sample_info$Sample_ID, drop = FALSE]
  colnames(expr_sig) <- sample_info$Prism_Colname
  
  # Build export table
  export_table <- data.frame(
    Ensembl_ID = sig_genes$Ensembl_ID,
    mgi_symbol = sig_genes$mgi_symbol,
    logFC = sig_genes$logFC,
    adj.P.Val = sig_genes$adj.P.Val,
    expr_sig,
    check.names = FALSE
  ) %>%
    dplyr::arrange(adj.P.Val)
  
  # Export
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "SigGenes_PlotData")
  openxlsx::writeData(wb, "SigGenes_PlotData", export_table)
  
  openxlsx::saveWorkbook(
    wb,
    file = output_file,
    overwrite = TRUE
  )
  
  return(export_table)
}

# 8. Functional analysis ----

# MitoCarta annotation + pathway summary for downregulated genes
contrast_name <- "NLFhTau_vs_NLF"

# 1. Read MitoCarta annotation
mito_annot <- readxl::read_excel(
  "Mouse.MitoCarta3.0.xls",
  sheet = "A Mouse MitoCarta3.0"
) %>%
  dplyr::select(
    Symbol,
    MitoCarta3.0_MitoPathways,
    MitoCarta3.0_SubMitoLocalization
  ) %>%
  dplyr::filter(!is.na(Symbol), Symbol != "") %>%
  dplyr::distinct(Symbol, .keep_all = TRUE)

# 2. Use existing limma result table
sig_down_results <- results_list_mito[[contrast_name]] %>%
  dplyr::filter(
    adj.P.Val < adjP_cutoff,
    logFC < -logFC_cutoff
  )

# 3. Annotate and clean pathway/localisation labels
final_table <- sig_down_results %>%
  dplyr::left_join(mito_annot, by = c("mgi_symbol" = "Symbol")) %>%
  dplyr::mutate(
    TopPathway = stringr::str_trim(
      stringr::str_extract(MitoCarta3.0_MitoPathways, "^[^>]+")
    ),
    TopPathway = dplyr::case_when(
      TopPathway == "Mitochondrial central dogma" ~
        "mt-DNA/RNA maintenance, metabolism & translation",
      TopPathway == "OXPHOS" ~
        "Oxidative phosphorylation",
      TRUE ~ TopPathway
    ),
    SubMito_clean = dplyr::case_when(
      MitoCarta3.0_SubMitoLocalization == "MIM" ~ "Inner membrane",
      MitoCarta3.0_SubMitoLocalization == "MOM" ~ "Outer membrane",
      MitoCarta3.0_SubMitoLocalization == "IMS" ~ "Intermembrane space",
      MitoCarta3.0_SubMitoLocalization == "Matrix" ~ "Matrix",
      MitoCarta3.0_SubMitoLocalization == "Membrane" ~ "Membrane",
      is.na(MitoCarta3.0_SubMitoLocalization) |
        MitoCarta3.0_SubMitoLocalization == "" ~ "Unknown",
      TRUE ~ MitoCarta3.0_SubMitoLocalization
    )
  ) %>%
  dplyr::select(
    Ensembl_ID,
    mgi_symbol,
    logFC,
    AveExpr,
    P.Value,
    adj.P.Val,
    TopPathway,
    MitoCarta3.0_SubMitoLocalization,
    SubMito_clean
  ) %>%
  dplyr::arrange(TopPathway, adj.P.Val)

# 4. Export gene-level table
openxlsx::write.xlsx(
  final_table,
  file = "NLFhTau_vs_NLF_downregulated_mitochondrial_genes_topPathway_with_SubMito.xlsx",
  rowNames = FALSE
)
# Summarise genes by pathway and sub-mito localisation
stacked_table <- final_table %>%
  dplyr::filter(
    !is.na(TopPathway),
    TopPathway != "",
    TopPathway != "0",
    !is.na(SubMito_clean),
    SubMito_clean != ""
  ) %>%
  dplyr::count(TopPathway, SubMito_clean, name = "n_genes")

pathway_order <- stacked_table %>%
  dplyr::group_by(TopPathway) %>%
  dplyr::summarise(total_genes = sum(n_genes), .groups = "drop") %>%
  dplyr::arrange(total_genes)

stacked_table <- stacked_table %>%
  dplyr::mutate(
    TopPathway = factor(TopPathway, levels = pathway_order$TopPathway)
  )

openxlsx::write.xlsx(
  stacked_table,
  file = "NLFhTau_vs_NLF_downregulated_mito_pathway_by_SubMito.xlsx",
  rowNames = FALSE
)

# Plot
submito_colors <- c(
  "Matrix"              = "#1D4E89",
  "Inner membrane"      = "#2A9D8F",
  "Outer membrane"      = "#4EA8DE",
  "Intermembrane space" = "#7B9ACC",
  "Membrane"            = "#A8DADC",
  "Unknown"             = "#D9D9D9"
)

p <- ggplot2::ggplot(
  stacked_table,
  ggplot2::aes(
    x = TopPathway,
    y = n_genes,
    fill = SubMito_clean
  )
) +
  ggplot2::geom_col(
    width = 0.82,
    colour = "white",
    linewidth = 0.3
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.06))
  ) +
  ggplot2::scale_fill_manual(values = submito_colors) +
  ggplot2::labs(
    x = NULL,
    y = "Number of genes",
    fill = "Sub-mito localisation"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 11),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 11),
    axis.title.x = ggplot2::element_text(size = 10),
    axis.title.y = ggplot2::element_text(size = 10),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    legend.position = "right",
    plot.margin = ggplot2::margin(5.5, 15, 5.5, 5.5)
  )

print(p)

ggplot2::ggsave(
  filename = "Downregulated_mitochondrial_genes_by_major_pathway_stacked_SubMito.pdf",
  plot = p,
  width = 8,
  height = 3,
  bg = "white"
)

# 9. Synaptic limma analysis ---- 
human_df <- read_excel("synaptic_genes.xlsx")

human_ensembl <- human_df$ensembl_id

# Remove NA, spaces, duplicates
human_ensembl <- unique(na.omit(trimws(human_ensembl)))

# Remove Ensembl version suffix if present (e.g. ENSG000001234.5 -> ENSG000001234)
human_ensembl <- sub("\\..*$", "", human_ensembl)

orth <- read_excel("mouse_human_orthologues_ensembl115.xls") %>%
  filter(`Human homology type` == "ortholog_one2one") %>%
  transmute(
    ensembl_id_mouse = sub("\\..*$", "", `Gene stable ID`),
    ensembl_id_human = sub("\\..*$", "", `Human gene stable ID`)
  ) %>%
  distinct()

# MAP HUMAN LIST TO MOUSE ORTHOLOGS 
human_syn_df <- data.frame(
  ensembl_id_human = human_ensembl,
  stringsAsFactors = FALSE
)

syn_orth <- human_syn_df %>%
  inner_join(orth, by = "ensembl_id_human")

# MATCH ORTHOLOGS TO HK-SCALED RNA-seq OBJECT

dge_24m_hk$genes$ensembl_id <- sub("\\..*$", "", dge_24m_hk$genes$ensembl_id)

mouse_syn_ensembl <- unique(syn_orth$ensembl_id_mouse)

matched_syn <- dge_24m_hk$genes$ensembl_id %in% mouse_syn_ensembl
table(matched_syn)

dge_24m_syn <- dge_24m_hk[matched_syn, , keep.lib.sizes = TRUE]
dge_24m_syn$genes <- dge_24m_hk$genes[matched_syn, , drop = FALSE]

# SAVE MAPPING TABLE

mapping_in_dataset <- syn_orth %>%
  dplyr::inner_join(
    dge_24m_syn$genes %>%
      dplyr::select(ensembl_id, mgi_symbol),
    by = c("ensembl_id_mouse" = "ensembl_id")
  )

openxlsx::write.xlsx(
  mapping_in_dataset,
  file = "SynGO_human_to_mouse_mapping_in_dataset.xlsx",
  rowNames = FALSE
)

# DESIGN MATRIX

genotype_24m <- factor(
  sample_list$Genotype,
  levels = c("WT", "NLF", "hTau", "NLFhTau")
)

design_24m_hk <- model.matrix(~ 0 + genotype_24m)
colnames(design_24m_hk) <- levels(genotype_24m)

cont.matrix_24m_hk <- makeContrasts(
  hTau_vs_WT      = hTau - WT,
  NLF_vs_WT       = NLF - WT,
  NLFhTau_vs_NLF  = NLFhTau - NLF,
  NLFhTau_vs_hTau = NLFhTau - hTau,
  NLFhTau_vs_WT   = NLFhTau - WT,
  levels = design_24m_hk
)

# VOOM + LIMMA
v_24m_syn_hk <- voom(dge_24m_syn, design_24m_hk, plot = TRUE)

fit_24m_syn_hk  <- lmFit(v_24m_syn_hk, design_24m_hk)
fit2_24m_syn_hk <- contrasts.fit(fit_24m_syn_hk, cont.matrix_24m_hk)
fit2_24m_syn_hk <- eBayes(fit2_24m_syn_hk)

# EXPORT
wb_syn <- createWorkbook()
results_list_syn <- list()

for (contrast_name in colnames(cont.matrix_24m_hk)) {
  
  results <- topTable(
    fit2_24m_syn_hk,
    coef          = contrast_name,
    number        = Inf,
    adjust.method = "BH",
    sort.by       = "none"
  )
  
  results$Ensembl_ID <- rownames(results)
  
  if ("ensembl_id" %in% colnames(results)) results$ensembl_id <- NULL
  if ("biotype" %in% colnames(results)) results$biotype <- NULL
  
  results$Status_BH <- "Not Significant"
  results$Status_BH[results$adj.P.Val < adjP_cutoff & results$logFC >  logFC_cutoff] <- "Upregulated"
  results$Status_BH[results$adj.P.Val < adjP_cutoff & results$logFC < -logFC_cutoff] <- "Downregulated"
  
  results$Status_BH <- factor(
    results$Status_BH,
    levels = c("Downregulated", "Upregulated", "Not Significant")
  )
  
  first_cols <- c("Ensembl_ID", "mgi_symbol")
  first_cols <- first_cols[first_cols %in% colnames(results)]
  other_cols <- setdiff(colnames(results), first_cols)
  results <- results[, c(first_cols, other_cols)]
  
  results_list_syn[[contrast_name]] <- results
  addWorksheet(wb_syn, contrast_name)
  writeData(wb_syn, contrast_name, results)
}

saveWorkbook(
  wb_syn,
  "RNAseq_DEG_24m_synaptic_all_contrasts_HKscaled.xlsx",
  overwrite = TRUE
)

# VOLCANO: NLFhTau vs NLF

results_volcano <- topTable(
  fit2_24m_syn_hk,
  coef = "NLFhTau_vs_NLF",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

results_volcano$Ensembl_ID <- rownames(results_volcano)

adjP_cutoff <- 0.05
label_logFC_cutoff <- 0

results_volcano <- results_volcano %>%
  mutate(
    adj.P.Val.plot = ifelse(adj.P.Val == 0, .Machine$double.xmin, adj.P.Val),
    negLog10AdjP = -log10(adj.P.Val.plot),
    Significance = case_when(
      adj.P.Val < adjP_cutoff & logFC < 0 ~ "Downregulated",
      adj.P.Val < adjP_cutoff & logFC > 0 ~ "Upregulated",
      TRUE ~ "Not Significant"
    ),
    Label = ifelse(
      adj.P.Val < adjP_cutoff & logFC < label_logFC_cutoff,
      mgi_symbol,
      NA
    )
  )

results_volcano$Significance <- factor(
  results_volcano$Significance,
  levels = c("Downregulated", "Upregulated", "Not Significant")
)

p_volcano <- ggplot(results_volcano, aes(x = logFC, y = negLog10AdjP)) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      "Downregulated" = "blue",
      "Upregulated" = "red",
      "Not Significant" = "grey75"
    )
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_text_repel(
    aes(label = Label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  theme_minimal(base_size = 10) +
  labs(
    title = "NLFhTau vs NLF",
    x = expression(Log[2]*"(fold change)"),
    y = expression(-Log[10]*"(adjusted P-value)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

print(p_volcano)

ggsave(
  filename = "Volcano_NLFhTau_vs_NLF_synaptic.pdf",
  plot = p_volcano,
  width = 6,
  height = 5
)

# 10. Export significant limma genes + expression values for Prism SYNAPSE ----
syn_contrast_name <- "NLFhTau_vs_NLF"
syn_output_file   <- "Prism_synaptic_NLFhTau_vs_NLF.xlsx"

# 1. Extract synapse limma results
syn_results_prism <- limma::topTable(
  fit2_24m_syn_hk,
  coef          = syn_contrast_name,
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "none"
) %>%
  tibble::rownames_to_column("Ensembl_ID")

# 2. Keep significant synaptic genes
syn_sig_genes <- syn_results_prism %>%
  dplyr::filter(adj.P.Val < adjP_cutoff)

# 3. Extract synapse expression matrix
syn_expr_mat <- v_24m_syn_hk$E

# 4. Keep genes present in expression matrix
syn_gene_ids_present <- intersect(
  syn_sig_genes$Ensembl_ID,
  rownames(syn_expr_mat)
)

cat(
  "Significant synaptic genes present in expression matrix:",
  length(syn_gene_ids_present),
  "\n"
)

syn_expr_sig <- syn_expr_mat[syn_gene_ids_present, , drop = FALSE]

# 5. Match limma table order to expression matrix order
syn_sig_genes <- syn_sig_genes %>%
  dplyr::filter(Ensembl_ID %in% syn_gene_ids_present) %>%
  dplyr::slice(match(rownames(syn_expr_sig), Ensembl_ID))

stopifnot(identical(syn_sig_genes$Ensembl_ID, rownames(syn_expr_sig)))

# 6. Build sample metadata
syn_sample_info <- data.frame(
  Sample_ID = colnames(syn_expr_sig),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    sample_list %>%
      dplyr::select(Sample_ID, Genotype),
    by = "Sample_ID"
  ) %>%
  dplyr::mutate(
    Genotype = factor(
      Genotype,
      levels = c("WT", "NLF", "hTau", "NLFhTau")
    )
  )

if (any(is.na(syn_sample_info$Genotype))) {
  stop("Some synapse samples could not be matched to sample_list$Sample_ID.")
}

# 7. Reorder samples and rename for Prism
syn_sample_info <- syn_sample_info %>%
  dplyr::arrange(Genotype, Sample_ID) %>%
  dplyr::mutate(
    Prism_Colname = paste(Genotype, Sample_ID, sep = "_")
  )

syn_expr_sig <- syn_expr_sig[, syn_sample_info$Sample_ID, drop = FALSE]
colnames(syn_expr_sig) <- syn_sample_info$Prism_Colname

# 8. Build export table
syn_export_table <- data.frame(
  Ensembl_ID = syn_sig_genes$Ensembl_ID,
  mgi_symbol = syn_sig_genes$mgi_symbol,
  logFC      = syn_sig_genes$logFC,
  adj.P.Val  = syn_sig_genes$adj.P.Val,
  syn_expr_sig,
  check.names = FALSE
) %>%
  dplyr::arrange(adj.P.Val)

# 9. Save Excel
syn_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(syn_wb, "SigGenes_PlotData")
openxlsx::writeData(syn_wb, "SigGenes_PlotData", syn_export_table)

openxlsx::saveWorkbook(
  syn_wb,
  file = syn_output_file,
  overwrite = TRUE
)

cat("Excel file saved:", syn_output_file, "\n")