# 1. Packages ----

library(readxl)
library(openxlsx)
library(tibble)
library(dplyr)
library(stringr)
library(ggplot2)
library(fgsea)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

library(clusterProfiler)
library(GOSemSim)
library(GO.db)
library(tidyr)
library(ggrepel)

# 2. Settings ----

human_cut_prot <- 0.05
mouse_cut_prot <- 0.10

clean_uniprot <- function(x) {
  x <- trimws(as.character(x))
  x <- ifelse(grepl("\\|", x), sub("^.*\\|([^|]+)\\|.*$", "\\1", x), x)
  x <- sub("-\\d+$", "", x)
  x <- sub("[ \\(].*$", "", x)
  x
}

# 3. Human proteomics: UniProt to HUMAN SYMBOL ----

human_prot <- readxl::read_excel(
  "human_and_mouse_PROTEOMICS.xlsx",
  sheet = "human_PROT"
) %>%
  dplyr::transmute(
    accession_human   = as.character(accession_human),
    uniprot_human     = clean_uniprot(accession_human),
    log2FC_human_prot = as.numeric(log2FC_human_prot),
    adjP_human_prot   = as.numeric(adjP_human_prot)
  )

hum_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(na.omit(human_prot$uniprot_human)),
  keytype = "UNIPROT",
  columns = "SYMBOL"
) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    uniprot_human = UNIPROT,
    gene_symbol_human = SYMBOL
  ) %>%
  dplyr::distinct(uniprot_human, .keep_all = TRUE)

human_prot <- human_prot %>%
  dplyr::left_join(hum_map, by = "uniprot_human")

cat(
  "Human mapped to symbols:",
  sum(!is.na(human_prot$gene_symbol_human)),
  "/",
  nrow(human_prot),
  "\n"
)

# 4. Mouse proteomics: mouse symbol to HUMAN SYMBOL ----

mouse_prot <- readxl::read_excel(
  "human_and_mouse_PROTEOMICS.xlsx",
  sheet = "mouse_PROT"
) %>%
  dplyr::transmute(
    mouse_symbol_raw  = stringr::str_trim(as.character(mouse_protein_name)),
    log2FC_mouse_prot = as.numeric(log2FC_mouse_prot),
    adjP_mouse_prot   = as.numeric(adjP_mouse_prot)
  ) %>%
  dplyr::mutate(
    sym1 = mouse_symbol_raw,
    sym2 = stringr::str_to_title(stringr::str_to_lower(mouse_symbol_raw)),
    sym3 = stringr::str_to_upper(mouse_symbol_raw)
  )

map_symbol_to_ens <- function(keys) {
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys      = unique(na.omit(keys)),
    keytype   = "SYMBOL",
    column    = "ENSEMBL",
    multiVals = "first"
  )
}

m1 <- map_symbol_to_ens(mouse_prot$sym1)
m2 <- map_symbol_to_ens(mouse_prot$sym2)
m3 <- map_symbol_to_ens(mouse_prot$sym3)

mouse_prot <- mouse_prot %>%
  dplyr::mutate(
    ensembl_mouse = unname(m1[sym1]),
    ensembl_mouse = ifelse(is.na(ensembl_mouse), unname(m2[sym2]), ensembl_mouse),
    ensembl_mouse = ifelse(is.na(ensembl_mouse), unname(m3[sym3]), ensembl_mouse)
  ) %>%
  dplyr::select(-sym1, -sym2, -sym3)

cat(
  "Mouse symbols mapped to mouse Ensembl:",
  sum(!is.na(mouse_prot$ensembl_mouse)),
  "/",
  nrow(mouse_prot),
  "\n"
)

orth <- readxl::read_excel("mouse_human_orthologues_ensembl115.xls") %>%
  dplyr::filter(`Human homology type` == "ortholog_one2one") %>%
  dplyr::transmute(
    ensembl_mouse = `Gene stable ID`,
    ensembl_human = `Human gene stable ID`
  )

mouse_prot <- mouse_prot %>%
  dplyr::left_join(orth, by = "ensembl_mouse")

hum_ens2sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = unique(na.omit(mouse_prot$ensembl_human)),
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)

mouse_prot <- mouse_prot %>%
  dplyr::mutate(
    gene_symbol_human = unname(hum_ens2sym[ensembl_human])
  )

cat(
  "Mouse mapped to HUMAN symbols:",
  sum(!is.na(mouse_prot$gene_symbol_human)),
  "/",
  nrow(mouse_prot),
  "\n"
)

# 5. Join mouse and human proteomics by HUMAN SYMBOL ----

pairs_prot <- mouse_prot %>%
  dplyr::filter(!is.na(gene_symbol_human), gene_symbol_human != "") %>%
  dplyr::inner_join(
    human_prot %>%
      dplyr::filter(!is.na(gene_symbol_human), gene_symbol_human != "") %>%
      dplyr::select(
        gene_symbol_human,
        log2FC_human_prot,
        adjP_human_prot
      ),
    by = "gene_symbol_human"
  )

cat("Matched proteins mouse-human:", nrow(pairs_prot), "\n")

# 6. Create human-ranked gene list ----

rank_tbl_prot <- pairs_prot %>%
  dplyr::filter(
    !is.na(gene_symbol_human),
    gene_symbol_human != "",
    is.finite(log2FC_human_prot),
    !is.na(adjP_human_prot),
    adjP_human_prot > 0
  ) %>%
  dplyr::mutate(
    rank_metric = log2FC_human_prot * -log10(adjP_human_prot)
  ) %>%
  dplyr::arrange(
    dplyr::desc(rank_metric),
    dplyr::desc(abs(log2FC_human_prot))
  ) %>%
  dplyr::distinct(gene_symbol_human, .keep_all = TRUE)

gene_list_human_prot <- rank_tbl_prot$rank_metric
names(gene_list_human_prot) <- rank_tbl_prot$gene_symbol_human
gene_list_human_prot <- sort(gene_list_human_prot, decreasing = TRUE)

# 7. Create mouse-defined gene sets ----

mouse_up_prot <- pairs_prot %>%
  dplyr::filter(adjP_mouse_prot < mouse_cut_prot, log2FC_mouse_prot > 0) %>%
  dplyr::pull(gene_symbol_human) %>%
  na.omit() %>%
  unique()

mouse_down_prot <- pairs_prot %>%
  dplyr::filter(adjP_mouse_prot < mouse_cut_prot, log2FC_mouse_prot < 0) %>%
  dplyr::pull(gene_symbol_human) %>%
  na.omit() %>%
  unique()

mouse_up_human_up_prot <- pairs_prot %>%
  dplyr::filter(
    adjP_mouse_prot < mouse_cut_prot,
    log2FC_mouse_prot > 0,
    adjP_human_prot < human_cut_prot,
    log2FC_human_prot > 0
  ) %>%
  dplyr::pull(gene_symbol_human) %>%
  na.omit() %>%
  unique()

mouse_down_human_down_prot <- pairs_prot %>%
  dplyr::filter(
    adjP_mouse_prot < mouse_cut_prot,
    log2FC_mouse_prot < 0,
    adjP_human_prot < human_cut_prot,
    log2FC_human_prot < 0
  ) %>%
  dplyr::pull(gene_symbol_human) %>%
  na.omit() %>%
  unique()

pathways_mouse_prot <- list(
  Mouse_DEG_Up        = mouse_up_prot,
  Mouse_DEG_Down      = mouse_down_prot,
  MouseUp_HumanUp     = mouse_up_human_up_prot,
  MouseDown_HumanDown = mouse_down_human_down_prot
)

# 8. Run GSEA ----

fg_prot <- fgsea::fgsea(
  pathways = pathways_mouse_prot,
  stats    = gene_list_human_prot,
  minSize  = 15,
  maxSize  = 5000,
  nperm    = 10000
) %>%
  dplyr::arrange(NES)

fg_prot

# 9. Save GSEA results ----

fg_prot_export <- fg_prot %>%
  dplyr::mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = "; ")
  )

wb_gsea <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_gsea, "GSEA_results")
openxlsx::writeData(wb_gsea, "GSEA_results", fg_prot_export)

openxlsx::saveWorkbook(
  wb_gsea,
  "GSEA_Proteomics_results.xlsx",
  overwrite = TRUE
)

# 10. Plot selected GSEA pathway ----

geneset_name <- "Mouse_DEG_Up"

geneset <- pathways_mouse_prot[[geneset_name]]

p0 <- fgsea::plotEnrichment(
  geneset,
  stats = gene_list_human_prot
)

geom_names <- vapply(
  p0$layers,
  function(x) class(x$geom)[1],
  character(1)
)

i_curve <- which(geom_names == "GeomLine")[1]
i_tick  <- which(geom_names == "GeomSegment")[1]

curve_df <- ggplot2::layer_data(p0, i_curve)
ticks_df <- ggplot2::layer_data(p0, i_tick)

peak_x <- curve_df$x[which.max(abs(curve_df$y))]

tick_ymin <- min(c(ticks_df$y, ticks_df$yend), na.rm = TRUE)
tick_ymax <- max(c(ticks_df$y, ticks_df$yend), na.rm = TRUE)

p_gsea <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = tick_ymin,
    ymax = tick_ymax,
    fill = "#0F80FF",
    alpha = 0.12
  ) +
  ggplot2::geom_segment(
    data = ticks_df,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    color = "#0F80FF",
    linewidth = 0.6,
    lineend = "butt"
  ) +
  ggplot2::geom_line(
    data = curve_df,
    ggplot2::aes(x = x, y = y),
    color = "#0F80FF",
    linewidth = 1.05
  ) +
  ggplot2::geom_vline(
    xintercept = peak_x,
    color = "#0F80FF",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  ggplot2::coord_cartesian(ylim = c(-0.1, 1), expand = FALSE) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    title = paste("GSEA:", geneset_name, "vs Human Proteomics Ranking"),
    x = expression("Human proteins ranked by " * log[2] * "FC × " * -log[10] * "(adjusted p-value)"),
    y = "Enrichment Score"
  )

ggplot2::ggsave(
  filename = paste0("GSEA_Proteomics_", geneset_name, ".pdf"),
  plot     = p_gsea,
  width    = 6,
  height   = 5,
  units    = "in",
  device   = if (capabilities("cairo")) cairo_pdf else "pdf"
)


# 11. Show NES, adjusted p-value, and leading-edge genes ----

gsea_summary <- fg_prot %>%
  dplyr::filter(pathway == geneset_name) %>%
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    size,
    leadingEdge
  )

gsea_summary

leading_edge_genes <- gsea_summary$leadingEdge[[1]]

leading_edge_genes

# 12. Heatmap: Shared DEPs with concordance ----

# Build table
hm_shared_tbl <- pairs_prot %>%
  dplyr::filter(
    !is.na(gene_symbol_human),
    gene_symbol_human != "",
    is.finite(log2FC_human_prot),
    is.finite(log2FC_mouse_prot),
    adjP_human_prot < human_cut_prot,
    adjP_mouse_prot < mouse_cut_prot
  ) %>%
  dplyr::transmute(
    Gene    = gene_symbol_human,
    Human   = log2FC_human_prot,
    NLFhTau = log2FC_mouse_prot,
    Concordance = dplyr::case_when(
      Human > 0 & NLFhTau > 0 ~ "Up-Up",
      Human < 0 & NLFhTau < 0 ~ "Down-Down",
      TRUE                    ~ "Discordant"
    ),
    score = abs(Human) + abs(NLFhTau)
  ) %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::arrange(
    factor(Concordance, levels = c("Up-Up", "Down-Down", "Discordant")),
    dplyr::desc(score)
  )

cat("Shared DEPs:", nrow(hm_shared_tbl), "\n")

# Matrix
hm_mat <- rbind(
  Human   = hm_shared_tbl$Human,
  NLFhTau = hm_shared_tbl$NLFhTau
)

colnames(hm_mat) <- hm_shared_tbl$Gene

# Concordance annotation
hm_top_ha <- ComplexHeatmap::HeatmapAnnotation(
  Concordance = hm_shared_tbl$Concordance,
  col = list(
    Concordance = c(
      "Up-Up"      = "#009E73",
      "Down-Down"  = "#E69F00",
      "Discordant" = "grey60"
    )
  ),
  show_legend = TRUE
)

# Colour scale
hm_lim <- max(abs(hm_mat), na.rm = TRUE)

hm_col_fun <- circlize::colorRamp2(
  c(-hm_lim, 0, hm_lim),
  c("#5A7DFF", "white", "#E35D5D")
)

# Save
pdf(
  "Heatmap_sharedDEPs.pdf",
  width  = max(12, 0.25 * ncol(hm_mat) + 4),
  height = 3.2,
  useDingbats = FALSE
)

ComplexHeatmap::draw(
  ComplexHeatmap::Heatmap(
    hm_mat,
    name = "log2FC",
    col = hm_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = hm_top_ha
  )
)

dev.off()
# 13. GO enrichment for human and mouse significant proteins ----
# Human significant proteins
human_sig_proteins <- human_prot %>%
  dplyr::filter(adjP_human_prot < human_cut_prot) %>%
  dplyr::filter(!is.na(gene_symbol_human), gene_symbol_human != "") %>%
  dplyr::pull(gene_symbol_human) %>%
  unique()

human_entrez <- clusterProfiler::bitr(
  human_sig_proteins,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Mouse significant proteins
mouse_sig_proteins <- mouse_prot %>%
  dplyr::filter(adjP_mouse_prot < mouse_cut_prot) %>%
  dplyr::mutate(
    mouse_symbol = stringr::str_to_title(stringr::str_to_lower(mouse_symbol_raw))
  ) %>%
  dplyr::pull(mouse_symbol) %>%
  unique()

mouse_entrez <- clusterProfiler::bitr(
  mouse_sig_proteins,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

run_go_bp <- function(entrez_ids, OrgDb) {
  clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  ) %>%
    as.data.frame() %>%
    dplyr::filter(Count >= 5) %>%
    dplyr::arrange(p.adjust)
}

human_go_df <- run_go_bp(human_entrez$ENTREZID, org.Hs.eg.db) %>%
  dplyr::mutate(Species = "Human")

mouse_go_df <- run_go_bp(mouse_entrez$ENTREZID, org.Mm.eg.db) %>%
  dplyr::mutate(Species = "Mouse")

go_tbl <- dplyr::bind_rows(human_go_df, mouse_go_df) %>%
  dplyr::rename(GO_ID = ID)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "GO_BP_Human_Mouse")
openxlsx::writeData(wb, "GO_BP_Human_Mouse", go_tbl)

openxlsx::saveWorkbook(
  wb,
  "GO_enrichment_mouse_human.xlsx",
  overwrite = TRUE
)


# 14. Cluster similar mouse GO terms ----

mouse_tbl <- go_tbl %>%
  dplyr::filter(Species == "Mouse") %>%
  dplyr::distinct(GO_ID, .keep_all = TRUE) %>%
  dplyr::select(
    Mouse_GO = GO_ID,
    Mouse_Description = Description,
    Mouse_adjP = p.adjust
  )

mouse_go <- mouse_tbl$Mouse_GO

hsGO <- GOSemSim::godata(
  annoDb = "org.Hs.eg.db",
  ont = "BP",
  computeIC = FALSE
)

mouse_sim_mat <- GOSemSim::mgoSim(
  GO1 = mouse_go,
  GO2 = mouse_go,
  semData = hsGO,
  measure = "Wang",
  combine = NULL
)

mouse_dist <- stats::as.dist(1 - mouse_sim_mat)

hc <- stats::hclust(
  mouse_dist,
  method = "average"
)

cluster_assignments <- stats::cutree(
  hc,
  h = 0.4
)

cluster_df <- tibble::tibble(
  Mouse_GO = names(cluster_assignments),
  Cluster_ID = paste0("C", cluster_assignments)
) %>%
  dplyr::left_join(mouse_tbl, by = "Mouse_GO")

cluster_summary <- cluster_df %>%
  dplyr::group_by(Cluster_ID) %>%
  dplyr::summarise(
    n_terms = dplyr::n(),
    min_mouse_adjP = min(Mouse_adjP, na.rm = TRUE),
    mouse_terms = paste(Mouse_GO, collapse = "; "),
    mouse_descriptions = paste(Mouse_Description, collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Cluster_ID)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "mouse_clusters")
openxlsx::writeData(wb, "mouse_clusters", cluster_summary)

openxlsx::saveWorkbook(
  wb,
  "Mouse_GO_Clusters_for_Curation.xlsx",
  overwrite = TRUE
)

cat("Saved: Mouse_GO_Clusters_for_Curation.xlsx\n")
# 17. Bubble plot: curated mouse GO clusters vs human GO terms ----

# Human GO table from existing go_tbl
human_tbl <- go_tbl %>%
  dplyr::filter(Species == "Human") %>%
  dplyr::distinct(GO_ID, .keep_all = TRUE) %>%
  dplyr::select(
    Human_GO = GO_ID,
    Human_Description = Description,
    Human_adjP = p.adjust
  )

human_go <- unique(human_tbl$Human_GO)

# Mouse GO table from existing go_tbl
mouse_tbl <- go_tbl %>%
  dplyr::filter(Species == "Mouse") %>%
  dplyr::distinct(GO_ID, .keep_all = TRUE) %>%
  dplyr::select(
    Mouse_GO = GO_ID,
    Mouse_Description = Description,
    Mouse_adjP = p.adjust
  )

# Read curated clusters
cluster_curated <- readxl::read_excel("Mouse_GO_Clusters_Curated.xlsx") %>%
  dplyr::select(
    Cluster_ID,
    n_terms,
    mouse_terms,
    mouse_descriptions,
    bubble_title
  ) %>%
  dplyr::filter(
    !is.na(Cluster_ID),
    !is.na(mouse_terms),
    !is.na(bubble_title),
    bubble_title != ""
  )

# Expand mouse terms
cluster_members <- cluster_curated %>%
  dplyr::mutate(mouse_terms = stringr::str_split(mouse_terms, ";\\s*")) %>%
  tidyr::unnest(mouse_terms) %>%
  dplyr::rename(Mouse_GO = mouse_terms) %>%
  dplyr::mutate(Mouse_GO = stringr::str_trim(Mouse_GO)) %>%
  dplyr::filter(Mouse_GO != "") %>%
  dplyr::left_join(mouse_tbl, by = "Mouse_GO")

# Compare each mouse GO term to all human GO terms
mouse_human_sim_mat <- GOSemSim::mgoSim(
  GO1 = unique(cluster_members$Mouse_GO),
  GO2 = human_go,
  semData = hsGO,
  measure = "Wang",
  combine = NULL
)

mouse_human_df <- as.data.frame(mouse_human_sim_mat) %>%
  tibble::rownames_to_column("Mouse_GO") %>%
  tidyr::pivot_longer(
    cols = -Mouse_GO,
    names_to = "Human_GO",
    values_to = "Wang"
  ) %>%
  dplyr::filter(!is.na(Wang)) %>%
  dplyr::left_join(human_tbl, by = "Human_GO")

# Keep best human match for each mouse GO term
best_human_match <- mouse_human_df %>%
  dplyr::filter(Mouse_GO != Human_GO) %>%
  dplyr::group_by(Mouse_GO) %>%
  dplyr::slice_max(order_by = Wang, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Plot
p_bubble <- ggplot2::ggplot(
  plot_df_plot,
  ggplot2::aes(
    x = mean_best_human_Wang,
    y = neglog10_min_mouse_adjP,
    size = n_terms
  )
) +
  ggplot2::geom_point(
    shape = 21,
    fill = "#BFDDF2",
    colour = "black",
    stroke = 0.5,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = bubble_label),
    size = 3.6,
    box.padding = 0.8,
    point.padding = 1.0,
    force = 3,
    force_pull = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey65",
    segment.size = 0.25,
    seed = 123
  ) +
  ggplot2::scale_size_continuous(
    range = c(4, 13),
    breaks = unique(sort(plot_df_plot$n_terms)),
    labels = function(x) as.integer(x),
    name = "Number of\nGO terms"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0.70, 1.00),
    breaks = seq(0.70, 1.00, by = 0.05),
    expand = ggplot2::expansion(mult = c(0.03, 0.08))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(1.6, 2.1),
    breaks = seq(1.6, 2.1, by = 0.1),
    expand = ggplot2::expansion(mult = c(0.05, 0.08))
  ) +
  ggplot2::labs(
      x = "Similarity between mouse 
        and human GO terms (Wang rating)",
      y = expression(-log[10]("adjusted p-value"))
    ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold", size = 12),
    axis.text = ggplot2::element_text(size = 12, colour = "black"),
    legend.title = ggplot2::element_text(face = "bold", size = 11),
    legend.text = ggplot2::element_text(size = 10),
    legend.position = c(0.02, 0.02),  # bottom-left (inside panel)
    legend.justification = c(0, 0),
    plot.margin = ggplot2::margin(25, 80, 25, 35)
  )
ggplot2::ggsave(
  filename = "Mouse_GO_Cluster_BubblePlot_vs_Human_GO.pdf",
  plot = p_bubble,
  width = 10,
  height = 8,
  units = "in",
  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
  bg = "white"
)