

##### set up env
####################################################################################################

# clear env
rm(list = ls())

# load packages
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(extrafont)
loadfonts(device = "pdf")

par(family = "Times New Roman")

# load event data from all samples
events <- read.delim("../dataset/sample_matrices/sample_matrix_Pan.txt", sep = " ", row.names = 1)

# load vector of target genes
genes <- scan("../dataset/gene_panel.txt", what = character(), sep = ",", quiet = TRUE)

# load group 11 theta matrix
theta <- as.matrix(read.table("../classification_results/thetas_oMHN/theta_group11.dat", header = FALSE, sep = " ", stringsAsFactors = FALSE, check.names = FALSE))
colnames(theta) <- genes; rownames(theta) <- c(genes, "Observation")

# load final group memberships
mships <- readLines("../classification_results/classification_oMHN_13groups.dat")
names(mships) <- rownames(events)

# load sample metadata
samples <- read.delim("../cbioportal/data_clinical_sample.txt", skip = 4)

# load patient metadata
patients <- read.delim("../cbioportal/data_clinical_patient.txt", skip = 4)

####################################################################################################


##### plot group 11 theta
####################################################################################################

# extract base rates and mask them in theta matrix
theta_to_plot <- theta
base_rates <- c()
for (r in rownames(theta_to_plot)) {
  for (c in colnames(theta_to_plot)) {
    if (r == c) {base_rates <- c(base_rates, theta_to_plot[r, c]); theta_to_plot[r, c] <- NA}
  }
}
theta_to_plot <- as.matrix(theta_to_plot)

# make color schemes
cf_baserates = colorRamp2(c(floor(min(base_rates)), ceiling(max(base_rates))), c("white", "#0072B2"))
cf_theta = colorRamp2(c(floor(min(theta_to_plot, na.rm = T)), 0, ceiling(max(theta_to_plot, na.rm = T))), c("#E69F00", "white", "#009E73"))

# make legends
lgd_theta = Legend(col_fun = cf_theta, title = "Dependencies", at = c(floor(min(theta_to_plot, na.rm = T)), 0, ceiling(max(theta_to_plot, na.rm = T))), direction = "horizontal", legend_width = unit(5, "cm"), grid_height = unit(1.5, "cm"), border = "black", title_gp = gpar(fontfamily = "Times New Roman", cex = 1.5, fontface = "bold"))
lgd_baserates = Legend(col_fun = cf_baserates, title = "Base Rates", at = c(floor(min(base_rates)), ceiling(max(base_rates))), direction = "horizontal", legend_width = unit(5, "cm"), grid_height = unit(1.5, "cm"), border = "black", title_gp = gpar(fontfamily = "Times New Roman", cex = 1.5, fontface = "bold"))
packed_lgd = packLegend(lgd_theta, lgd_baserates, direction = "horizontal", column_gap = unit(2, "cm"))

# set up heatmap to plot
HM <- Heatmap(theta_to_plot, col = cf_theta,
              cluster_rows = F, cluster_columns = F,
              show_heatmap_legend = F,
              rect_gp = gpar(fontfamily = "Times New Roman", col = "lightgrey", lwd = 1),
              row_names_side = "left",
              row_names_gp = gpar(fontfamily = "Times New Roman", fontface = "bold", cex = 1.25),
              column_names_side = "top",
              column_names_rot = 45,
              column_names_gp = gpar(fontfamily = "Times New Roman", fontface = "bold", cex = 1.25),
              row_gap = unit(4, "mm"), column_gap = unit(4, "mm"),
              border = c("black"),
              row_title = NULL,
              column_title = NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (i != j) {
                  if (abs(theta_to_plot[i, j]) > 0.3) {
                    grid.text(sprintf("%.1f", theta_to_plot[i, j]), x, y, gp = gpar(fontfamily = "Times New Roman", fontsize = 10))
                  }
                }
              },
              left_annotation = rowAnnotation(BaseRates = c(base_rates, NA),
                                              col = list(BaseRates = cf_baserates),
                                              border = T,
                                              show_legend = F,
                                              annotation_label = c(""),
                                              annotation_name_gp = gpar(fontfamily = "Times New Roman", fontface = "bold", cex = 1.5),
                                              simple_anno_size=unit(1, "cm")))


# save as pdf
pdf(paste0("theta_group11.pdf"), w=16, h=16)

draw(HM, padding = unit(c(3, 2, 1, 2), "cm"))
draw(packed_lgd, x = unit(25, "cm"), y = unit(0.2, "cm"), just = c("right", "bottom"))

dev.off()

####################################################################################################


##### plot group 11 oncoplot
####################################################################################################

# get events from group 11
g11 <- mships[mships == 11]
g11_events <- events[names(g11), setdiff(colnames(events), "Altered")]

# unify sample ID format
names(g11) <- gsub("msk_chord_2024:", "", names(g11))
rownames(g11_events) <- gsub("msk_chord_2024:", "", rownames(g11_events))

# get samples' metadata
rownames(samples) <- samples$SAMPLE_ID
g11_samples <- samples[rownames(g11_events), ]

# check oncotree cancer type dist
# sort(table(g11_samples$CANCER_TYPE))
g11_samples$NSCLC <- as.factor(ifelse(g11_samples$CANCER_TYPE == "Non-Small Cell Lung Cancer", "NSCLC", "OTHERS"))
# table(g11_samples$NSCLC)

# sanity check
if(!identical(rownames(g11_samples), rownames(g11_events))) {stop("wrong sorting")}

# get 2-type factor
type <- droplevels(g11_samples[rownames(g11_events), "NSCLC"])

# get mutation prevalence by gene
gene_prev <- colSums(g11_events == 1, na.rm = TRUE)

# get the 10 events with the strongest interaction (absolute sum of both directions) with STK11
int_strengths <- c()
for (i in setdiff(colnames(theta), "STK11")) {
    is <- abs(theta[i, "STK11"]) + abs(theta["STK11", i])
    int_strengths <- c(int_strengths, is)
}
names(int_strengths) <- setdiff(colnames(theta), "STK11")

# get the prevalences of the 10 most interacting events
int_prev <- gene_prev[names(sort(int_strengths, decreasing = TRUE))[1:10]]

# sort them by prevalence, with STK11 first
gene_order <- c("STK11", names(sort(int_prev, decreasing = TRUE)))

# make color scheme
mut_col <- "#0072B2"
alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w, h, gp = gpar(fill = "#F2F2F2", col = NA)),
  Mut = function(x, y, w, h)
    grid.rect(x, y, w, h, gp = gpar(fill = mut_col, col = NA))
)

# oncoplot column order function (also works for a prior group split)
ord_within_subset <- function(pats, events_mat, gene_order) {
  E <- events_mat[pats, gene_order, drop = FALSE]  
  args <- lapply(as.data.frame(E), function(x) -as.integer(x))
  pats[do.call(order, args)]
}

# function to finalise OP for a specific group label and print it to pdf
save_oncoprint_one_group <- function(group_label, out_file, w = 8, h = 8) {
  pats <- rownames(g11_events)[type == group_label]
  pats <- ord_within_subset(pats, g11_events, gene_order)

  events_sub <- g11_events[pats, , drop = FALSE]

  # genes x patients (only the selected rows/genes, but burden barplot stays "ALL genes")
  mat01_sub  <- t(events_sub[, gene_order, drop = FALSE])
  mat_chr_sub <- ifelse(mat01_sub == 1, "Mut", "")

  ht_sub <- oncoPrint(
    mat_chr_sub,
    pct_gp = gpar(fontfamily = "Times New Roman"),
    alter_fun = alter_fun,
    col = c(Mut = mut_col),
    row_order = gene_order,
    column_order = pats,
    top_annotation = NULL,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(title = "Mutation", at = "Mut", labels = "Mut"),
    row_names_gp = gpar(fontfamily = "Times New Roman"),
    column_names_gp = gpar(fontfamily = "Times New Roman"),
  )

  pdf(out_file, width = w, height = h, useDingbats = FALSE)
  draw(ht_sub, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# write to pdfs of same size
save_oncoprint_one_group(
  group_label = "NSCLC",
  out_file = "oncoprint_group11_NSCLC.pdf",
  w = 4, h = 4
)

save_oncoprint_one_group(
  group_label = "OTHERS",
  out_file = "oncoprint_group11_OTHERS.pdf",
  w = 4, h = 4
)

####################################################################################################
