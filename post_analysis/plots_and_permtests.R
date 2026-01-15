

##### set up env
####################################################################################################

# clear env
rm(list = ls())

# load packages
library(ComplexHeatmap)
library(circlize)

# load event data from all samples
events <- read.delim("dataset/sample_matrices/sample_matrix_Pan.txt", sep = " ", row.names = 1)

# load vector of target genes
genes <- scan("dataset/gene_panel.txt", what = character(), sep = ",", quiet = TRUE)

# load group 11 theta matrix
theta <- as.matrix(read.table("classification_results/thetas_oMHN/theta_group11.dat", header = FALSE, sep = " ", stringsAsFactors = FALSE, check.names = FALSE))
colnames(theta) <- genes; rownames(theta) <- c(genes, "Observation")

# load final group memberships
mships <- readLines("classification_results/classification_oMHN_13groups.dat")
names(mships) <- rownames(events)

# load sample metadata
samples <- read.delim("dataset/data_clinical_sample.txt", skip = 4)

# load patient metadata
patients <- read.delim("dataset/data_clinical_patient.txt", skip = 4)

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
lgd_theta = Legend(col_fun = cf_theta, title = "Dependencies", at = c(floor(min(theta_to_plot, na.rm = T)), 0, ceiling(max(theta_to_plot, na.rm = T))), direction = "horizontal", legend_width = unit(5, "cm"), grid_height = unit(1.5, "cm"), border = "black", title_gp = gpar(cex = 1.5, fontface = "bold"))
lgd_baserates = Legend(col_fun = cf_baserates, title = "Base Rates", at = c(floor(min(base_rates)), ceiling(max(base_rates))), direction = "horizontal", legend_width = unit(5, "cm"), grid_height = unit(1.5, "cm"), border = "black", title_gp = gpar(cex = 1.5, fontface = "bold"))
packed_lgd = packLegend(lgd_theta, lgd_baserates, direction = "horizontal", column_gap = unit(2, "cm"))

# set up heatmap to plot
HM <- Heatmap(theta_to_plot, col = cf_theta,
              cluster_rows = F, cluster_columns = F,
              show_heatmap_legend = F,
              rect_gp = gpar(col = "lightgrey", lwd = 1),
              row_names_side = "left",
              row_names_gp = gpar(fontface = "bold", cex = 1.25),
              column_names_side = "top",
              column_names_rot = 45,
              column_names_gp = gpar(fontface = "bold", cex = 1.25),
              row_gap = unit(4, "mm"), column_gap = unit(4, "mm"),
              border = c("black"),
              row_title = NULL,
              column_title = NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (i != j) {
                  if (abs(theta_to_plot[i, j]) > 0.3) {
                    grid.text(sprintf("%.1f", theta_to_plot[i, j]), x, y, gp = gpar(fontsize = 10))
                  }
                }
              },
              left_annotation = rowAnnotation(BaseRates = c(base_rates, NA),
                                              col = list(BaseRates = cf_baserates),
                                              border = T,
                                              show_legend = F,
                                              annotation_label = c(""),
                                              annotation_name_gp = gpar(fontface = "bold", cex = 1.5),
                                              simple_anno_size=unit(1, "cm")))


# save as pdf
pdf(paste0("post_analysis/theta_group11.pdf"), w=16, h=16)

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
    alter_fun = alter_fun,
    col = c(Mut = mut_col),
    row_order = gene_order,
    column_order = pats,
    top_annotation = NULL,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(title = "Mutation", at = "Mut", labels = "Mut")
  )

  pdf(out_file, width = w, height = h, useDingbats = FALSE)
  draw(ht_sub, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# write to pdfs of same size
save_oncoprint_one_group(
  group_label = "NSCLC",
  out_file = "post_analysis/oncoprint_group11_NSCLC.pdf",
  w = 4, h = 4
)

save_oncoprint_one_group(
  group_label = "OTHERS",
  out_file = "post_analysis/oncoprint_group11_OTHERS.pdf",
  w = 4, h = 4
)

####################################################################################################



##### prep data for permutation testing
####################################################################################################

# get events without "Altered" col from cbioportal
events <- events[, setdiff(colnames(events), "Altered")]

# clean sample IDs
rownames(events) <- gsub("msk_chord_2024:", "", rownames(events))

# get sample metadata for all PTs used in analysis
samples <- samples[rownames(events), ]
if (any(duplicated(samples$PATIENT_ID))) {stop("we exclude patients with >1 sample")}

# get the patient metadata matching these samples (we have strictly 1 sample per patient)
rownames(patients) <- patients$PATIENT_ID
patients <- patients[samples$PATIENT_ID, ]

# sanity check
if (!identical(samples$PATIENT_ID, rownames(patients))) {stop("wrong order")}

# get cols of interest from both per-sample and per-patient metadata - here I chose essentially all applicable variables that have no or only few NAs
perm_data <- cbind(samples[, c("SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "TMB_NONSYNONYMOUS", "CLINICAL_SUMMARY")], patients[, c("GENDER", "CURRENT_AGE_DEID", "STAGE_HIGHEST_RECORDED", "ADRENAL_GLANDS", "BONE", "CNS_BRAIN", "INTRA_ABDOMINAL", "LIVER", "LUNG", "LYMPH_NODES", "PLEURA", "REPRODUCTIVE_ORGANS", "SMOKING_PREDICTIONS_3_CLASSES")])

# replace unknowns with actual NAs
perm_data[perm_data == "Unknown"] <- NA
perm_data[perm_data == "Do not report"] <- NA
perm_data[perm_data == "Other"] <- NA

# recode clinical summary
perm_data$CLINICAL_SUMMARY[which(perm_data$CLINICAL_SUMMARY %in% c("Unknown/Unstaged", "Unstaged Unknown"))] <- NA
perm_data$CLINICAL_SUMMARY[which(perm_data$CLINICAL_SUMMARY %in% c("In Situ", "Localized"))] <- "localized"
perm_data$CLINICAL_SUMMARY[which(perm_data$CLINICAL_SUMMARY %in% c("Regional Both 2 and 3", "Regional by Direct Extension", "Regional Nos", "Regional To Lymph Nodes", "Regional,Direct Extension", "Regional,Extension and Nodes", "Regional,Lymph Nodes Only"))] <- "regional"
perm_data$CLINICAL_SUMMARY[which(perm_data$CLINICAL_SUMMARY %in% c("Distant", "Distant Metastases/Systemic Disease"))] <- "distant"

# remove entries with any NA
perm_data <- perm_data[which(rowSums(is.na(perm_data)) == 0), ]
str(perm_data)

# clean membership sample IDs
names(mships) <- gsub("msk_chord_2024:", "", names(mships))

# add to data as group 11 or others
perm_data_mships <- mships[perm_data$SAMPLE_ID]
perm_data$GROUP <- ifelse(perm_data_mships == "11", "11", "OTHERS")


####################################################################################################



##### define permutation testing function
####################################################################################################

# function summary: 
# for each applicable variable, the hypothesis is that group 11 has a more extreme distribution
# than randomly sampled groups from outside group 11 that have the same size and cancer type composition.
# for that, we first obtain B randomly sampled groups that fulfil the criterion
# then, for numerical vars we get all B means, for categorical vars we get the dist over labels
# then, for numerical vars we compare the dist over means with the observed mean to get p values - the difference in means could be interpreted as effect size
# for cat vars, we test the enrichment of each individual label in group 11 against all others (for 2 labels only 1 test is performed)
# the two-sided p-value calc should follow the formula in https://cran.r-project.org/web/packages/exactRankTests/exactRankTests.pdf (p. 9)

compare_group11_to_matched_others <- function(
  df, # prepped data
  x, # cols to test
  group_col  = "GROUP", # which col splits the data
  cancer_col = "CANCER_TYPE", # which cols distribution should be used in permutation
  group11_value = "11", # which value of split to test on one side
  others_value  = "OTHERS", # which value of split to test on the other side
  B = 10000L, # how many perms
  min_prev = 0.10, # all labels below this prevalence will be merged into "Other" group
  seed = NULL # set seed
) {

  # perform some input checks
  stopifnot(is.data.frame(df))
  stopifnot(is.character(x), length(x) >= 1)
  stopifnot(group_col %in% names(df), cancer_col %in% names(df))
  stopifnot(all(x %in% names(df)))
  stopifnot(is.numeric(B), length(B) == 1, B >= 1)

  # get group and type vectors
  grp    <- as.character(df[[group_col]])
  cancer <- as.character(df[[cancer_col]])
  
  # get group indices
  idx11 <- which(grp == group11_value)
  idxO  <- which(grp == others_value)

  # counts per cancer type (based on GROUP==11)
  n11_by_type <- table(cancer[idx11])
  types <- names(n11_by_type)
  n11_total <- sum(as.integer(n11_by_type))

  # Check OTHERS has sufficient rows in each required cancer type
  nO_by_type <- table(cancer[idxO])
  missing_types <- setdiff(types, names(nO_by_type))
  if (length(missing_types) > 0) {
    stop("OTHERS is missing required cancer types present in GROUP==11: ",
         paste(missing_types, collapse = ", "))
  }
  insufficient <- types[as.integer(nO_by_type[types]) < as.integer(n11_by_type[types])]
  if (length(insufficient) > 0) {
    stop("Insufficient OTHERS rows for exact matching without replacement in: ",
         paste(insufficient, collapse = ", "))
  }

  # get indices in OTHERS by cancer type
  others_idx_by_type <- lapply(types, function(t) idxO[cancer[idxO] == t])
  names(others_idx_by_type) <- types

  # get segment positions within the matched cohort vector
  counts <- as.integer(n11_by_type[types])
  ends <- cumsum(counts)
  starts <- c(1L, head(ends, -1L) + 1L)

  # set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # build an index matrix of matched OTHERS samples: n11_total x B
  null_idx <- matrix(NA_integer_, nrow = n11_total, ncol = as.integer(B))
  for (b in seq_len(B)) {
    for (i in seq_along(types)) {
      t <- types[i]
      k <- counts[i]
      seg <- starts[i]:ends[i]
      null_idx[seg, b] <- sample(others_idx_by_type[[t]], size = k, replace = FALSE)
    }
  }

  # define helper func for mpirical two-sided p-value centered at expected (mean of null)
  emp_p <- function(obs, null_vec) {
    expv <- mean(null_vec)
    (1 + sum(abs(null_vec - expv) >= abs(obs - expv))) / (length(null_vec) + 1)
  }

  # init stuff for reporting
  num_out <- list()
  cat_out <- list()
  ni <- 1L
  ci <- 1L

  # iterate over vars of interest
  for (coi in x) {

    v <- df[[coi]]

    # if numeric
    if (is.numeric(v)) {

      # get observed mean
      obs_mean <- mean(v[idx11])
      
      # make matrix of non-group-11 cases x permutations for var of interest
      v_null <- matrix(v[as.vector(null_idx)], nrow = n11_total, ncol = as.integer(B))

      # get permutation means
      null_means <- colMeans(v_null)

      # mean of means (expected value)
      exp_mean <- mean(null_means)

      # perform p value calc and report
      num_out[[ni]] <- data.frame(
        COI      = coi,
        n11      = n11_total,
        obs_mean = obs_mean,
        exp_mean = exp_mean,
        effect   = obs_mean - exp_mean,
        p        = emp_p(obs_mean, null_means),
        stringsAsFactors = FALSE
      )
      ni <- ni + 1L

    } else {
      f <- as.factor(v)

      # collapse rare levels based on overall prevalence
      lvl_tab <- table(f)
      lvl_prev <- as.numeric(lvl_tab) / sum(lvl_tab)
      rare_lvls <- names(lvl_tab)[lvl_prev < min_prev]

      f2 <- as.character(f)
      if (length(rare_lvls) > 0) f2[f2 %in% rare_lvls] <- "Other"
      f2 <- factor(f2)
      lvls <- levels(f2)

      # if no variability, nothing to test
      if (length(lvls) < 2) next

      # integer codes for fast counting
      codes <- as.integer(f2)
      codes_null <- matrix(codes[as.vector(null_idx)], nrow = n11_total, ncol = as.integer(B))

      # if exactly 2 levels, test only one (the other is redundant)
      levels_to_test <- lvls
      if (length(lvls) == 2) {
        if ("Other" %in% lvls) {
          levels_to_test <- setdiff(lvls, "Other")
          if (length(levels_to_test) == 0) levels_to_test <- lvls[1]
        } else {
          levels_to_test <- lvls[1]  
        }
      }

      # if more than 2 levels, do one test for each
      for (L in levels_to_test) {
        k <- match(L, lvls)

        obs_count <- sum(codes[idx11] == k)
        null_counts <- colSums(codes_null == k)
        exp_count <- mean(null_counts)

        cat_out[[ci]] <- data.frame(
          COI          = coi,
          level        = L,
          n11          = n11_total,
          obs_count    = obs_count,
          exp_count    = exp_count,
          effect_count = obs_count - exp_count,
          effect_prop  = (obs_count / n11_total) - (exp_count / n11_total),
          p            = emp_p(obs_count, null_counts),
          stringsAsFactors = FALSE
        )
        ci <- ci + 1L
      }
    }
  }

  # pack stuff for reporting, dependent on if num/cat
  numeric_df <- if (length(num_out)) do.call(rbind, num_out) else
    data.frame(COI=character(), n11=numeric(), obs_mean=numeric(), exp_mean=numeric(),
               effect=numeric(), p=numeric(), stringsAsFactors = FALSE)

  categorical_df <- if (length(cat_out)) do.call(rbind, cat_out) else
    data.frame(COI=character(), level=character(), n11=numeric(), obs_count=numeric(),
               exp_count=numeric(), effect_count=numeric(), effect_prop=numeric(),
               p=numeric(), stringsAsFactors = FALSE)

  # enforce numeric columns
  for (cc in c("n11","obs_mean","exp_mean","effect","p")) {
    if (cc %in% names(numeric_df)) numeric_df[[cc]] <- as.numeric(numeric_df[[cc]])
  }
  for (cc in c("n11","obs_count","exp_count","effect_count","effect_prop","p")) {
    if (cc %in% names(categorical_df)) categorical_df[[cc]] <- as.numeric(categorical_df[[cc]])
  }

  list(numeric = numeric_df, categorical = categorical_df)
}


####################################################################################################



##### execute permutation testing and write results
####################################################################################################

# for the entire group 11 - 100k perms (takes < 20 sec on my macbook)
res <- compare_group11_to_matched_others(df = perm_data, x = setdiff(colnames(perm_data), c("SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "GROUP")), seed = 7, B = 100000L)

# subset group11 for non NSCLC only
perm_data_nonlung <- perm_data[which(perm_data$CANCER_TYPE != "Non-Small Cell Lung Cancer"), ]

# rerun the same tests on subset - takes < 10s for 100k perms
res_nonlung <- compare_group11_to_matched_others(df = perm_data_nonlung, x = setdiff(colnames(perm_data_nonlung), c("SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "GROUP")), seed = 7, B = 100000L)

# View(res$numeric)
# View(res_nonlung$numeric)

# View(res$categorical)
# View(res_nonlung$categorical)

# write to files
write.table(res$numeric, "post_analysis/perm_test_results_g11_numvars.txt", sep = " ", row.names = F)
write.table(res$categorical, "post_analysis/perm_test_results_g11_catvars.txt", sep = " ", row.names = F)
write.table(res_nonlung$numeric, "post_analysis/perm_test_results_g11_nonlung_numvars.txt", sep = " ", row.names = F)
write.table(res_nonlung$categorical, "post_analysis/perm_test_results_g11_nonlung_catvars.txt", sep = " ", row.names = F)


####################################################################################################


