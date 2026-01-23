### PART A: Load the data

setwd("/Users/pujasanjaythorat/Downloads")
install.packages("readxl")
library(readxl)

# then:
gene_expression <- read_excel("Gene_Expression_Data.xlsx")
gene_info       <- read.csv("Gene_Information.csv", stringsAsFactors=FALSE)
sample_info     <- read.table("Sample_Information.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

### PART B: Change sample names

# Make a copy of sample_info for safety
samples <- sample_info

# Clean patient column: extract just the number
samples$patient <- gsub("patient: ", "", samples$patient)

# Add unique sample names combining phenotype + patient number
samples$new_name <- paste0(samples$group, "_", samples$patient)

# Extract original GSM sample names from gene_expression
old_names <- colnames(gene_expression)[-1]   # first column is Probe_ID

# Make sure the order matches: sample_info must be in same order as the columns
if(!all(old_names == rownames(samples))) {
  # If rownames not set, set them
  rownames(samples) <- old_names
}

# Create the vector of new names
new_names <- c("Probe_ID", samples$new_name)

# Rename columns
colnames(gene_expression) <- new_names

# Check result
colnames(gene_expression)


### PART C: merge gene_expression with gene_info using Probe_ID

merged_data <- merge(gene_expression, gene_info, by = "Probe_ID", all.x = TRUE)

# Check the merged result
str(merged_data)

# Identify tumor and normal columns based on their names
tumor_cols  <- grep("^tumor_", colnames(merged_data),  value = TRUE)
normal_cols <- grep("^normal_", colnames(merged_data), value = TRUE)

# Keep Probe_ID as first column
probe_col <- "Probe_ID"

# Create tumor-only dataset
tumor_data <- merged_data[, c(probe_col, tumor_cols)]

# Create normal-only dataset
normal_data <- merged_data[, c(probe_col, normal_cols)]

# Inspect results
str(tumor_data)
str(normal_data)


### PART D: Compute average expression per gene for tumor and normal datasets

# --- QC1: Check tumor and normal datasets exist ---
if (!exists("tumor_data")) stop("ERROR: tumor_data missing.")
if (!exists("normal_data")) stop("ERROR: normal_data missing.")

# --- QC2: Check Probe_ID column exists in both ---
if (!"Probe_ID" %in% colnames(tumor_data)) stop("ERROR: Probe_ID missing in tumor_data.")
if (!"Probe_ID" %in% colnames(normal_data)) stop("ERROR: Probe_ID missing in normal_data.")

# --- Compute average tumor expression ---
avg_tumor <- rowMeans(tumor_data[ , -1], na.rm = TRUE)

# --- Compute average normal expression ---
avg_normal <- rowMeans(normal_data[ , -1], na.rm = TRUE)

# --- Assign probe names to vectors ---
names(avg_tumor) <- tumor_data$Probe_ID
names(avg_normal) <- normal_data$Probe_ID

# QC3: Ensure gene/probe ordering matches
if (!all(names(avg_tumor) == names(avg_normal))) {
  stop("ERROR: Probe_ID mismatch between tumor and normal datasets.")
}

# Check a few values
head(avg_tumor)
head(avg_normal)


### PART E: Determine the log2 fold change for each Probe 

# --- QC1: Ensure average expression vectors exist ---
if (!exists("avg_tumor")) stop("ERROR: avg_tumor does not exist. Run Part 1(d).")
if (!exists("avg_normal")) stop("ERROR: avg_normal does not exist. Run Part 1(d).")

# --- QC2: Ensure both vectors have identical lengths ---
if (length(avg_tumor) != length(avg_normal)) {
  stop("ERROR: Tumor and Normal vectors are not the same length.")
}

# --- QC3: Ensure gene names match exactly ---
if (!all(names(avg_tumor) == names(avg_normal))) {
  stop("ERROR: Gene names do not match between avg_tumor and avg_normal.")
}

# --- Compute fold change using YOUR EXACT FORMULA ---
# FC = (Tumor – Control) / Control
fc_raw <- (avg_tumor - avg_normal) / avg_normal

# --- QC4: Check for division-by-zero problems ---
if (any(avg_normal == 0)) {
  warning("WARNING: Some Control (Normal) values are zero ??? FC may produce Inf or NaN.")
}

# Now compute log2 fold change
log2_fc <- log2(fc_raw)

# --- QC5: Detect problematic values (NaN, Inf) ---
if (any(is.na(log2_fc))) {
  warning("WARNING: NaN values created during log2FC calculation.")
}
if (any(is.infinite(log2_fc))) {
  warning("WARNING: Infinite values detected ??? division by zero possible.")
}

# --- Output tidy results ---
fold_change_data <- data.frame(
  Probe_ID   = names(log2_fc),
  log2_FC    = log2_fc,
  Tumor_Avg  = avg_tumor,
  Normal_Avg = avg_normal,
  stringsAsFactors = FALSE
)

# Show the first few rows
head(fold_change_data)

### PART F: Identify genes with |log2FC| > 5 

# --- Make avg_data and normalize column names expected later ----
avg_data <- fold_change_data

# Rename log2_FC -> log2FC if present
if ("log2_FC" %in% colnames(avg_data) && !("log2FC" %in% colnames(avg_data))) {
  names(avg_data)[names(avg_data) == "log2_FC"] <- "log2FC"
}

# Also create lowercase tumor_mean / normal_mean aliases if later code expects them
# (This keeps your original column names and only adds aliases)
if ("Tumor_Avg" %in% colnames(avg_data) && !("tumor_mean" %in% colnames(avg_data))) {
  avg_data$tumor_mean <- avg_data$Tumor_Avg
}
if ("Normal_Avg" %in% colnames(avg_data) && !("normal_mean" %in% colnames(avg_data))) {
  avg_data$normal_mean <- avg_data$Normal_Avg
}


# Step 1: QC — Ensure required columns exist
stopifnot("Probe_ID" %in% colnames(avg_data))
stopifnot("log2FC" %in% colnames(avg_data))
stopifnot("Probe_ID" %in% colnames(gene_info))

# Step 2: Calculate absolute value of log2FC
avg_data$abs_log2FC <- abs(avg_data$log2FC)

# Step 3: Filter genes with magnitude > 5
fc_filtered <- subset(avg_data, abs_log2FC > 5)

# QC: Check if any genes were returned
if(nrow(fc_filtered) == 0) {
  message("No genes found with |log2FC| > 5.")
}

# Step 4: Merge with gene annotation data
fc_annotated <- merge(fc_filtered, gene_info, by = "Probe_ID", all.x = TRUE)

# Step 5: QC — Check for NA gene annotations
fc_annotated$Missing_Annotation <- is.na(fc_annotated$Symbol)

# Step 6: Inspect final filtered + annotated table
head(fc_annotated)
summary(fc_annotated$abs_log2FC)
table(fc_annotated$Missing_Annotation)

### PART G: Add column showing whether gene is higher in Tumor or Normal

# QC 1: Ensure needed columns exist
required_cols <- c("Probe_ID", "tumor_mean", "normal_mean", "log2FC")
missing_cols <- setdiff(required_cols, colnames(fc_annotated))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# QC 2: Ensure numeric columns are indeed numeric
stopifnot(is.numeric(fc_annotated$tumor_mean))
stopifnot(is.numeric(fc_annotated$normal_mean))
stopifnot(is.numeric(fc_annotated$log2FC))

# Step: Add expression direction column
fc_annotated$Higher_Expression <- ifelse(
  fc_annotated$tumor_mean > fc_annotated$normal_mean,
  "Tumor",
  "Normal"
)

# QC 3: Verify column creation
stopifnot("Higher_Expression" %in% colnames(fc_annotated))

# QC 4: Table summary of how many genes higher in each phenotype
table(fc_annotated$Higher_Expression)

# View first part of the final dataset
head(fc_annotated)

### PART 2: A EDA ON fc_annotated
### ANALYSIS ONLY (NO GRAPHS)

### ---------- QC 1: Check object exists ----------
if (!exists("fc_annotated")) {
  stop("ERROR: 'fc_annotated' dataset not found. Please run Part 1f and 1g first.")
}

### ---------- QC 2: Check required columns ----------
required_cols <- c("Probe_ID", "Symbol", "tumor_mean", "normal_mean", "log2FC", "Higher_Expression")
missing_cols <- setdiff(required_cols, colnames(fc_annotated))

if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing required columns:",
             paste(missing_cols, collapse = ", ")))
}

### ---------- QC 3: Check for NA values ----------
cat("\n--- QC: Missing Values Check ---\n")
print(colSums(is.na(fc_annotated)))


### ---------- EDA SUMMARY STATISTICS ----------
cat("\n--- SUMMARY: Tumor Mean Expression ---\n")
print(summary(fc_annotated$tumor_mean))

cat("\n--- SUMMARY: Normal Mean Expression ---\n")
print(summary(fc_annotated$normal_mean))

cat("\n--- SUMMARY: log2 Fold Change ---\n")
print(summary(fc_annotated$log2FC))

cat("\n--- SUMMARY: Up vs Down Regulated Gene Counts ---\n")
print(table(fc_annotated$Higher_Expression))


### ---------- Identify Extreme Fold Change Genes ----------
threshold <- 5
extreme_genes <- subset(fc_annotated, abs(log2FC) > threshold)

cat("\n--- EXTREME FC GENES (|log2FC| > 5) ---\n")
cat("Number of genes:", nrow(extreme_genes), "\n\n")

cat("--- Top 10 Upregulated Genes ---\n")
print(head(extreme_genes[order(-extreme_genes$log2FC), ], 10))

cat("\n--- Top 10 Downregulated Genes ---\n")
print(head(extreme_genes[order(extreme_genes$log2FC), ], 10))


### Part 2 B: Histogram of DEGs by Chromosome

### QC 1: ensure dataset exists
if (!exists("fc_annotated")) {
  stop("ERROR: fc_annotated not found. Run Part 1f and 1g first.")
}

### QC 2: ensure chromosome column exists
if (!"Chromosome" %in% colnames(fc_annotated)) {
  stop("ERROR: Chromosome column missing from fc_annotated.")
}

### Clean and extract numeric chromosome

# Remove empty or NA chromosome values
fc_clean <- subset(fc_annotated, Chromosome != "" & !is.na(Chromosome))

# Extract numeric chromosome from "3p21", "12q", etc.
fc_clean$Chromosome_Num <- as.numeric(sub("^([0-9]+).*", "\\1", fc_clean$Chromosome))

### QC 3: Check extracted values
cat("Extracted chromosome numbers:\n")
print(sort(unique(fc_clean$Chromosome_Num)))

### Define differentially expressed genes (DEGs)
### log2FC threshold = 1

deg_data <- subset(fc_clean, abs(log2FC) > 1)

### QC 4: Check DEG count
cat("Number of DEGs identified:", nrow(deg_data), "\n")

### HISTOGRAM (with chromosome labels fixed)

# Get unique chromosome numbers
chroms <- sort(unique(deg_data$Chromosome_Num))

# Make histogram with controlled breaks
hist(deg_data$Chromosome_Num,
     breaks = seq(min(chroms)-0.5, max(chroms)+0.5, by = 1),
     xaxt = "n",                     # disable default x-axis
     main = "Distribution of Differentially Expressed Genes by Chromosome",
     xlab = "Chromosome",
     col = "pink")

# Add chromosome labels manually
axis(1, at = chroms, labels = chroms)

### Part 2 C: DEG Histogram by Chromosome & Sample Type

### QC 1: dataset exists
if (!exists("fc_annotated")) {
  stop("ERROR: 'fc_annotated' not found.")
}

### QC 2: required columns exist
required_cols <- c("Chromosome", "log2FC", "Higher_Expression")
missing_cols <- setdiff(required_cols, colnames(fc_annotated))

if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing required columns:", 
             paste(missing_cols, collapse=", ")))
}

### Remove missing chromosome entries
fc_clean <- subset(fc_annotated, Chromosome != "" & !is.na(Chromosome))

### Extract numeric chromosome portion
fc_clean$Chromosome_Num <- as.numeric(sub("^([0-9]+).*", "\\1", fc_clean$Chromosome))

### QC 3: check chromosome extraction
if (all(is.na(fc_clean$Chromosome_Num))) {
  stop("ERROR: Chromosome number extraction failed — check input format.")
}

### Define DEGs
deg_data <- subset(fc_clean, abs(log2FC) > 1)

### QC 4: check DEGs found
if (nrow(deg_data) == 0) {
  stop("ERROR: No DEGs found at log2FC > 1. Cannot plot histograms.")
}

### Split DEGs by group
deg_normal <- subset(deg_data, Higher_Expression == "Normal")
deg_tumor  <- subset(deg_data, Higher_Expression == "Tumor")

### Helper function: generate fixed breaks for histogram
chrom_breaks <- function(x) {
  ux <- sort(unique(x))
  seq(min(ux)-0.5, max(ux)+0.5, by = 1)
}

### Combined Histogram: Tumor vs Normal (Side-by-Side)

if (nrow(deg_tumor) == 0 & nrow(deg_normal) == 0) {
  warning("No DEGs found — cannot plot combined histogram.")
} else {
  
  # Tabulate counts for each chromosome
  tab_normal <- table(deg_normal$Chromosome_Num)
  tab_tumor  <- table(deg_tumor$Chromosome_Num)
  
  # Ensure same chromosome names
  all_chroms <- sort(unique(c(names(tab_normal), names(tab_tumor))))
  
  counts_normal <- as.numeric(tab_normal[all_chroms])
  counts_tumor  <- as.numeric(tab_tumor[all_chroms])
  
  # Replace NAs with zeros
  counts_normal[is.na(counts_normal)] <- 0
  counts_tumor[is.na(counts_tumor)] <- 0
  
  # Create matrix
  count_matrix <- rbind(Tumor = counts_tumor,
                        Normal = counts_normal)
  
  # Side-by-side barplot
  barplot(count_matrix,
          beside = TRUE,
          col = c("lightcoral", "lightblue"),
          names.arg = all_chroms,
          main = "Distribution of DEGs by Chromosome (Normal vs Tumor)",
          xlab = "Chromosome",
          ylab = "Frequency")
  
  legend("topright",
         fill = c("lightcoral", "lightblue"),
         legend = c("Tumor DEGs", "Normal DEGs"))
}

### Part 2 D: Bar Chart of DEG Percentages 

### QC 1 — Check dataset exists
if (!exists("fc_annotated")) {
  stop("ERROR: 'fc_annotated' not found. Run earlier steps first.")
}

### QC 2 — Check required columns
required_cols <- c("log2FC", "Higher_Expression")
missing_cols <- setdiff(required_cols, colnames(fc_annotated))

if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing required columns:", 
             paste(missing_cols, collapse = ", ")))
}

### Define DEGs (|log2FC| > 1)
deg_data <- subset(fc_annotated, abs(log2FC) > 1)

### QC 3 — Ensure DEGs exist
if (nrow(deg_data) == 0) {
  stop("ERROR: No DEGs found at |log2FC| > 1 — cannot compute percentages.")
}

### Count DEGs that are higher in Tumor vs Normal
deg_counts <- table(deg_data$Higher_Expression)

### QC 4 — Ensure both groups exist
if (!("Tumor" %in% names(deg_counts))) {
  deg_counts["Tumor"] <- 0
}
if (!("Normal" %in% names(deg_counts))) {
  deg_counts["Normal"] <- 0
}

### Compute percentages
total_degs <- sum(deg_counts)
percent_up   <- (deg_counts["Tumor"]  / total_degs) * 100
percent_down <- (deg_counts["Normal"] / total_degs) * 100

percentages <- c(percent_up, percent_down)
names(percentages) <- c("Up in Tumor", "Down in Tumor")

### BAR CHART
barplot(
  percentages,
  main = "Percentages of DEGs Up/Down in Tumor Samples",
  ylab = "Percentage (%)",
  col = c("steelblue", "tomato"),
  ylim = c(0, 100)
)

### Optional: Legend
legend("topright",
       legend = c("Upregulated in Tumor", "Downregulated in Tumor"),
       fill = c("steelblue", "tomato"))


### Part 2 E: – Raw data Heatmap to visualize gene expression by sample

# QC1: Ensure dataset exists
if (!exists("gene_expression")) {
  stop("ERROR: 'gene_expression' not found. Load it first.")
}

# QC2: Detect probe/gene ID column
probe_col <- names(gene_expression)[sapply(gene_expression, function(x) !is.numeric(x))][1]

if (is.na(probe_col)) {
  stop("ERROR: No non-numeric probe/gene column found in dataset.")
}

cat("Detected probe column:", probe_col, "\n")

# Remove probe column and store IDs
probe_ids <- gene_expression[[probe_col]]
expr_matrix <- gene_expression[, colnames(gene_expression) != probe_col]

# QC3: Ensure numeric sample columns
if (any(sapply(expr_matrix, function(x) !is.numeric(x)))) {
  stop("ERROR: Some sample columns contain non-numeric values.")
}

# Convert to matrix
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- probe_ids

# Check for missing values
if (any(is.na(expr_matrix))) {
  warning("WARNING: NA values detected in expression matrix.")
}

# Scale rows
expr_scaled <- scale(expr_matrix)

# High-resolution PNG output
png("heatmap_highres.png", width = 2500, height = 2000, res = 300)

# Custom color palette: white ??? yellow ??? light red ??? red
my_cols <- colorRampPalette(c("white", "yellow", "orange", "red"))(256)

heatmap(expr_scaled,
        main = "High-Resolution Heatmap of Raw Gene Expression",
        xlab = "Samples",
        ylab = "Genes (Probes)",
        col = my_cols,
        margins = c(8, 8))

dev.off()

### Part 2 F:– Clustermap of Raw Gene Expression

# QC1: Ensure dataset exists
if (!exists("gene_expression")) {
  stop("ERROR: 'gene_expression' dataset not found. Load it first.")
}

# QC2: Detect probe column (first non-numeric)
probe_col <- names(gene_expression)[sapply(gene_expression, function(x) !is.numeric(x))][1]

if (is.na(probe_col)) {
  stop("ERROR: Could not detect probe column. Make sure gene IDs exist.")
}

cat("Detected Probe Column:", probe_col, "\n")

# Store probe IDs and keep numeric matrix
probe_ids <- gene_expression[[probe_col]]
expr_matrix <- gene_expression[, colnames(gene_expression) != probe_col]

# QC3: Ensure remaining columns are numeric
if (any(sapply(expr_matrix, function(x) !is.numeric(x)))) {
  stop("ERROR: Some sample columns contain non-numeric values.")
}

# Convert to matrix
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- probe_ids

# QC4: Missing values check
if (any(is.na(expr_matrix))) {
  warning("WARNING: NA values detected in expression matrix.")
}

# Scale rows for better visualization
expr_scaled <- scale(expr_matrix)

# --- Create hierarchical clustering ---
row_clust <- hclust(dist(expr_scaled))    # genes
col_clust <- hclust(dist(t(expr_scaled))) # samples

# --- Generate clustermap ---
heatmap(
  expr_scaled,
  Rowv = as.dendrogram(row_clust),
  Colv = as.dendrogram(col_clust),
  main = "Clustered Heatmap of Gene Expression",
  xlab = "Samples",
  ylab = "Genes (Probes)",
  col = heat.colors(256)
)


### PART G: Explaining the findings of analysis

# The differential expression analysis revealed a clear separation between tumor and normal samples, indicating strong underlying biological differences.
# A substantial number of genes showed large fold-change values (|log2FC| > 5), suggesting significant transcriptional dysregulation in the tumor group.
# The majority of DEGs were found to be upregulated in tumor samples, indicating increased transcriptional activity associated with tumorigenic processes.
# Chromosome-level summaries showed that certain chromosomes contained higher densities of DEGs, suggesting potential chromosome-specific abnormalities or regulatory hotspots.
# Histograms separating DEGs by phenotype showed that tumor samples contributed more to DEG counts, further supporting the observed gene dysregulation pattern.
# Heatmap and clustermap visualizations demonstrated strong clustering by phenotype, with tumor and normal samples forming distinct groups based on gene expression patterns.
# Gene-level clustering also revealed groups of genes with shared expression behavior, indicating coordinated regulatory programs that differ between tumor and normal tissues.


