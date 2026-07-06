#
#
#
#--------------------------------------------------------
# Part 1: loading packages and functions
#--------------------------------------------------------
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots, affy, oligoClasses, testit, NanoStringNCTools, NanoStringDiff
)
# function is log2transformed
isLog2Transformed <- function(data) {
    qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    shouldBeLogged <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    return(!shouldBeLogged)
}
#
#
#
#--------------------------------------------------------
# Part 2: Obtaining metadata and raw data
#--------------------------------------------------------
# First obtain metadata
id <- "GSE143382" # GPL27981

# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}

# Obtaining Metadata
metaff <- getGEO(id, GSEMatrix = TRUE, getGPL = FALSE, destdir = ".temp")
length(metaff) #how many datasets

meta <- metaff[[1]]
colnames(exprs(meta))
head(pData(meta))

pData(meta)$"channel_count"

# Metadata wrangling
metadata <- pData(meta) |>
    dplyr::rename(
        condition = "condition:ch1",
        caseID = "title"
    ) |>
    dplyr::select(
        geo_accession,
        supplementary_file,
        condition,
        caseID
    ) |>
    dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1))

#
#
#
#--------------------------------------------------------
# Part 3: Fetching data
#--------------------------------------------------------
## Not following this
options(download.file.method = "curl")
raw <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225759/suppl/GSE225759_RAW.tar"
matrix <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225759/matrix/GSE225759_series_matrix.txt.gz"
download.file(matrix, destfile = "./.temp/GSE225759.gz", method = "curl")
dir.create("./.temp", showWarnings = FALSE, recursive = TRUE)
options(GEOquery.integrate.https = FALSE)
###########################

## Download all the files in the .temp of this environment
for (i in 1:length(metadata$supplementary_file)) {
    url <- metadata$supplementary_file[i]
    destfile <- file.path(paste0(".temp/", id, "/", metadata$file[i]))
    # Download the file
    tryCatch(
        {
            download.file(url, destfile, mode = "wb", method = "curl")
        },
        error = function(e) {
            # Fallback to default method if curl fails
            download.file(url, destfile, mode = "wb")
        }
    )
    # Optional: Extract if it's a tar file
    if (grepl("\\.tar(\\.gz)?$", destfile)) {
        untar(destfile, exdir = paste0(".temp/", id))
    }
}

#####################################################
# testing subtypes, it has supplementary files links
unique(metadata$final_subtype)
datatable(metadata)
#atl_ptcl_geo <- metadata |>
#    filter(final_subtype %in% c("ATLL", "PTCL-TBX21", "PTCL-GATA3", "ENKTL", "AITL", "ALCL-pos", "ALCL-neg", "Intermediate AITL/TBX21", "Intermediate.AITL/GATA3", "Intermediate AITL/GATA3")) |>
#    pull(geo_accession)
#
#
#
#
#--------------------------------------------------------
# Part 3: Reading nanostring data
#--------------------------------------------------------
# Read each file individually and check feature names
feature_lists <- lapply(file.path(".temp", id, metadata$file), function(f) {
  dat <- readNanoStringRccSet(rccFiles = f)
  rownames(dat)
})

# Compare all files against the first file - Validation set has been readed in two different files
mismatches <- sapply(feature_lists, function(feats) !identical(feats, feature_lists[[1]]))
bad_files <- metadata$file[mismatches]
print(bad_files)  # Files that don't match the first file
#
#
#
#--------------------------------------------------------
# Part 3A: Reading dataset
#--------------------------------------------------------
dataset <- readNanoStringRccSet(
  rccFiles = file.path(".temp", id, metadata$file),
  )
dataset
######### Part 3 - merging metadata

pData(dataset) = pData(dataset) |>
  tibble::rownames_to_column("ID") |>
  dplyr::left_join(metadata %>% dplyr::select(-supplementary_file), by=c("ID"="file")) |>
  tibble::column_to_rownames("geo_accession")

colnames(exprs(dataset)) <- rownames(pData(dataset))
colnames(exprs(dataset))
#--------------------------------------------------------
# Part 3B: QC
#--------------------------------------------------------
######### Part 4 - QC
# Basic metrics
qc_metrics <- sData(dataset)
barplot(qc_metrics$FovCount, main = "Fields of View", las = 2)
barplot(qc_metrics$BindingDensity, main = "Binding Density", las = 2)
abline(h = c(0.05, 2.25), col = "red")  # Recommended range <2.25

# Check positive controls (should show linearity)
pos_controls_features <- fData(dataset)[fData(dataset)$CodeClass == "Positive", ]
# Get corresponding expression data
pos_exprs <- exprs(dataset)[rownames(dataset), ]
boxplot(pos_exprs, las=2, log="y", main="Positive Controls")

# Linear plot with expected concentrations
# Extract concentrations from control names (e.g., POS_A(128))
pos_exprs
conc <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", rownames(pos_exprs)))
pos_means <- rowMeans(pos_exprs)

# Plot expected vs observed
plot(log2(conc+1), log2(pos_means), 
     main="Positive Control Linearity",
     xlab="Log2 Expected Concentration", 
     ylab="Log2 Observed Counts",
     pch=16)
abline(lm(log2(pos_means) ~ log2(conc+1)), col="red")

# Calculate R-squared for linearity assessment
model <- lm(log2(pos_means) ~ log2(conc+1))
r_squared <- summary(model)$r.squared
text(min(log2(conc+1)), max(log2(pos_means)), 
     paste("R² =", round(r_squared, 3)),
     pos=4)

# Check negative controls (should show linearity)
neg_controls_features <- fData(dataset)[fData(dataset)$CodeClass == "Negative", ]
neg_exprs <- exprs(dataset)[rownames(neg_controls_features), ]
# Get corresponding expression data
boxplot(neg_exprs+1, las=2, log="y", main="Negative Controls")

# Check housekeeping genes (should be consistent)
housekeeping_features <- fData(dataset)[fData(dataset)$CodeClass == "Housekeeping", ]
house_exprs <- exprs(dataset)[rownames(housekeeping_features), ]
boxplot(house_exprs +1, las = 2, log="y", main = "Housekeeping Genes")
```
#
#--------------------------------------------------------
# Part 3D: Normalization
#--------------------------------------------------------
colnames(exprs(dataset)) <- rownames(pData(dataset))
normalized_matrix <- normalize(dataset, 
"Housekeeping-Log2")

exprs(normalized_matrix) = assayData(normalized_matrix)$exprs_norm

exprs(normalized_matrix) <- log2(exprs(normalized_matrix) + 1)

# Verify normalization
boxplot(exprs(normalized_matrix), main = "After Log2 Transform", las = 2)

# PCA analysis
pca_result <- prcomp(t(exprs(normalized_matrix)))

pData(normalized_matrix)
plot(pca_result$x[,1], pca_result$x[,2], 
     col = as.factor(pData(normalized_matrix)$condition),
     pch = 16, main = "PCA Plot")
legend("topright", legend = levels(as.factor(pData(normalized_matrix)$condition)),
       col = 1:length(levels(as.factor(pData(normalized_matrix)$condition))),
       pch = 16)

biological_genes <- fData(normalized_matrix)[fData(normalized_matrix)$CodeClass == "Endogenous", ]
filtered_matrix <- normalized_matrix[rownames(biological_genes), ]
exprs(filtered_matrix)

#
#
#

#--------------------------------------------------------
# Part 4: Downstream
#--------------------------------------------------------
# Matrix Design
design <- model.matrix(~ 0 + factor(pData(normalized_matrix)$condition))

rownames(design) <- rownames(pData(normalized_matrix))
colnames(design) <- levels(factor(pData(normalized_matrix)$condition))
head(design, 20)
colnames(design) = c("D","eMF","HS","MF")

# Make contrasts (adjust based on your comparisons)
colnames(design)
contrast.matrix <- makeContrasts(
  MF_vs_eMF = MF - eMF,
  MF_vs_D = MF - D,
  eMF_vs_D = eMF - D,
  eMF_vs_HS = eMF - HS,
  MF_vs_HS = MF - HS,
  levels = design
)

fit <- lmFit(filtered_matrix, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

top_probes <- topTable(fit, number = Inf,  adjust.method = "BH", p.value = 0.001)
top_probes

# Add ProbeID as a column and gene names
top_probes_export <- top_probes %>%
  tibble::rownames_to_column("ProbeID") %>%
  left_join(
    fData(normalized_matrix) %>%
      tibble::rownames_to_column("ProbeID") %>%
      select(ProbeID, GeneName),
    by = "ProbeID"
  )
top_probes_export
# Save the enriched version
write.csv(top_probes_export, file = "./nanoString/2707_top_probes_annotated_results.csv", row.names = FALSE)

top_probes |>
  arrange(desc(MF_vs_D)) |>
  dplyr::filter(adj.P.Val < 0.01) |>
  head(50)

#
#
#
#
#
#
#
#
#
