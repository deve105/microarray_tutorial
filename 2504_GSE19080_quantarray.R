#--------------------------------------------------------
# Title: Microarray Data Analysis with GEOquery and limma
"""
#- Gene expression profiling in patients infected with HTLV-1: Identification of ATL and HAM/TSP-specific genetic profiles

Customized CHIP
gene expression profiles of CD4+ T-cells isolated from 7 ATL, 12 HAM/TSP and 11 AC (asymptomatic carriers).
Agilent two colors, but HD were read in a different protocol.

"""
#--------------------------------------------------------
# Part 1: loading packages and functions
#--------------------------------------------------------
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots, affy, oligoClasses, testit
)

# function is log2transformed
isLog2Transformed <- function(data) {
    qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    shouldBeLogged <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    return(!shouldBeLogged)
}

#--------------------------------------------------------
# Part 2: Obtaining metadata and raw data
#--------------------------------------------------------

# First obtain metadata
id <- "GSE19080" # GPL4133

# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}

# Obtaining Metadata
meta <- getGEO(id, GSEMatrix = TRUE, destdir = ".temp")
meta <- meta[[1]]

#Basic exploration of the metadata
head(fData(meta), 2)
head(exprs(meta), 2)
pData(meta)
dim(fData(meta))
colnames(exprs(meta))
dim(meta)
head(pData(meta), 5)
head(fData(meta))

annotx = fData(meta) #>
annotx |>
    filter(GENE_NAME=="")

# Metadata wrangling
pd <- pData(meta) |>
    tibble::rownames_to_column("ID") |>
    dplyr::select("characteristics_ch1", "characteristics_ch1.3", ID, supplementary_file) |>
    dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
    dplyr::rename(atl_subtype = "characteristics_ch1", gender = "characteristics_ch1.3") |>
    dplyr::mutate(across(
        c(atl_subtype, gender),
        ~ stringr::str_remove(., "^.*: ")
    )) |>
    mutate(atl_subtype2 = case_when(
        str_detect(atl_subtype, "ATL") ~ "ATL",
        str_detect(atl_subtype, "HAM") ~ "HAMTSP",
        str_detect(atl_subtype, "AC") ~ "AC",
        str_detect(atl_subtype, "Healthy") ~ "HD",
        .default = NA
    ))

pd7 = pd 

## Download all the files in the .temp of this environment
for (i in 1:length(pd7$supplementary_file)) {
    url <- pd7$supplementary_file[i]
    destfile <- file.path(paste0(".temp/", id, "/", pd7$file[i]))
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

#--------------------------------------------------------
# Part 3: Capturing targets
#--------------------------------------------------------
# Exploring the structure of one file
con <- gzfile(file.path(".temp/GSE19080/GSM472356_HISH0317.txt.gz"))
file_lines <- readLines(con, n=50)
close(con)
file_lines

# Subsettting files (here we have two different platforms) 
pd7 = pd |>
    filter(ID%in%c(paste0("GSM4723", c(56:73, 82:93))))

# Re-read the raw files, this was a quantarray (described in metadata)
agilent_data <- read.maimages(
    files = file.path(".temp", id, pd7$file),
    source = "quantarray",
    green.only = FALSE,
    names = pd7$ID,
    other.columns = list(Flag = "Ignore Filter"))


#--------------------------------------------------------
# Part 4: Manipulating targets
#--------------------------------------------------------
# Convert targets to tibble for easy joining
agilent_data$targets = agilent_data$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")


#### gene Annotation
gpl2 <- getGEO("GPL9686")

annot <- Table(gpl2)[, c("SYMBOL", "GB_ACC", "ID")]
annot$dataset <- "GPL9686"
"""
BLANK: Empty spots with no DNA
Arabidopsis: Plant DNA (shouldn't hybridize to human samples)
SSC: Salmon Sperm DNA (unrelated to human targets)
"""
controls=c("Arabidopsis", "SSC", "BLANK")


head(agilent_data$genes)
# Add gene symbols to agilent_data

agilent_data$genes <- agilent_data$genes |>
    dplyr::select(Name) |>
    dplyr::left_join(annot, by = c("Name" = "GB_ACC")) |>
    dplyr::mutate(is_control = factor(
        case_when(
            Name %in% controls & is.na(dataset) ~ "control",
            !(Name %in% controls) & is.na(dataset) ~ "UNK",
            dataset == "GPL9686" ~ "gene"
        )
    ))

agilent_data$
table(agilent_data$genes$is_control, useNA="always")

control_status
#--------------------------------------------------------
# Part 5: QC
#--------------------------------------------------------
# (A) Use negative controls for background correction
agilent_two_color_qc <- function(agilent_data,
                                verbose = TRUE) {
    # 1. Input Validation --------------------------------------------------------
    required_components <- c("R", "G", "genes", "targets")
    missing_comps <- setdiff(required_components, names(agilent_data))
    if(length(missing_comps)) {
        stop("Missing required components: ", paste(missing_comps, collapse=", "))
    }
    # 2. Dimensional Checks ----------------------------------------------------
    if(verbose) message("\nVerifying dimensions...")
    dim_checks <- list(
        list(nrow(agilent_data$R), nrow(agilent_data$G), "Channel row mismatch"),
        list(ncol(agilent_data$R), ncol(agilent_data$G), "Channel column mismatch"),
        list(nrow(agilent_data$R), nrow(agilent_data$genes), "Gene annotation row mismatch"),
        list(ncol(agilent_data$R), nrow(agilent_data$targets), "Sample annotation mismatch")
    )
    
    for(chk in dim_checks) {
        if(chk[[1]] != chk[[2]]) stop(chk[[3]])
    }

    # 3. Missing Value Analysis ------------------------------------------------
    if(verbose) message("\nAnalyzing missing values...")
    na_stats <- list(
        R = list(
            probes = rowMeans(is.na(agilent_data$R)),
            samples = colMeans(is.na(agilent_data$R))
        ),
        G = list(
            probes = rowMeans(is.na(agilent_data$G)),
            samples = colMeans(is.na(agilent_data$G))
        )
    )
    
    # 4. Background correction
    control_status <- as.numeric(agilent_data$genes$is_control == "control")
    agilent_data <- limma::backgroundCorrect(agilent_data, method="normexp", offset=20) 
    # 5. Normalization
    agilent_data <- limma::normalizeWithinArrays(agilent_data, method="loess")
    agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")

    # 6. technical replicates reduction by probename
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$Name)
    agilent_data <- agilent_data[agilent_data$genes$is_control=="gene", ]
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$SYMBOL)
    agilent_data$genes <- droplevels(agilent_data$genes)
    return(agilent_data)
}

## Apply the function to the data
gse1 = agilent_two_color_qc(agilent_data)

##### PCA
dim(gse1)

boxplot(gse1$M, main="Normalized Negative Controls")
pca <- prcomp(t(gse1$M))
pca = pca$x |>
    as.data.frame() 
pca |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(pd7, by="ID") |> 
    tidyplots::tidyplot(x=PC1, y=PC2, color=atl_subtype) |>
    tidyplots::add_data_points()
#--------------------------------------------------------
# Part 6: Read the additiional data *Healthy donors*
#--------------------------------------------------------
con <- gzfile(file.path(".temp/GSE19080/GSM472379_HISH0562.txt.gz"))
file_lines <- readLines(con, n = 50)
close(con)
file_lines
?
pd7_2 <- pd |>
    filter(ID %in% c(paste0("GSM4723", 74:81)))

agilent_data2 <- read.maimages(
    files = file.path(".temp", id, pd7_2$file),
    source = "genepix.custom",
    green.only = FALSE,
    names = pd7_2$ID,
       other.columns = list(
        Flags = "Autoflag",
        ControlType = "ControlType",   # Specific to GenePix
        Block = "Block"
))

agilent_data2
head(agilent_data$printer)
table(agilent_data2$other$Autoflag, useNA="always")

# For GenePix data (you're already using this)
agilent_data2$genes$is_control <- grepl("^(CONTROL|CTL|NC_|PC_)", agilent_data2$genes$Annot)
head(agilent_data2$genes)


# Convert targets to tibble for easy joining
agilent_data2$targets <- agilent_data2$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7_2, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")

#### gene Annotation
gpl2 <- getGEO("GPL9686")
annot <- Table(gpl2)[, c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")]
head(Table(gpl2))
# Add gene symbols to agilent_data
head(agilent_data2$genes)
agilent_data2$genes <- agilent_data2$genes |>
    tidyr::separate(
        Name,
        into = c("GENE_SYMBOL", "Annot", "Name"),
        sep = ":"
    ) |>
    dplyr::left_join(annot, by = "ID") #|>
    mutate(is_control = factor(ifelse(is.na(SYMBOL) & is.na(GENE_NAME), "control", "gene")))

agilent_data2$genes$is_control <- grepl("^(CONTROL|CTL|NC_|PC_)", agilent_data2$genes$Annot)


table(agilent_data2$genes$is_control)    

gse2 = agilent_two_color_qc(agilent_data2)
dim(gse2)
dim(gse1)
gse1$targets$dataset= "dataset1"
gse2$targets$dataset= "dataset2"

#--------------------------------------------------------
# Part 7: Joining both datasets
#--------------------------------------------------------
common_genes = intersect(gse1$genes$Name, gse2$genes$Name)

# Extract normalized M- and A-values for both datasets
M1 <- gse1$M[gse1$genes$Name %in% common_genes, ]
M2 <- gse2$M[gse2$genes$Name %in% common_genes, ]
# For A values (absolute expression)
A1 <- gse1$A[match(common_genes, gse1$genes$Name), ]
A2 <- gse2$A[match(common_genes, gse2$genes$Name), ]

# Add row names
rownames(M1) <- rownames(A1) <- common_genes
rownames(M2) <- rownames(A2) <- common_genes

# Combine M values (log ratios)
combined_M <- cbind(M1, M2)

# Combine A values (average intensities)
combined_A <- cbind(A1, A2)

# Batch correct before combining
pacman::p_load(sva)
# Create batch information
batch <- c(rep(1, ncol(M1)), rep(2, ncol(M2)))

# Create model matrix with intercept for ComBat
mod <- model.matrix(~1, data = data.frame(Intercept = rep(1, ncol(combined_M))))

combined_M_batch_corrected <- ComBat(
    dat = combined_M,
    batch = batch,
    mod = mod
)

# Same for A values if needed
combined_A_batch_corrected <- ComBat(
    dat = combined_A,
    batch = batch,
    mod = mod
)

gse1$targets = gse1$targets |>
    dplyr::select(atl_subtype2)

gse2$targets <- gse2$targets |>
    dplyr::select(atl_subtype2)
# Combine phenotype data
combined_targets <- rbind(
    gse1$targets,
    gse2$targets
)
# Verify disease status column is consistent
table(combined_targets$atl_subtype2)

gse_combined <- list(
  M = combined_M_batch_corrected,
  A = combined_A_batch_corrected,
  targets = combined_targets,
  genes = data.frame(
    SYMBOL = common_genes,
    stringsAsFactors = FALSE
  ))

class(gse_combined) <- c("MAList", "list")

# Now you can use for differential expression
design_combined <- model.matrix(~ 0 + factor(combined_targets$atl_subtype2)) |>
    `colnames<-`(levels(factor(combined_targets$atl_subtype2))) |>
    `rownames<-`(rownames(combined_targets))

contrast.matrix <- makeContrasts(
    ATL_AC = ATL - AC,
    HAM_AC = HAMTSP - AC,
    ATL_HD = ATL - HD,
    HAM_HD = HAMTSP - HD,
    AC_HD = AC - HD,
    levels = design_combined
)

fit2 <- lmFit(gse_combined, design_combined) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

top_probes <- topTable(fit2, number = Inf, adjust.method = "BH" ) |>#,  p.value = 0.01) |>
    as.data.frame() |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(annot[, c("GB_ACC", "SYMBOL")], by = c("ID" = "GB_ACC")) |>
    dplyr::select(-ID) |>
    dplyr::arrange(desc(adj.P.Val))

top_probes$new_SYMBOL <- update_gene_symbols(top_probes$SYMBOL.y)
head(top_probes)

top_probes |>
filter(SYMBOL.y %in% htlv_apc | new_SYMBOL %in% htlv_apc) #|>
    arrange(desc(HAM_AC)) |>
    head(30)
#--------------------------------------------------------
# Part 7: DE
#--------------------------------------------------------
# Create design matrix based on final_subtype
pd7
design <- pd7 %>%
    mutate(
        atl_factor = factor(atl_subtype2)
    ) %>%
    modelr::model_matrix(~ 0 + atl_factor) %>%
    as.data.frame() %>%
    `rownames<-`(pd7$ID) |>
    `colnames<-`(levels(factor(pd7$atl_subtype2)))


# Make contrasts (adjust based on your comparisons)
array_weights <- arrayWeights(gse1)

levels(factor(pd7$atl_subtype2))
contrast.matrix <- makeContrasts(
    ATL_AC = ATL - AC,
    HAM_AC = HAMTSP - AC,
    ATL_HAM = ATL - HAMTSP,
    levels = design
)
fit2 <- lmFit(gse1, design, weights=array_weights) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf, adjust.method = "BH") #,  p.value = 0.01)
head(top_probes)


# Load required library
library(org.Hs.eg.db)

# Function to update gene symbols
update_gene_symbols <- function(gene_symbols) {
    # Map old gene symbols to new ones using org.Hs.eg.db
    updated_symbols <- mapIds(
        org.Hs.eg.db,
        keys = gene_symbols,
        column = "SYMBOL",
        keytype = "ALIAS",
        multiVals = function(x) paste(unique(x), collapse = ";")
    )
    # Replace NA values with original symbols if no mapping is found
    updated_symbols[is.na(updated_symbols)] <- gene_symbols[is.na(updated_symbols)]
    return(factor(updated_symbols))
}
head(top_probes)
# Example: Update gene symbols in the annotation data
top_probes$new_SYMBOL <- update_gene_symbols(top_probes$SYMBOL)
top_probes$new_SYMBOL
# Verify updated gene symbols
setdiff(table$new_SYMBOL, table$SYMBOL)

head(table)


#--------------------------------------------------------
# Part 7: Manipulation
#--------------------------------------------------------
"""
Proximal TCR signaling: CTLA-4, PD-1, SHP-1, SHP-2, LYP, Cbl-b, GRAIL, SIT, PAG, Dok family members, Drak2, and CD5. 
These regulators often target kinases like Lck and ZAP-70, adaptor proteins like LAT and SLP-76, or the CD3ζ chain itself.   
Calcium signaling: Some negative regulators, such as CTLA-4 and Dok family members, can influence calcium mobilization and downstream signaling events.   
MAPK pathway: Negative regulation of the MAPK pathway, which is critical for gene expression and effector functions, is mediated by CTLA-4, PD-1, SHP-2, Dok family members, and RASA2.   
NF-κB pathway: Several negative regulators, including CTLA-4, PD-1, Cbl-b, Itch, A20, and Peli1, target the NF-κB signaling pathway, which controls the expression of genes involved in inflammation, apoptosis, and immune responses.   
NFAT pathway: CTLA-4, PD-1, and MDM2 are among the negative regulators that can modulate the activity of NFAT transcription factors, which are critical for cytokine production.   
PI3K/Akt/mTOR pathway: The PI3K/Akt/mTOR pathway, involved in cell survival, growth, and metabolism, is targeted by negative regulators such as PD-1, Cbl-b, Peli1, and TSC1/TSC2.
"""
early_tcr=c("LCK", "CD3D", "CD3E", "CD3G", "LAT", "LCP2") # ZAP70, SLP76, ITK

neg_tcr=c("CD5", "PTPN6", "PTPN22","SOCS1", "CBLB", "DUSP14", "CD6", "PDCD1LG2", "PTPN12", "DOK1", "MDM2", "JUNB", "SHC1", "UBD", "ATM", "BTK", "IL1RL1", "IFNGR2", "IFNGR1", "GNAI3", "SOCS2") #DRAK2 , PTPN11

htlv_tcr <- c("MYD88", "CD247", "TIGIT", "CD101", "IL2RA", "CD58", "ITK", "TXK", "GATA3", "IRF4", "CBLB", "PTPN12", "ATM", "BTK", "SOCS2")

htlv_apc <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "ERAP1", "NLRC5", "PSMB1", "PSMB2", "PSMB8", "PSMB9", "PSMB10", "TAP1", "TAP2", "TAPBP")


il17df = gse1$A["T75556", ] |>
    as.data.frame() |>
    tibble::rownames_to_column("IL17RA") |>
    pivot_longer(-IL17RA, names_to="ID", values_to="IL17RAvalue") |>
    dplyr::select(-ID) 

top_probes |>
filter(SYMBOL %in% htlv_apc | new_SYMBOL %in% htlv_apc) #|>
    pull(Name)

annotation_df <- data.frame(
    row.names = pd7$ID,
    Subtype = pd7$atl_subtype2
)

tcrneg_heatmapgenes <- gse1$A[gse1$genes$SYMBOL %in% htlv_apc, ] |>
    as.data.frame() |>
    tibble::rownames_to_column("GB_ACC") |>
    left_join(annot[, c("GB_ACC", "SYMBOL")], by = c("GB_ACC")) |>
    tibble::column_to_rownames("SYMBOL") |>
    dplyr::select(-GB_ACC) |>
    as.matrix()

pheatmap(
    mat = heatmapgenes,
    scale = "row",
    annotation_col = annotationdf,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellheight = 15,
    main = "TCR Negative Regulators",
    filename = "2504_negativeTCRregulators.png",
    width = 10, # Nature standard single-column
    height = 7,
    units = "in",
    family = "Arial" # Embed font
)


head(gse1)
top_probes |>  
  dplyr::filter(SYMBOL =="Il17ra") |>
  tibble::rownames_to_column("ID") |>
  pull(ID)
    #dplyr::filter(SYMBOL %in% early_tcr) |>
    #as_tibble() 

top_probes |>
    dplyr::select(-c(1:5, 8)) |>
    #dplyr::filter(GENE_SYMBOL =="ZNF856B") #|>
    dplyr::filter(SYMBOL %in% neg_tcr) |>
    as_tibble() #|>
    #arrange(desc(ATL_AC)) #|>

datatable(probesx |>
    dplyr::select(-c(1:5, 7:8, 10:12)) |>
    dplyr::filter(adj.P.Val < 0.01) |>
        arrange(desc(ATL_AC))) #|>)

negtcr_gbacc = annot |>
    dplyr::filter(SYMBOL %in% neg_tcr) |>
    pull(GB_ACC)

#--------------------------------------------------------
# Part 8: Heatmap
#--------------------------------------------------------
# Create a heatmap of the top differentially expressed genes
# Filter the top genes based on your criteria
annotationdf =data.frame(
    row.names = pd7$ID,
    Subtype=pd7$atl_subtype2)
class(heatmapgenes)
heatmapgenes = gse1$A[gse1$genes$SYMBOL %in% neg_tcr,] |>
    as.data.frame() |>
    tibble::rownames_to_column("GB_ACC") |>
    left_join(annot[,c("GB_ACC", "SYMBOL")], by=c("GB_ACC")) |>
    tibble::column_to_rownames("SYMBOL") |>
    dplyr::select(-GB_ACC) |>
    as.matrix()

pheatmap(
    mat = heatmapgenes,
    scale = "row",
    annotation_col = annotationdf,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellheight=15,
    main = "TCR Negative Regulators",
    filename = "2504_negativeTCRregulators.png",
    width = 10,  # Nature standard single-column
    height = 7,
    units = "in",
    family = "Arial"  # Embed font
)

#--------------------------------------------------------
# Part 9: Gene Annotation
#--------------------------------------------------------
# Gene annotation
library(org.Hs.eg.db)
library(AnnotationDbi)
neg_tcr
current_symbol <- mapIds(
    org.Hs.eg.db,
    keys = neg_tcr,        # Input symbol/alias
    column ="SYMBOL",    # Output column (current symbol)
    keytype = "ALIAS",    # Input type
    multiVals = "CharacterList"   # Return first match if multiple exist
)

#--------------------------------------------------------
# Part 10: Volcano plot
#--------------------------------------------------------

top_probes |>
    mutate(sign = ifelse(adj.P.Val < 0.01 & abs(ATL_AC) > 1.5, "Sign", "No")) |>
    dplyr::mutate(adp = -log10(adj.P.Val)) |>
    tidyplots::tidyplot(x = ATL_AC, y = adp, color = sign) |>
    tidyplots::add_data_points(alpha = .5) |>
    tidyplots::add_data_labels_repel(label = SYMBOL, data = filter_rows(sign == "Sign"), color = "black") |>
    tidyplots::adjust_x_axis_title("Log2(Fold Change)") |>
    tidyplots::adjust_y_axis_title("-Log10(Adjusted Pvalue)") |>
    tidyplots::remove_legend() |>
    tidyplots::save_plot("ATLvATLc.png",
        bg = "transparent")



#--------------------------------------------------------
# Part 4B: Manipulating targets
#--------------------------------------------------------


#--------------------------------------------------------
# Part 5: QC
#--------------------------------------------------------
# Convert targets to tibble for easy joining

## Apply the function to the data
gse2 = agilent_two_color_qc(agilent_data2)

##### PCA
dim(gse2)
boxplot(gse2$M, main="Normalized Negative Controls")

#--------------------------------------------------------
# Part 6B: Merging pending, it is necessary to merge and then to sva - batch, 
#Extracting 

merged_data <- merge(
  as.data.frame(gse1$M),
  as.data.frame(gse2$M),
  by.x = "genes$Name",
  by.y = "genes$probename",
  all = FALSE  # Keep only matching genes
)
class(gse1)

head(gse2$genes)
