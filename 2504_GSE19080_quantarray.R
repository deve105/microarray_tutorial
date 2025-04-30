#--------------------------------------------------------
# Title: Microarray Data Analysis with GEOquery and limma

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

methods(class=class(meta))

meta <- meta[[1]]
meta
head(fData(meta), 2)
head(exprs(meta), 2)
pData(meta)
fData(meta)
colnames(exprs(meta))
dim(meta)

head(pData(meta), 5)
names(pData(meta))
pData(meta)$characteristics_ch1

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
    other.columns = list(
        Flag = "Ignore Filter"))

#--------------------------------------------------------
# Part 4: Manipulating targets
#--------------------------------------------------------
# Convert targets to tibble for easy joining
agilent_data$targets = agilent_data$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")
agilent_data$targets

#### gene Annotation
gpl2 <- getGEO("GPL9686")
head(fData(meta)[,c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")])
annot <- Table(gpl2)[, c("SYMBOL", "GENE_NAME", "GB_ACC")]

# Add gene symbols to agilent_data
head(agilent_data$gene)
agilent_data$genes <- agilent_data$genes |>
  dplyr::left_join(annot, by=c("Name"="GB_ACC")) |>
  mutate(is_control=factor(ifelse(is.na(SYMBOL) & is.na(GENE_NAME), "control", "gene")))

table(agilent_data$gene$is_control)

#agilent_data$genes$GENE_SYMBOL <- annot$SYMBOL[match(agilent_data$genes$Name, annot$GB_ACC)]

table(agilent_data$genes$is_control)

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
    agilent_data <- limma::backgroundCorrect(agilent_data, method="normexp", offset=20) 
    # 5. Normalization
    agilent_data <- limma::normalizeWithinArrays(agilent_data, method="loess")
    agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")

    # 6. technical replicates reduction by probename
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$Name)
    agilent_data <- agilent_data[agilent_data$genes$is_control=="gene", ]
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$GENE_SYMBOL)
    agilent_data$genes <- droplevels(agilent_data$genes)
    return(agilent_data)
}

## Apply the function to the data
gse1 = agilent_two_color_qc(agilent_data)

##### PCA
dim(gse1)

boxplot(gse1$M, main="Normalized Negative Controls")

pd7

pca <- prcomp(t(gse1$M))
pca
rownames(pca$x)

pca = pca$x |>
    as.data.frame() 
pca |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(pd7, by="ID") |> 
    tidyplots::tidyplot(x=PC1, y=PC2, color=atl_subtype) |>
    tidyplots::add_data_points()

#--------------------------------------------------------
# Part 6: DE
#--------------------------------------------------------
# Create design matrix based on final_subtype
pd7$atl_subtype2

design <- model.matrix(~ 0 + factor(pd7$atl_subtype2))
design
rownames(design) <- pd7$ID
colnames(design) <- levels(factor(pd7$atl_subtype2))
head(design, 20)

# Make contrasts (adjust based on your comparisons)
levels(factor(pd7$atl_subtype2))
contrast.matrix <- makeContrasts(
    ATL_AC = ATL - AC,
    levels = design
)
fit2 <- lmFit(gse1, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf, adjust.method = "BH") #,  p.value = 0.01)
head(top_probes)

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

colnames(top_probes)

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
# Part 3B: HD - Capturing targets
#--------------------------------------------------------
# Subsettting files (here we have two different platforms) 

pd7_2 = pd |>
    filter(ID%in%c(paste0("GSM4723", 74:81)))

agilent_data2 <- read.maimages(
    files = file.path(".temp", id, pd7_2$file),
    source = "genepix",
    green.only = FALSE,
    names = pd7_2$ID)
#--------------------------------------------------------
# Part 4B: Manipulating targets
#--------------------------------------------------------
# Convert targets to tibble for easy joining
agilent_data2$targets = agilent_data2$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7_2, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")

#### gene Annotation
gpl2 <- getGEO("GPL9686")
annot <- Table(gpl2)[, c("SYMBOL", "GENE_NAME", "GB_ACC")]

# Add gene symbols to agilent_data
head(agilent_data2$gene)
agilent_data2$genes <- agilent_data2$genes |>
    tidyr::separate(
    Name,
    into = c("GENE_SYMBOL", "Annot", "probename"),
    sep = ":")  |> 
  dplyr::left_join(annot, by=c("probename"="GB_ACC")) |>
  mutate(is_control=factor(ifelse(is.na(SYMBOL) & is.na(GENE_NAME), "control", "gene")))

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
