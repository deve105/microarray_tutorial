#--------------------------------------------------------
# Title: Microarray Data Analysis with GEOquery and limma


#GSE65823 Reanalyzed by: GSE86362 GSE119087

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
id <- "GSE19069" # GPL4133

# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}
# Obtaining Metadata
meta <- getGEO(id, GSEMatrix = TRUE, destdir = ".temp")
meta <- meta[[1]]

tail(pData(meta), 5)
dim(meta)

# Metadata wrangling
tail(pData(meta)$"age:ch1", 5)
names(pData(meta))
pd <- pData(meta) |>
    tibble::rownames_to_column("ID") |>
    dplyr::select("age:ch1", "gender:ch1", "cell type:ch1", ID, supplementary_file, title) |>
    dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
    dplyr::rename(age = "age:ch1", gender = "gender:ch1", type="cell type:ch1")  |>
    mutate(ptcl = case_when(
        str_detect(title, "Adult") ~ "ATL",
        str_detect(title, "ALK-NEG") ~ "ALKneg",
        str_detect(title, "ALK-POS") ~ "ALKpos",
        str_detect(title, "Angioimmunoblastic") ~ "AITL",
        str_detect(title, "Peripheral") ~ "PTCL",
        str_detect(title, "Activated CD4") ~ "CD4a",
        str_detect(title, "Activated CD8") ~ "CD8a",
        str_detect(title, "Resting CD4") ~ "CD4r",
         str_detect(title, "Resting CD8") ~ "CD8r",
        str_detect(title, "JURKAT") ~ "JURKAT",
         str_detect(title, "MOLT") ~ "MOLT3",
         str_detect(title, "PEER") ~ "PEER",
         str_detect(title, "Sezary") ~ "Sezary",
        str_detect(title, "prolymphocytic") ~ "PLL",
         str_detect(title, "Lymphoid") ~ "Lymph",
        .default = NA)) 


## Download all the files in the .temp of this environment
for (i in 1:length(pd$supplementary_file)) {
    url <- pd$supplementary_file[i]
    destfile <- file.path(paste0(".temp/", id, "/", pd$file[i]))
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
con <- gzfile(file.path("XXX"))
file_lines <- readLines(con, n=1000)
close(con)

# Subsettting files (here we have two different platforms) 
table(pd$ptcl)

PTCL = c("ALKneg", "ALKpos", "ATL", "CD4a", "AITL", "PTCL")

pd7 = pd |>
    filter(ptcl%in%PTCL) 
    
dim(pd7)

#### gene Annotation
gpl2 <- getGEO("GPL570")
annot <- Table(gpl2)[,c("ID", "GB_ACC", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID")] #|>
    tibble::column_to_rownames("ID")

# Re-read the raw files, this was a quantarray (described in metadata)
affy_data = ReadAffy(
    filenames=pd7$file,
    sampleNames = pd7$ID,
    celfile.path = file.path(".temp", id)
)

affy_rma = rma(affy_data)

methods(class=class(affy_data))
protocolData(affy_data)


head(affy_data)
#--------------------------------------------------------
# Part 6: DE
#--------------------------------------------------------
# Create design matrix based on final_subtype
design <- model.matrix(~ 0 + factor(pd7$ptcl))
rownames(design) <- pd7$ID
colnames(design) <- levels(factor(pd7$ptcl))
head(design, 20)

# Make contrasts (adjust based on your comparisons)
levels(factor(pd7$ptcl))
contrast.matrix <- makeContrasts(
    ATL_AITL = ATL - AITL,
    ATL_ALKpos = ATL - ALKpos,
    ATL_ALKneg = ATL - ALKneg,
    ATL_PTCL = ATL - PTCL,
    ATL_CD4a = ATL - CD4a,
    levels = design
)
fit2 <- lmFit(affy_rma, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

top_probes <- topTable(fit2, number = Inf, adjust.method = "BH" ,  p.value = 0.01)



head(annot$ID)
top_probes |>
    tibble::rownames_to_column("ID") |>
    left_join(annot, by="ID") |>
    arrange((ATL_CD4a)) |>
    head(50)

#--------------------------------------------------------
# Part 4: Manipulating targets
#--------------------------------------------------------
# Convert targets to tibble for easy joining
agilent_data$targets = agilent_data$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")



# Add gene symbols to agilent_data

head(agilent_data$gene)
agilent_data$genes <- agilent_data$genes |>
  dplyr::left_join(annot, by=c("Name"="GB_ACC")) |>
  mutate(is_control=factor(ifelse(is.na(SYMBOL) & is.na(GENE_NAME), "control", "gene")))

table(agilent_data$gene$is_control)

#agilent_data$genes$GENE_SYMBOL <- annot$SYMBOL[match(agilent_data$genes$Name, annot$GB_ACC)]

table(agilent_data$genes$is_control)
