## Loading packages
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots, affy, oligoClasses, testit,
    hgu133plus2.db,
    org.Hs.eg.db, reshape2, survminer, arrayQualityMetrics, testit, biomaRt, AgiMicroRna, Biobase, marray, matrixStats
)
## Function
isLog2Transformed <- function(data) {
    qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    shouldBeLogged <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    return(!shouldBeLogged)
}

## First obtain metadata
id="GSE40639"
GSE40639_meta <- getGEO(id, GSEMatrix = FALSE, destdir = ".temp")
GSE40639_meta <- GSE40639_meta[[1]]
Meta(GSE40639_meta)
methods(class=class(GSE40639_meta))
table(GSMList(GSE40639_meta)[1])
# Explore columns
class(GSE40639_meta)
varLabels(GSE40639_meta)
head(pData(GSE40639_meta))
head(fData(GGSE40639_meta))

# Data wrangling
pd <- pData(GSE40639_meta) |>
    tibble::rownames_to_column("ID") |>
    dplyr::select(title, "characteristics_ch1.1", "characteristics_ch1.3", "characteristics_ch1.4", "characteristics_ch1.5", ID) |>
    #dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
    dplyr::rename(
        region = "characteristics_ch1.1", 
        malignancy = "title", 
        gender= "characteristics_ch1.3",
        age = "characteristics_ch1.4",
        stage="characteristics_ch1.5") |>
        mutate(
            across(
        c(region, gender, age, stage),
        ~ stringr::str_remove(., "^.*: ")
        ),
        malignancy=ifelse(str_detect(malignancy, "mycosis"), "MF", "ATLL"))
    
pd    

pacman::p_load(beadarray)
beadData <- read.ilmn(path = ".temp/GSE40639/")
beadData

agilent_data <- read.maimages(
    files = file.path(".temp", id, pd7$file),
    source = "quantarray",
    green.only = FALSE,
    names = pd7$ID,
    other.columns = list(
        Flag = "Ignore Filter"))