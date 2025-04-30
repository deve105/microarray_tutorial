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
GSE40639_meta <- getGEO(id, GSEMatrix = TRUE, destdir = ".temp")

GSE40639_meta <- GSE40639_meta[[1]]
methods(class=class(GSE40639_meta))

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


## Download Supp files
# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}

## Downloading supplementary raw data files
getGEOSuppFiles(id, baseDir = paste0(".temp/", id), makeDirectory = FALSE)


## Loop for unzipping the files
files_to_process <- rownames(raw)[2:3]
for (file in files_to_process) {
  dest_path <- gsub("\\.gz$", "", paste0(".temp/", id, "/", basename(file)))
  gunzip(file, remove = TRUE, overwrite = TRUE, destname = dest_path)
}

## untar the raw data:
untar(paste0(rownames(raw)[1]), exdir = paste0(".temp/", id, "/GSE40639"))
gunzip(".temp/GSE40639/GSE40639/GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz")

## Read the data
illumina_data <- read.ilmn(
    ".temp/GSE40639/GSE40639_non-normalized_dermis.txt",
    expr = "SAMPLE",
    probeid = "ID_REF",
    other.columns = c("Detection"),

)

illumina_data_e <- read.ilmn(
    ".temp/GSE40639/GSE40639_non-normalized_epidermis.txt",
    expr = "SAMPLE",
    probeid = "ID_REF",
    other.columns = c("Detection"),

)


colnames(illumina_data) = c("GSM998616","GSM998617","GSM998618","GSM998619","GSM998620","GSM998621")
colnames(illumina_data_e) = c("GSM998606","GSM998607","GSM998608","GSM998609","GSM998610","GSM998611","GSM998612","GSM998613","GSM998614","GSM998615")
illumina_data
illumina_data_e

all(rownames(illumina_data) == rownames(illumina_data_e))
finaldb=cbind(illumina_data, illumina_data_e)
######### Joining
combined_exprs <- cbind(illumina_data$E, illumina_data_e$E)
head(combined_exprs)
table(rownames(illumina_data$E) == rownames(illumina_data_e$E))
combined_detection <- cbind(illumina_data$other$Detection, illumina_data_e$other$Detection)
table(rownames(illumina_data$other$Detectio) == rownames(illumina_data_e$other$Detectio))
dim(combined_exprs) == dim(combined_detection)


sample_ids = colnames(combined_exprs)
pd_matched = pd[match(sample_ids, pd$ID), ]
rownames(pd_matched) = sample_ids
pd_matched
pheno_data = new("AnnotatedDataFrame", data = pd_matched)

annot_matched = annot[match(rownames(combined_exprs), annot$ID), ]
rownames(annot_matched) = annot_matched$ID
head(annot_matched)
head(combined_exprs)
head(combined_detection)


feature_data = new("AnnotatedDataFrame", data = annot_matched)

# First, check your probe types
table(fData(gse)$Probe_Type, useNA = "always")

### Creating an expression set object
gse = new("ExpressionSet",
    exprs = combined_exprs,
    phenoData = pheno_data,
    featureData = feature_data,
    annotation = "GPL6884"
)

assayDataElement(gse, "Detection") <- combined_detection
rownames(gse)

# Normalizing
gse_normalized <- neqc(
  x = exprs(gse),
  detection.p = assayDataElement(gse, "Detection"),
  offset = 20
)
head(assayDataElement(gse, "Detection"))

anyNA(combined_exprs$others$Detection)
dim(combined_exprs)

# Platform annotation
platform <- "GPL6884"
GPL6884= getGEO(platform, destdir = ".temp")
annot = Table(GPL6884)[,c("ID","Probe_Type", "Symbol", "Transcript", "ILMN_Gene")]
colnames(Table(GPL6884))

# Check if the annotation is correct
combined_exprs$Annot = data.frame(
    ID=rownames(combined_exprs)) |>
    left_join(annot, by = "ID") |>
    column_to_rownames("ID") 


annot_matched =  annot_matched |>
    dplyr::mutate(ProbeType=ifelse(Probe_Type == "A", TRUE, FALSE))

probess = data.frame(
    row.names=rownames(annot_matched),
    ProbeType=annot_matched$ProbeType)
dim(combined_exprs)
dim(combined_detection)
dim(probess)
combined_qc = neqc(x=combined_exprs,  offset = 20, detection.p=combined_detection, status=probess)

combined_qc
table(annot$Probe_Type)

failed_probes <- rowSums(combined_exprs$other$Detection > 0.05) > ncol(combined_exprs)/2
sum(failed_probes) 

dim(combined_exprs$other$Detection)


#Using Lumi
pacman::p_load(lumi)

gse
raw_data_summarized <- lumiExpresso(gse, 
                                    bg.correct = TRUE, 
                                    bgcorrect.param = list(method = 'bgAdjust'), 
                                    variance.stabilize = FALSE,  # No scaling
                                    varianceStabilize.param = list(), 
                                    normalize = FALSE, 
                                    normalize.param = list(), 
                                    QC.evaluation = TRUE, 
                                    QC.param = list(), 
                                    verbose = TRUE)
lumiB(exprs(gse), method = "bgAdjust")


library(beadarray)


sampleSheetFile <- paste0(chipPath,"/sampleSheet.csv")
data <- readIllumina(
    dir="./output",
    sampleSheet="GSE40639_non-normalized_dermis.txt",
    illuminaAnnotation="Humanv3"
    )
