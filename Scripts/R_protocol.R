####################################################################################################
### Statistical analysis of MS-MS data
####################################################################################################
# version January 2020



## Proteomic libraries
library(ProteoMM)
library(MSnbase)
library(mzID)
library(MSnID)
library(mzR)

## Additional libraries
library(readr)
library(Homo.sapiens)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(Rtsne)
library(DESeq2)
library(openxlsx)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(viridis)




######## -------------------------------------------------------------------------------------------
## 0: Global variables and sample information
######## -------------------------------------------------------------------------------------------

path_data <- "processed/Align_withCrap_MSGF/"
path_analysis_osiris <- "/home/fjc38/DATA/mntosiris3/MS/SoniaMelo_TMT/"
protein_fasta <- "analysis/fasta/UP000005640_9606.fasta" ## EMBL-EBI
contaminant_fasta <- "analysis/fasta/crap.fasta" ## cRAP repository

## Sample information
sample_info <- read_csv("Sample_info.csv", col_types = "ifffc")[,-1]
sample_info <- sample_info[with(sample_info, order(Cell_line, Subpopulation)),]




######## -------------------------------------------------------------------------------------------
## 1: Protein information
######## -------------------------------------------------------------------------------------------

## Protein information
protein_info <- fread(paste0("grep '^>' ", protein_fasta), 
                      header = FALSE, sep = "|", col.names = c("class", "UNIPROT", "Description"), 
                      colClasses = rep("character", 3))
protein_info <- DataFrame(protein_info[,-3], do.call(rbind, strsplit(
  protein_info$Description, "MAN | ..=")))
colnames(protein_info)[-c(1:2)] <- c("UNIPROT_ID", "Description", "OS", "OX", "GN", "PE", "SV")
protein_info$UNIPROT_ID <- paste0(protein_info$UNIPROT_ID, "MAN")

## Adding Entrez IDs
a <- AnnotationDbi::select(Homo.sapiens, keys = protein_info$UNIPROT, 
                           keytype = "UNIPROT", columns = "ENTREZID")
protein_info$ENTREZID <- a[match(protein_info$UNIPROT, a$UNIPROT), "ENTREZID"]


# Subcellular location from Human Protein Atlas (downloaded 08-05-2020)
subcell <- read_tsv("analysis/additional_data/subcellular_location.tsv")
protein_info <- cbind(protein_info, subcell[match(protein_info$GN, subcell$`Gene name`), 
                                            c(4:5)])
protein_info$class <- sub(">", "", protein_info$class)


## Contaminants information
contaminants_info <- fread(paste0("grep '^>' ", contaminant_fasta), 
                           header = FALSE, sep = "|", 
                           col.names = c("class", "UNIPROT", "Description"), 
                           colClasses = rep("character", 3))
contaminants_info <- DataFrame(contaminants_info[,-3])




######## -------------------------------------------------------------------------------------------
## 2: Reading MGSF+ files from searching/alignment
######## -------------------------------------------------------------------------------------------

## Defining path to 39 files prev. aligned with MSGF+ 
outl <- !(sample_info$Subpopulation == "EpCAM+" & sample_info$Cell_line == "MIAPaCa2") ## outliers
names <- c(sample_info[,"File"])[[1]][outl]
conv <- c(sample_info[,"Cell_line"])[[1]][outl]
files <- paste0(path_data, names, ".mzid")

## 2.1 -- Reading files separately: to define protein lists
msnid_l <- lapply(files, function(x) {
  print(x)
  read_mzIDs(msnid, x)
})
names(msnid_l) <- names

## 2.2 -- Reading in a single object: for quantitation and differential analysis
msnid <- MSnID(".") # provide working directory
msnid_m <- read_mzIDs(msnid, files, backend = "mzR")




######## -------------------------------------------------------------------------------------------
## 3: Quality checks
######## -------------------------------------------------------------------------------------------

ppm <- lapply(msnid_l, function(x) {
  ppm <- apply_filter(x, "abs(mass_measurement_error(x)) < 25")
  mass_measurement_error(ppm)
})
names(ppm) <- names
ppm <- melt(ppm)
colnames(ppm) <- c("value", "sample")
ppm <- data.frame(ppm, sample_info[match(ppm$sample, names), -4])

qc_plot <- ggplot(ppm, aes(x=value, fill = Subpopulation, col = Subpopulation)) + 
  scale_fill_viridis_d() + scale_color_viridis_d() + 
  xlab("mass_measurement_error") + ylab("number peptides") + 
  geom_histogram(alpha = 0.5, bins = 100) + theme_bw() + 
  geom_vline(xintercept = 0, lty = 2) +
  ggtitle("mass_measurement_error") + 
  facet_wrap(~Fraction+Cell_line, scales = "free", ncol = 4)

## Plot
# qc_plot


######## -------------------------------------------------------------------------------------------
## 4: Filter for protein lists: FDR 1% at peptide level
######## -------------------------------------------------------------------------------------------


## 4.1 -- Scoring individual samples
scored <- lapply(msnid_l, function(x) {
  x$absParentMassErrorPPM <- abs(mass_measurement_error(x))
  x$msmsScore <- -log10(x$`MS-GF:SpecEValue`)
  x$contaminant <- ifelse(unlist(lapply(strsplit(x$accession, '\\|'), length)) == 3, 
                          "protein", "contaminant")
  x
})

filtObj <- MSnIDFilter(scored[[1]]) 
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)

# Define level (PSM, peptide, accession - protein)
msnid_filtered_l <- lapply(scored, function(x) {
  filtObj.grid <- optimize_filter(filtObj, x, fdr.max=0.01,
                                  method="Nelder-Mead", level="peptide",
                                  n.iter=500)
  # applying fine filter
  x <- apply_filter(x, filtObj.grid)
  # removing decoys
  x <- apply_filter(x, "isDecoy == FALSE")
  # # removing keratins
  x <- apply_filter(x, "!grepl('[Kk]eratin', description)")
  # # removing contaminants
  x <- apply_filter(x, "contaminant != 'contaminant'")
  # # removing serum
  apply_filter(x, "!grepl('Serum albumin', description)")
})


## 4.2 -- Scoring and filtering all samples
msnid_m$absParentMassErrorPPM <- abs(mass_measurement_error(msnid_m))
msnid_m$msmsScore <- -log10(msnid_m$MS.GF.SpecEValue)
msnid_m$contaminant <- ifelse(unlist(lapply(strsplit(msnid_m$accession, '\\|'), length)) == 3, 
                              "protein", "contaminant")

filtObj <- MSnIDFilter(msnid_m) 
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)

# Define level (PSM, peptide, accession - protein)
filtObj.grid <- optimize_filter(filtObj, msnid_m, fdr.max=0.01,
                                method="Nelder-Mead", level="peptide",
                                n.iter=500)
# applying fine filter
msnid_filtered_m <- apply_filter(msnid_m, filtObj.grid)
# removing decoys
msnid_filtered_m <- apply_filter(msnid_filtered_m, "isDecoy == FALSE")
# removing keratins
msnid_filtered_m <- apply_filter(msnid_filtered_m, "!grepl('[Kk]eratin', description)")
# removing contaminants
msnid_filtered_m <- apply_filter(msnid_filtered_m, "contaminant != 'contaminant'")
# removing serum
msnid_filtered_m <- apply_filter(msnid_filtered_m, "!grepl('Serum albumin', description)")




######## -------------------------------------------------------------------------------------------
## 5: Detected proteins: FDR 1% - peptide level
######## -------------------------------------------------------------------------------------------

## Proteins per sample at 1% (peptide level)
proteins_sample <- lapply(msnid_filtered_l, function(x) 
  unique(sapply(proteins(x), function(y) 
    str_match(y, pattern = "(?<=\\|)[^|]++(?=\\|)")[1,1])))

prot_sample_info <- sapply(proteins_sample, function(x) 
  protein_info[match(x, protein_info$UNIPROT),])

## Writing file with all identified proteins in each samples
write.xlsx(prot_sample_info, 
           file = "IdentifiedProteins_FDR0.01_F.xlsx")


## Matrix defining presence or absence of every protein found per sample
co <- data.frame(acast(melt(lapply(proteins_sample, cbind)), 
                       value ~ L1))

coocurrency <- co
coocurrency[!is.na(coocurrency)] <- 1
coocurrency[is.na(coocurrency)] <- 0
coocurrency <- apply(coocurrency, 2, as.numeric)
rownames(coocurrency) <- rownames(co)
colnames(coocurrency) <- sort(names(proteins_sample))

write.table(data.frame(UNIPROT = rownames(coocurrency), coocurrency), 
            file = "BinaryMatrix_MS.tsv", row.names = FALSE)




######## -------------------------------------------------------------------------------------------
## 6: Quantitation - Differential analysis on protein counts (DESeq2)
######## -------------------------------------------------------------------------------------------


## 6.1 -- Summarising peptide to protein counts
msnset <- as(msnid_filtered_m, "MSnSet")
sampleNames(msnset) <- sub("\\.mzid$", "", sampleNames(msnset))

msnset_prot <- combineFeatures(msnset,
                               fData(msnset)$accession,
                               redundancy.handler="unique",
                               method="sum",
                               cv=FALSE)

k <- data.frame(sample_info[match(sampleNames(msnset_prot), sample_info$File),])
k$contrast <- factor(ifelse(k$Subpopulation %in% c("CD133+", "CD24+CD44+"), "B", "A"))
rownames(k) <- sampleNames(msnset_prot)

pData(msnset_prot) <- k


## 6.2 -- DESeq2 analysis on Exosomes fraction

e <- exprs(msnset_prot)
pData(msnset_prot)$Subpopulation <- gsub("\\+", "p", pData(msnset_prot)$Subpopulation)
pData(msnset_prot)$Subpopulation <- gsub("\\-", "n", pData(msnset_prot)$Subpopulation)
pData(msnset_prot)$Subpopulation <- gsub("\\:", "_", pData(msnset_prot)$Subpopulation)
pData(msnset_prot)$Subpopulation <- factor(pData(msnset_prot)$Subpopulation)

## DESeq2 object for all samples
dds <- DESeqDataSetFromMatrix(countData = e[rowSums(e) > 0, ], 
                              colData = pData(msnset_prot), 
                              design = ~ 0 + Subpopulation + Fraction + Cell_line)


## 6.3 -- CSCs vs NCSCs exososomes
exo <- pData(msnset_prot)$Fraction == "Exosomes"
dds_exo <- DESeqDataSetFromMatrix(countData = e[rowSums(e[,exo]) > 0, exo], 
                              colData = pData(msnset_prot)[exo,-1], 
                              design = ~ contrast + Cell_line)
dds_exo <- DESeq(dds_exo, fitType='local')

## Fraction 
res_dds <- DESeq2::results(dds_exo, contrast = c("contrast", "B", "A"))
res_dds <- data.frame(data.frame(ID = rownames(res_dds), res_dds), 
                      protein_info[match(sub(".*\\|", "", rownames(res_dds)), 
                                         protein_info$UNIPROT_ID),])

write.table(res_dds, file = "CSCvsNCSCexo_DESeq2.tsv", row.names = FALSE)


                           
# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] MSnID_1.22.0        mzID_1.26.0         MSnbase_2.14.2      ProtGenerics_1.20.0 S4Vectors_0.26.1    mzR_2.22.0          Biobase_2.48.0      BiocGenerics_0.34.0
# [9] ProteoMM_1.6.0      Rcpp_1.0.5         
# 
# loaded via a namespace (and not attached):
#   [1] R.utils_2.10.1              tidyselect_1.1.0            RSQLite_2.2.1               AnnotationDbi_1.50.3        htmlwidgets_1.5.2          
# [6] grid_4.0.3                  BiocParallel_1.22.0         munsell_0.5.0               codetools_0.2-17            preprocessCore_1.50.0      
# [11] DT_0.16                     colorspace_1.4-1            knitr_1.30                  rstudioapi_0.11             ggsignif_0.6.0             
# [16] GenomeInfoDbData_1.2.3      bit64_4.0.5                 coda_0.19-4                 vctrs_0.3.4                 generics_0.0.2             
# [21] xfun_0.18                   BiocFileCache_1.12.1        R6_2.4.1                    doParallel_1.0.15           GenomeInfoDb_1.24.2        
# [26] clue_0.3-57                 locfit_1.5-9.4              bitops_1.0-6                fgsea_1.14.0                DelayedArray_0.14.1        
# [31] assertthat_0.2.1            scales_1.1.1                gtable_0.3.0                affy_1.66.0                 rlang_0.4.8                
# [36] genefilter_1.70.0           scatterplot3d_0.3-41        GlobalOptions_0.1.2         splines_4.0.3               rtracklayer_1.48.0         
# [41] rstatix_0.6.0               impute_1.62.0               broom_0.7.1                 BiocManager_1.30.10         yaml_2.2.1                 
# [46] reshape2_1.4.4              abind_1.4-5                 made4_1.61.0                GenomicFeatures_1.40.1      backports_1.1.10           
# [51] tools_4.0.3                 statnet.common_4.4.1        ggplot2_3.3.2               affyio_1.58.0               deco_1.4.0                 
# [56] ellipsis_0.3.1              gplots_3.1.0                RColorBrewer_1.1-2          ggridges_0.5.2              plyr_1.8.6                 
# [61] progress_1.2.2              zlibbioc_1.34.0             purrr_0.3.4                 RCurl_1.98-1.2              prettyunits_1.1.1          
# [66] ggpubr_0.4.0                openssl_1.4.3               GetoptLong_1.0.3            sfsmisc_1.1-7               SummarizedExperiment_1.18.2
# [71] haven_2.3.1                 ggrepel_0.8.2               cluster_2.1.0               magrittr_1.5                data.table_1.13.0          
# [76] sna_2.6                     openxlsx_4.2.2              circlize_0.4.10             pcaMethods_1.80.0           R.cache_0.14.0             
# [81] matrixStats_0.57.0          hms_0.5.3                   evaluate_0.14               xtable_1.8-4                XML_3.99-0.5               
# [86] rio_0.5.16                  readxl_1.3.1                IRanges_2.22.2              gridExtra_2.3               shape_1.4.5                
# [91] compiler_4.0.3              biomaRt_2.45.9              tibble_3.0.4                KernSmooth_2.23-17          ncdf4_1.17                 
# [96] crayon_1.3.4                R.oo_1.24.0                 htmltools_0.5.0             segmented_1.3-0             tidyr_1.1.2                
# [101] geneplotter_1.66.0          DBI_1.1.0                   dbplyr_1.4.4                ComplexHeatmap_2.4.3        MASS_7.3-53                
# [106] rappdirs_0.3.1              BiocStyle_2.16.1            Matrix_1.2-18               ade4_1.7-15                 car_3.0-10                 
# [111] vsn_3.56.0                  R.methodsS3_1.8.1           gdata_2.18.0                igraph_1.2.6                GenomicRanges_1.40.0       
# [116] forcats_0.5.0               pkgconfig_2.0.3             GenomicAlignments_1.24.0    foreign_0.8-80              MALDIquant_1.19.3          
# [121] xml2_1.3.2                  foreach_1.5.0               annotate_1.66.0             XVector_0.28.0              viper_1.22.0               
# [126] stringr_1.4.0               digest_0.6.25               rle_0.9.2                   Biostrings_2.56.0           rmarkdown_2.4              
# [131] cellranger_1.1.0            fastmatch_1.1-0             curl_4.3                    kernlab_0.9-29              Rsamtools_2.4.0            
# [136] gtools_3.8.2                rjson_0.2.20                lifecycle_0.2.0             carData_3.0-4               network_1.16.1             
# [141] askpass_1.1                 limma_3.44.3                pillar_1.4.6                lattice_0.20-41             httr_1.4.2                 
# [146] survival_3.2-7              glue_1.4.2                  zip_2.1.1                   png_0.1-7                   iterators_1.0.12           
# [151] bit_4.0.4                   class_7.3-17                stringi_1.5.3               mixtools_1.2.0              blob_1.2.1                 
# [156] DESeq2_1.28.1               caTools_1.18.0              memoise_1.1.0               dplyr_1.0.2                 e1071_1.7-3
                           
                           
