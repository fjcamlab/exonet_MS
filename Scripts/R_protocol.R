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
library(DT)
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


