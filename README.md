# ExoNet - MS/MS analysis

Alignment and statistical analysis of LC/ESI-MS/MS data from human PDAC subpopulations.

## Project description

Characterization of the protein composition of exosomes and cells of five distinct subpopulations from four human PDAC cell lines (BxPC-3, PANC-1, T3M4, and MIA PaCa-2), via liquid chromatography electrospray ionization tandem mass spectrometry (LC/ESI–MS/MS). 

Independent Data Acquisition (IDA) using Analyst TF 1.7 software from a Triple TOF 5600 System (SCIEX, USA). The data were processed using PeakView® 2.2 Software (SCIEX, Foster City, CA), and [MS-GF+ software](https://github.com/MSGFPlus/msgfplus) as a search engine.

## Data availability

The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD023529 and 10.6019/PXD023529. 

## Data analysis

Acquired data were aligned using MSGF+ against human reference proteome (UP000005640_9606.fasta) obtained from [EMBL-EBI repository](https://www.ebi.ac.uk/reference_proteomes/), as described in [msgfplus_alignment.sh](Scripts/msgfplus_alignment.sh) script. QC check, filtering, and final sets were obtained as described at [data_analysis.R](Scripts/data_analysis.R). Here two main sets were generated derived from same alignment strategy:

- Co-ocurrency matrix, indicating which proteins have been indepedently detected in each sample.

- Quantitative matrix, peptide quantitation of proteins found in all the set of samples.

## Citation

Ruivo CF, Bastos N, Adem B, et al Extracellular Vesicles from Pancreatic Cancer Stem Cells Lead an Intratumor Communication Network (EVNet) to fuel tumour progression Gut Published Online First: 10 January 2022. doi: 10.1136/gutjnl-2021-324994


