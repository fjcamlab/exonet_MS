###################################################################### 
### MS-MS analysis to Sonia Melo group
###
# Author: FJ Campos-Laborie
# MS-GF+ Release (v2020.03.14) (14 March 2020)

## Data folder
DATA=./

##### MSGF+ search
## See tutorial in https://msgfplus.github.io/msgfplus/MSGFPlus.html
MSGF=msgfplus/MSGFPlus.jar

cat crap.fasta UP000005640_9606.fasta > combined_9606_crap.fasta
DB=combined_9606_crap.fasta 

## create a configuration file (same for all files)
CONF=Qual_conf_file.txt

# Search Parameters:
#         PrecursorMassTolerance: 25.0 ppm
#         IsotopeError: -1,2
#         TargetDecoyAnalysis: true
#         FragmentationMethod: As written in the spectrum or CID if no info
#         Instrument: TOF
#         Enzyme: Tryp
#         Protocol: TMT
#         NumTolerableTermini: 2
#         MinPepLength: 6
#         MaxPepLength: 50
#         MinCharge: 2
#         MaxCharge: 5
#         NumMatchesPerSpec: 1
#         MaxMissedCleavages: 2
#         MaxNumModsPerPeptide: 3
#         ChargeCarrierMass: 1.00727649 (proton)
#         MinNumPeaksPerSpectrum: 10
#         NumIsoforms: 128
# Post translational modifications in use:
#         Fixed (static):     Carbamidomethyl on C (+57.0215)
#         Variable (dynamic): Glu->pyro-Glu on E at the peptide N-terminus (-18.0106)
#         Variable (dynamic): Gln->pyro-Glu on Q at the peptide N-terminus (-17.0265)
#         Variable (dynamic): Oxidation on M (+15.9949)


## Align samples to database
mkdir $DATA/aligned

for file in $DATA/processed/*.mgf; do
  
  echo "Searching with crap matches in file: $file"
  java -Xmx3500M -jar $MSGF -s $file -d $DB -ti -1,2 -tda 1 \
    -maxMissedCleavages 2 -protocol 4 -conf $CONF
  
done

mv *.mzid $DATA/aligned
