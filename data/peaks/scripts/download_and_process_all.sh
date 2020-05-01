# On Linux
export PROMPT_COMMAND="echo -n \[\$(date +%H:%M:%S)\]\ "
cd /mnt/d/projects/HOP-HMM/data/peaks/scripts
yes | sudo apt update
yes | sudo apt install bedtools


# On Windows
cygwin
cd /cygdrive/d/projects/HOP-HMM/data/peaks/scripts


# On Windows and Linux

################################
# Bedtools2 Installation
################################
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
pushd bedtools2
make
popd


################################
# JASPAR
################################
mkdir ../../Jaspar
mkdir ../../Jaspar/raw
pushd ../../Jaspar/raw
# PWMs
wget 'http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt' --user-agent "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.22 (KHTML, like Gecko) Ubuntu Chromium/25.0.1364.160 Chrome/25.0.1364.160 Safari/537.22"
# expression of TFs
wget 'https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.N.pc.gz'
gunzip -vf 57epigenomes.N.pc.gz
# names of TFs
wget 'https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/Ensembl_v65.Gencode_v10.ENSG.gene_info'

/cygdrive/c/Python38/python.exe pwm_expressions.py

popd

################################
# ChromHMM
################################
mkdir ../../ChromHMM
pushd ../../ChromHMM
wget 'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz'
tar zxvf all.mnemonics.bedFiles.tgz
gunzip -vf E*
popd
################################
# BigWigToBedGraph (ubuntu, no windows support)
################################
wget -O bigWigToBedGraph "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph"


mkdir ../raw_bed
cd ../raw_bed

################################
# HG19 genome size
################################
wget -O hg19.chrom.sizes http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

################################
# Known Genes
################################

# failed on slop:
# wget -O hg19.KnownGenes.bed.gz hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
# gunzip -f hg19.KnownGenes.bed.gz

# failed on sort:
# wget -O hg19.KnownGenes.bed "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=597259395_zm48O4B23IAG8yLNYr5ekoweAIGe&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED"

# failed to connect:
# wget --no-check-certificate -O hg19.KnownGenes.bed.gz https://10gbps-io.dl.sourceforge.net/project/rseqc/BED/Human_Homo_sapiens/hg19_UCSC_knownGene.bed.gz

# success:
wget --no-check-certificate -O hg19.KnownGenes.bed.gz https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_UCSC_knownGene.bed.gz
gunzip -f hg19.KnownGenes.bed.gz

bedtools slop -i ./hg19.KnownGenes.bed -g ./hg19.chrom.sizes -b 5000 > ./hg19.KnownGenes.slopped.bed
ls -la

mkdir ../processed
bedtools sort -faidx hg19.chrom.sizes -i hg19.KnownGenes.slopped.bed > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > background.narrowPeak

cat ../scripts/download_and_process_one.sh | dos2unix -u > ../scripts/download_and_process_one2.sh
mv ../scripts/download_and_process_one2.sh ../scripts/download_and_process_one.sh



# 44 tissue data with all needed tracks
# 1: E003 - H1 Cells [1596]
../scripts/download_and_process_one.sh E003
# 2: E004 - H1 BMP4 Derived Mesendoderm Cultured Cells [2228]
../scripts/download_and_process_one.sh E004
# 3: E005 - H1 BMP4 Derived Trophoblast Cultured Cells [1640]
../scripts/download_and_process_one.sh E005
# 4: E006 - H1 Derived Mesenchymal Stem Cells [217]
../scripts/download_and_process_one.sh E006
# 5: E007 - H1 Derived Neuronal Progenitor Cultured Cells [1005]
../scripts/download_and_process_one.sh E007
# 6: E008 - H9 Cells [917]
../scripts/download_and_process_one.sh E008
# 7: E017 - IMR90 fetal lung fibroblasts Cell Line [710]
../scripts/download_and_process_one.sh E017
# 8: E021 - iPS DF 6.9 Cells [400]
../scripts/download_and_process_one.sh E021
# 9: E022 - iPS DF 19.11 Cells [797]
../scripts/download_and_process_one.sh E022
# 10: E029 - Primary monocytes from peripheral blood [705]
../scripts/download_and_process_one.sh E029
# 11: E032 - Primary B cells from peripheral blood [635]
../scripts/download_and_process_one.sh E032
# 12: E034 - Primary T cells from peripheral blood [383]
../scripts/download_and_process_one.sh E034
# 13: E046 - Primary Natural Killer cells from peripheral blood [388]
../scripts/download_and_process_one.sh E046
# 14: E050 - Primary hematopoietic stem cells G-CSF-mobilized Female [826]
../scripts/download_and_process_one.sh E050
# 15: E055 - Foreskin Fibroblast Primary Cells skin01 [883]
../scripts/download_and_process_one.sh E055
# 16: E056 - Foreskin Fibroblast Primary Cells skin02 [347]
../scripts/download_and_process_one.sh E056
# 17: E059 - Foreskin Melanocyte Primary Cells skin01 [398]
../scripts/download_and_process_one.sh E059
# 18: E080 - Fetal Adrenal Gland [3182]
../scripts/download_and_process_one.sh E080
# 19: E084 - Fetal Intestine Large [2186]
../scripts/download_and_process_one.sh E084
# 20: E085 - Fetal Intestine Small [1556]
../scripts/download_and_process_one.sh E085
# 21: E089 - Fetal Muscle Trunk [526]
../scripts/download_and_process_one.sh E089
# 22: E090 - Fetal Muscle Leg [1901]
../scripts/download_and_process_one.sh E090
# 23: E091 - Placenta [1594]
../scripts/download_and_process_one.sh E091
# 24: E092 - Fetal Stomach [874]
../scripts/download_and_process_one.sh E092
# 25: E093 - Fetal Thymus [1076]
../scripts/download_and_process_one.sh E093
# 26: E094 - Gastric [189]
../scripts/download_and_process_one.sh E094
# 27: E097 - Ovary [855]
../scripts/download_and_process_one.sh E097
# 28: E098 - Pancreas [384]
../scripts/download_and_process_one.sh E098
# 29: E100 - Psoas Muscle [181]
../scripts/download_and_process_one.sh E100
# 30: E109 - Small Intestine [218]
../scripts/download_and_process_one.sh E109
# 31: E114 - A549 EtOH 0.02pct Lung Carcinoma Cell Line [481]
../scripts/download_and_process_one.sh E114
# 32: E116 - GM12878 Lymphoblastoid Cells [638]
../scripts/download_and_process_one.sh E116
# 33: E117 - HeLa-S3 Cervical Carcinoma Cell Line [703]
../scripts/download_and_process_one.sh E117
# 34: E118 - HepG2 Hepatocellular Carcinoma Cell Line [478]
../scripts/download_and_process_one.sh E118
# 35: E119 - HMEC Mammary Epithelial Primary Cells [1626]
../scripts/download_and_process_one.sh E119
# 36: E120 - HSMM Skeletal Muscle Myoblasts Cells [875]
../scripts/download_and_process_one.sh E120
# 37: E121 - HSMM cell derived Skeletal Muscle Myotubes Cells [2035]
../scripts/download_and_process_one.sh E121
# 38: E122 - HUVEC Umbilical Vein Endothelial Primary Cells [1420]
../scripts/download_and_process_one.sh E122
# 39: E123 - K562 Leukemia Cells [643]
../scripts/download_and_process_one.sh E123
# 40: E124 - Monocytes-CD14+ RO01746 Primary Cells [871]
../scripts/download_and_process_one.sh E124
# 41: E125 - NH-A Astrocytes Primary Cells [1325]
../scripts/download_and_process_one.sh E125
# 42: E126 - NHDF-Ad Adult Dermal Fibroblast Primary Cells [309]
../scripts/download_and_process_one.sh E126
# 43: E127 - NHEK-Epidermal Keratinocyte Primary Cells [1207]
../scripts/download_and_process_one.sh E127
# 44: E128 - NHLF Lung Fibroblast Primary Cells [92]
../scripts/download_and_process_one.sh E128

# without full data
# ../scripts/download_and_process_one.sh E001
# ../scripts/download_and_process_one.sh E002
# ../scripts/download_and_process_one.sh E003
# ../scripts/download_and_process_one.sh E004
# ../scripts/download_and_process_one.sh E005
# ../scripts/download_and_process_one.sh E006
# ../scripts/download_and_process_one.sh E007
# ../scripts/download_and_process_one.sh E008
# ../scripts/download_and_process_one.sh E009
# ../scripts/download_and_process_one.sh E010
# ../scripts/download_and_process_one.sh E011
# ../scripts/download_and_process_one.sh E012
# ../scripts/download_and_process_one.sh E013
# ../scripts/download_and_process_one.sh E014
# ../scripts/download_and_process_one.sh E015
# ../scripts/download_and_process_one.sh E016
# ../scripts/download_and_process_one.sh E017
# ../scripts/download_and_process_one.sh E018
# ../scripts/download_and_process_one.sh E019
# ../scripts/download_and_process_one.sh E021
# ../scripts/download_and_process_one.sh E022
# ../scripts/download_and_process_one.sh E023
# ../scripts/download_and_process_one.sh E024
# ../scripts/download_and_process_one.sh E025
# ../scripts/download_and_process_one.sh E026
# ../scripts/download_and_process_one.sh E027
# ../scripts/download_and_process_one.sh E028
# ../scripts/download_and_process_one.sh E029
# ../scripts/download_and_process_one.sh E030
# ../scripts/download_and_process_one.sh E031
# ../scripts/download_and_process_one.sh E032
# ../scripts/download_and_process_one.sh E033
# ../scripts/download_and_process_one.sh E034
# ../scripts/download_and_process_one.sh E035
# ../scripts/download_and_process_one.sh E036
# ../scripts/download_and_process_one.sh E037
# ../scripts/download_and_process_one.sh E038
# ../scripts/download_and_process_one.sh E039
# ../scripts/download_and_process_one.sh E040
# ../scripts/download_and_process_one.sh E041
# ../scripts/download_and_process_one.sh E042
# ../scripts/download_and_process_one.sh E043
# ../scripts/download_and_process_one.sh E044
# ../scripts/download_and_process_one.sh E045
# ../scripts/download_and_process_one.sh E046
# ../scripts/download_and_process_one.sh E047
# ../scripts/download_and_process_one.sh E048
# ../scripts/download_and_process_one.sh E049
# ../scripts/download_and_process_one.sh E050
# ../scripts/download_and_process_one.sh E051
# ../scripts/download_and_process_one.sh E052
# ../scripts/download_and_process_one.sh E053
# ../scripts/download_and_process_one.sh E054
# ../scripts/download_and_process_one.sh E055
# ../scripts/download_and_process_one.sh E056
# ../scripts/download_and_process_one.sh E057
# ../scripts/download_and_process_one.sh E058
# ../scripts/download_and_process_one.sh E059
# ../scripts/download_and_process_one.sh E060
# ../scripts/download_and_process_one.sh E061
# ../scripts/download_and_process_one.sh E062
# ../scripts/download_and_process_one.sh E063
# ../scripts/download_and_process_one.sh E064
# ../scripts/download_and_process_one.sh E065
# ../scripts/download_and_process_one.sh E066
# ../scripts/download_and_process_one.sh E067
# ../scripts/download_and_process_one.sh E068
# ../scripts/download_and_process_one.sh E069
# ../scripts/download_and_process_one.sh E070
# ../scripts/download_and_process_one.sh E071
# ../scripts/download_and_process_one.sh E072
# ../scripts/download_and_process_one.sh E073
# ../scripts/download_and_process_one.sh E074
# ../scripts/download_and_process_one.sh E075
# ../scripts/download_and_process_one.sh E076
# ../scripts/download_and_process_one.sh E077
# ../scripts/download_and_process_one.sh E078
# ../scripts/download_and_process_one.sh E079
# ../scripts/download_and_process_one.sh E080
# ../scripts/download_and_process_one.sh E081
# ../scripts/download_and_process_one.sh E082
# ../scripts/download_and_process_one.sh E083
# ../scripts/download_and_process_one.sh E084
# ../scripts/download_and_process_one.sh E085
# ../scripts/download_and_process_one.sh E086
# ../scripts/download_and_process_one.sh E087
# ../scripts/download_and_process_one.sh E088
# ../scripts/download_and_process_one.sh E089
# ../scripts/download_and_process_one.sh E090
# ../scripts/download_and_process_one.sh E091
# ../scripts/download_and_process_one.sh E092
# ../scripts/download_and_process_one.sh E093
# ../scripts/download_and_process_one.sh E094
# ../scripts/download_and_process_one.sh E095
# ../scripts/download_and_process_one.sh E096
# ../scripts/download_and_process_one.sh E097
# ../scripts/download_and_process_one.sh E098
# ../scripts/download_and_process_one.sh E099
# ../scripts/download_and_process_one.sh E100
# ../scripts/download_and_process_one.sh E101
# ../scripts/download_and_process_one.sh E102
# ../scripts/download_and_process_one.sh E103
# ../scripts/download_and_process_one.sh E104
# ../scripts/download_and_process_one.sh E105
# ../scripts/download_and_process_one.sh E106
# ../scripts/download_and_process_one.sh E107
# ../scripts/download_and_process_one.sh E108
# ../scripts/download_and_process_one.sh E109
# ../scripts/download_and_process_one.sh E110
# ../scripts/download_and_process_one.sh E111
# ../scripts/download_and_process_one.sh E112
# ../scripts/download_and_process_one.sh E113
# ../scripts/download_and_process_one.sh E114
# ../scripts/download_and_process_one.sh E115
# ../scripts/download_and_process_one.sh E116
# ../scripts/download_and_process_one.sh E117
# ../scripts/download_and_process_one.sh E118
# ../scripts/download_and_process_one.sh E119
# ../scripts/download_and_process_one.sh E120
# ../scripts/download_and_process_one.sh E121
# ../scripts/download_and_process_one.sh E122
# ../scripts/download_and_process_one.sh E123
# ../scripts/download_and_process_one.sh E124
# ../scripts/download_and_process_one.sh E125
# ../scripts/download_and_process_one.sh E126
# ../scripts/download_and_process_one.sh E127
# ../scripts/download_and_process_one.sh E128


cp -f background.narrowPeak ../processed/background.cleaned.narrowPeak
mv -f hg19.KnownGenes.bed ../processed/genes.cleaned.narrowPeak
cd ../scripts
python filter_background.py ../processed/background.cleaned.narrowPeak 5000

################################
# Genome FASTA
################################

mkdir ../raw_genome
cd ../raw_genome

wget -O chr1.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr1.fa.gz
wget -O chr2.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr2.fa.gz
wget -O chr3.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr3.fa.gz
wget -O chr4.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr4.fa.gz
wget -O chr5.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr5.fa.gz
wget -O chr6.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr6.fa.gz
wget -O chr7.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr7.fa.gz
wget -O chr8.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr8.fa.gz
wget -O chr9.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr9.fa.gz
wget -O chr10.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr10.fa.gz
wget -O chr11.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr11.fa.gz
wget -O chr12.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr12.fa.gz
wget -O chr13.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr13.fa.gz
wget -O chr14.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr14.fa.gz
wget -O chr15.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr15.fa.gz
wget -O chr16.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr16.fa.gz
wget -O chr17.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr17.fa.gz
wget -O chr18.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr18.fa.gz
wget -O chr19.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr19.fa.gz
wget -O chr20.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr20.fa.gz
wget -O chr21.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr21.fa.gz
wget -O chr22.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr22.fa.gz
wget -O chrX.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chrX.fa.gz
wget -O chrY.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chrY.fa.gz

gunzip -vf *.gz
rm -f *.gz


################################
# DNase + H3K27ac BigWigs
################################

mkdir ../raw_bigwig
cd ../raw_bigwig
wget -O E003-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E003-DNase.imputed.pval.signal.bigwig"
wget -O E004-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E004-DNase.imputed.pval.signal.bigwig"
wget -O E005-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E005-DNase.imputed.pval.signal.bigwig"
wget -O E006-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E006-DNase.imputed.pval.signal.bigwig"
wget -O E007-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E007-DNase.imputed.pval.signal.bigwig"
wget -O E008-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E008-DNase.imputed.pval.signal.bigwig"
wget -O E017-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E017-DNase.imputed.pval.signal.bigwig"
wget -O E021-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E021-DNase.imputed.pval.signal.bigwig"
wget -O E022-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E022-DNase.imputed.pval.signal.bigwig"
wget -O E029-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E029-DNase.imputed.pval.signal.bigwig"
wget -O E032-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E032-DNase.imputed.pval.signal.bigwig"
wget -O E034-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E034-DNase.imputed.pval.signal.bigwig"
wget -O E046-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E046-DNase.imputed.pval.signal.bigwig"
wget -O E050-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E050-DNase.imputed.pval.signal.bigwig"
wget -O E055-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E055-DNase.imputed.pval.signal.bigwig"
wget -O E056-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E056-DNase.imputed.pval.signal.bigwig"
wget -O E059-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E059-DNase.imputed.pval.signal.bigwig"
wget -O E080-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E080-DNase.imputed.pval.signal.bigwig"
wget -O E084-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E084-DNase.imputed.pval.signal.bigwig"
wget -O E085-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E085-DNase.imputed.pval.signal.bigwig"
wget -O E089-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E089-DNase.imputed.pval.signal.bigwig"
wget -O E090-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E090-DNase.imputed.pval.signal.bigwig"
wget -O E091-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E091-DNase.imputed.pval.signal.bigwig"
wget -O E092-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E092-DNase.imputed.pval.signal.bigwig"
wget -O E093-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E093-DNase.imputed.pval.signal.bigwig"
wget -O E094-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E094-DNase.imputed.pval.signal.bigwig"
wget -O E097-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E097-DNase.imputed.pval.signal.bigwig"
wget -O E098-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E098-DNase.imputed.pval.signal.bigwig"
wget -O E100-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E100-DNase.imputed.pval.signal.bigwig"
wget -O E109-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E109-DNase.imputed.pval.signal.bigwig"
wget -O E114-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E114-DNase.imputed.pval.signal.bigwig"
wget -O E116-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E116-DNase.imputed.pval.signal.bigwig"
wget -O E117-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E117-DNase.imputed.pval.signal.bigwig"
wget -O E118-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E118-DNase.imputed.pval.signal.bigwig"
wget -O E119-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E119-DNase.imputed.pval.signal.bigwig"
wget -O E120-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E120-DNase.imputed.pval.signal.bigwig"
wget -O E121-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E121-DNase.imputed.pval.signal.bigwig"
wget -O E122-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E122-DNase.imputed.pval.signal.bigwig"
wget -O E123-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E123-DNase.imputed.pval.signal.bigwig"
wget -O E124-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E124-DNase.imputed.pval.signal.bigwig"
wget -O E125-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E125-DNase.imputed.pval.signal.bigwig"
wget -O E126-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E126-DNase.imputed.pval.signal.bigwig"
wget -O E127-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E127-DNase.imputed.pval.signal.bigwig"
wget -O E128-DNase.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/DNase/E128-DNase.imputed.pval.signal.bigwig"

wget -O E003-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E003-H3K27ac.imputed.pval.signal.bigwig"
wget -O E004-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E004-H3K27ac.imputed.pval.signal.bigwig"
wget -O E005-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E005-H3K27ac.imputed.pval.signal.bigwig"
wget -O E006-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E006-H3K27ac.imputed.pval.signal.bigwig"
wget -O E007-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E007-H3K27ac.imputed.pval.signal.bigwig"
wget -O E008-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E008-H3K27ac.imputed.pval.signal.bigwig"
wget -O E017-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E017-H3K27ac.imputed.pval.signal.bigwig"
wget -O E021-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E021-H3K27ac.imputed.pval.signal.bigwig"
wget -O E022-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E022-H3K27ac.imputed.pval.signal.bigwig"
wget -O E029-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E029-H3K27ac.imputed.pval.signal.bigwig"
wget -O E032-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E032-H3K27ac.imputed.pval.signal.bigwig"
wget -O E034-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E034-H3K27ac.imputed.pval.signal.bigwig"
wget -O E046-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E046-H3K27ac.imputed.pval.signal.bigwig"
wget -O E050-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E050-H3K27ac.imputed.pval.signal.bigwig"
wget -O E055-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E055-H3K27ac.imputed.pval.signal.bigwig"
wget -O E056-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E056-H3K27ac.imputed.pval.signal.bigwig"
wget -O E059-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E059-H3K27ac.imputed.pval.signal.bigwig"
wget -O E080-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E080-H3K27ac.imputed.pval.signal.bigwig"
wget -O E084-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E084-H3K27ac.imputed.pval.signal.bigwig"
wget -O E085-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E085-H3K27ac.imputed.pval.signal.bigwig"
wget -O E089-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E089-H3K27ac.imputed.pval.signal.bigwig"
wget -O E090-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFilfeType/signal/consolidatedImputed/H3K27ac/E090-H3K27ac.imputed.pval.signal.bigwig"
wget -O E091-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E091-H3K27ac.imputed.pval.signal.bigwig"
wget -O E092-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E092-H3K27ac.imputed.pval.signal.bigwig"
wget -O E093-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E093-H3K27ac.imputed.pval.signal.bigwig"
wget -O E094-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E094-H3K27ac.imputed.pval.signal.bigwig"
wget -O E097-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E097-H3K27ac.imputed.pval.signal.bigwig"
wget -O E098-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E098-H3K27ac.imputed.pval.signal.bigwig"
wget -O E100-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E100-H3K27ac.imputed.pval.signal.bigwig"
wget -O E109-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E109-H3K27ac.imputed.pval.signal.bigwig"
wget -O E114-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E114-H3K27ac.imputed.pval.signal.bigwig"
wget -O E116-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E116-H3K27ac.imputed.pval.signal.bigwig"
wget -O E117-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E117-H3K27ac.imputed.pval.signal.bigwig"
wget -O E118-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E118-H3K27ac.imputed.pval.signal.bigwig"
wget -O E119-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E119-H3K27ac.imputed.pval.signal.bigwig"
wget -O E120-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E120-H3K27ac.imputed.pval.signal.bigwig"
wget -O E121-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E121-H3K27ac.imputed.pval.signal.bigwig"
wget -O E122-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E122-H3K27ac.imputed.pval.signal.bigwig"
wget -O E123-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E123-H3K27ac.imputed.pval.signal.bigwig"
wget -O E124-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E124-H3K27ac.imputed.pval.signal.bigwig"
wget -O E125-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E125-H3K27ac.imputed.pval.signal.bigwig"
wget -O E126-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E126-H3K27ac.imputed.pval.signal.bigwig"
wget -O E127-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E127-H3K27ac.imputed.pval.signal.bigwig"
wget -O E128-H3K27ac.pval.signal.bigwig "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E128-H3K27ac.imputed.pval.signal.bigwig"


# bedgraph of ac and dnase on enhancer sequences
mkdir ../raw_bedgraphs
mkdir ../processed_bedgraphs
mkdir ../tmp

rm ../tmp/enhancers.bed
cat ../processed/E003-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E004-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E005-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E006-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E007-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E008-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E017-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E021-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E022-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E029-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E032-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E034-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E046-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E050-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E055-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E056-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E059-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E080-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E084-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E085-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E089-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E090-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E091-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E092-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E093-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E094-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E097-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E098-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E100-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E109-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E114-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E116-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E117-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E118-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E119-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E120-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E121-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E122-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E123-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E124-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E125-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E126-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E127-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/E128-H3K27ac.cleaned.narrowPeak >> ../tmp/enhancers.bed
cat ../processed/background.cleaned.narrowPeak >> ../tmp/enhancers.bed
bedtools sort -i ../tmp/enhancers.bed > ../tmp/enhancers.sorted.bed
bedtools merge -i ../tmp/enhancers.sorted.bed > ../tmp/enhancers.merged.bed
bedtools slop -i ../tmp/enhancers.merged.bed -g ../raw_bed/hg19.chrom.sizes -b 3000 > ../tmp/enhancers.slopped.bed


################################
# BigWigToBedGraph (on linux only)
################################

export PROMPT_COMMAND="echo -n \[\$(date +%H:%M:%S)\]\ "

../scripts/bigWigToBedGraph ../raw_bigwig/E003-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E003-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E004-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E004-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E005-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E005-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E006-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E006-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E007-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E007-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E008-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E008-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E017-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E017-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E021-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E021-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E022-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E022-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E029-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E029-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E032-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E032-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E034-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E034-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E046-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E046-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E050-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E050-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E055-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E055-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E056-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E056-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E059-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E059-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E080-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E080-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E084-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E084-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E085-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E085-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E089-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E089-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E090-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E090-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E091-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E091-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E092-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E092-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E093-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E093-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E094-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E094-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E097-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E097-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E098-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E098-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E100-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E100-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E109-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E109-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E114-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E114-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E116-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E116-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E117-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E117-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E118-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E118-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E119-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E119-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E120-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E120-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E121-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E121-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E122-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E122-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E123-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E123-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E124-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E124-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E125-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E125-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E126-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E126-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E127-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E127-H3K27ac.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E128-H3K27ac.pval.signal.bigwig ../raw_bedgraphs/E128-H3K27ac.pval.bedgraph

../scripts/bigWigToBedGraph ../raw_bigwig/E003-DNase.pval.signal.bigwig ../raw_bedgraphs/E003-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E004-DNase.pval.signal.bigwig ../raw_bedgraphs/E004-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E005-DNase.pval.signal.bigwig ../raw_bedgraphs/E005-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E006-DNase.pval.signal.bigwig ../raw_bedgraphs/E006-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E007-DNase.pval.signal.bigwig ../raw_bedgraphs/E007-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E008-DNase.pval.signal.bigwig ../raw_bedgraphs/E008-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E017-DNase.pval.signal.bigwig ../raw_bedgraphs/E017-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E021-DNase.pval.signal.bigwig ../raw_bedgraphs/E021-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E022-DNase.pval.signal.bigwig ../raw_bedgraphs/E022-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E029-DNase.pval.signal.bigwig ../raw_bedgraphs/E029-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E032-DNase.pval.signal.bigwig ../raw_bedgraphs/E032-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E034-DNase.pval.signal.bigwig ../raw_bedgraphs/E034-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E046-DNase.pval.signal.bigwig ../raw_bedgraphs/E046-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E050-DNase.pval.signal.bigwig ../raw_bedgraphs/E050-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E055-DNase.pval.signal.bigwig ../raw_bedgraphs/E055-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E056-DNase.pval.signal.bigwig ../raw_bedgraphs/E056-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E059-DNase.pval.signal.bigwig ../raw_bedgraphs/E059-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E080-DNase.pval.signal.bigwig ../raw_bedgraphs/E080-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E084-DNase.pval.signal.bigwig ../raw_bedgraphs/E084-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E085-DNase.pval.signal.bigwig ../raw_bedgraphs/E085-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E089-DNase.pval.signal.bigwig ../raw_bedgraphs/E089-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E090-DNase.pval.signal.bigwig ../raw_bedgraphs/E090-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E091-DNase.pval.signal.bigwig ../raw_bedgraphs/E091-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E092-DNase.pval.signal.bigwig ../raw_bedgraphs/E092-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E093-DNase.pval.signal.bigwig ../raw_bedgraphs/E093-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E094-DNase.pval.signal.bigwig ../raw_bedgraphs/E094-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E097-DNase.pval.signal.bigwig ../raw_bedgraphs/E097-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E098-DNase.pval.signal.bigwig ../raw_bedgraphs/E098-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E100-DNase.pval.signal.bigwig ../raw_bedgraphs/E100-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E109-DNase.pval.signal.bigwig ../raw_bedgraphs/E109-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E114-DNase.pval.signal.bigwig ../raw_bedgraphs/E114-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E116-DNase.pval.signal.bigwig ../raw_bedgraphs/E116-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E117-DNase.pval.signal.bigwig ../raw_bedgraphs/E117-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E118-DNase.pval.signal.bigwig ../raw_bedgraphs/E118-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E119-DNase.pval.signal.bigwig ../raw_bedgraphs/E119-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E120-DNase.pval.signal.bigwig ../raw_bedgraphs/E120-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E121-DNase.pval.signal.bigwig ../raw_bedgraphs/E121-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E122-DNase.pval.signal.bigwig ../raw_bedgraphs/E122-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E123-DNase.pval.signal.bigwig ../raw_bedgraphs/E123-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E124-DNase.pval.signal.bigwig ../raw_bedgraphs/E124-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E125-DNase.pval.signal.bigwig ../raw_bedgraphs/E125-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E126-DNase.pval.signal.bigwig ../raw_bedgraphs/E126-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E127-DNase.pval.signal.bigwig ../raw_bedgraphs/E127-DNase.pval.bedgraph
../scripts/bigWigToBedGraph ../raw_bigwig/E128-DNase.pval.signal.bigwig ../raw_bedgraphs/E128-DNase.pval.bedgraph

bedtools intersect -a ../raw_bedgraphs/E003-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E003-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E004-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E004-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E005-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E005-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E006-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E006-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E007-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E007-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E008-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E008-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E017-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E017-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E021-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E021-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E022-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E022-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E029-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E029-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E032-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E032-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E034-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E034-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E046-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E046-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E050-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E050-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E055-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E055-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E056-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E056-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E059-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E059-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E080-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E080-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E084-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E084-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E085-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E085-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E089-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E089-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E090-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E090-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E091-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E091-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E092-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E092-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E093-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E093-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E094-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E094-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E097-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E097-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E098-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E098-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E100-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E100-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E109-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E109-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E114-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E114-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E116-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E116-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E117-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E117-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E118-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E118-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E119-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E119-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E120-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E120-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E121-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E121-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E122-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E122-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E123-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E123-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E124-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E124-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E125-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E125-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E126-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E126-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E127-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E127-H3K27ac.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E128-H3K27ac.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E128-H3K27ac.enh.bedgraph
rm ../raw_bedgraphs/*H3K27ac.pval.bedgraph


bedtools intersect -a ../raw_bedgraphs/E003-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E003-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E004-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E004-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E005-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E005-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E006-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E006-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E007-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E007-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E008-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E008-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E017-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E017-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E021-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E021-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E022-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E022-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E029-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E029-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E032-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E032-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E034-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E034-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E046-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E046-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E050-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E050-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E055-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E055-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E056-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E056-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E059-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E059-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E080-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E080-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E084-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E084-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E085-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E085-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E089-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E089-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E090-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E090-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E091-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E091-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E092-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E092-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E093-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E093-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E094-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E094-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E097-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E097-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E098-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E098-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E100-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E100-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E109-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E109-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E114-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E114-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E116-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E116-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E117-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E117-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E118-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E118-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E119-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E119-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E120-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E120-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E121-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E121-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E122-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E122-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E123-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E123-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E124-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E124-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E125-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E125-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E126-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E126-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E127-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E127-DNase.enh.bedgraph
bedtools intersect -a ../raw_bedgraphs/E128-DNase.pval.bedgraph -b ../tmp/enhancers.slopped.bed > ../processed_bedgraphs/E128-DNase.enh.bedgraph

rm ../raw_bedgraphs/*DNase.pval.bedgraph


cat ../processed_bedgraphs/E003-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E003-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E004-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E004-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E005-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E005-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E006-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E006-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E007-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E007-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E008-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E008-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E017-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E017-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E021-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E021-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E022-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E022-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E029-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E029-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E032-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E032-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E034-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E034-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E046-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E046-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E050-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E050-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E055-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E055-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E056-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E056-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E059-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E059-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E080-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E080-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E084-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E084-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E085-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E085-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E089-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E089-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E090-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E090-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E091-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E091-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E092-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E092-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E093-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E093-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E094-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E094-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E097-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E097-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E098-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E098-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E100-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E100-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E109-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E109-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E114-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E114-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E116-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E116-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E117-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E117-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E118-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E118-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E119-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E119-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E120-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E120-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E121-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E121-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E122-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E122-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E123-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E123-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E124-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E124-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E125-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E125-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E126-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E126-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E127-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E127-H3K27ac.enh.u.bedgraph
cat ../processed_bedgraphs/E128-H3K27ac.enh.bedgraph | uniq > ../processed_bedgraphs/E128-H3K27ac.enh.u.bedgraph

rm ../processed_bedgraphs/*-H3K27ac.enh.bedgraph

cat ../processed_bedgraphs/E003-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E003-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E004-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E004-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E005-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E005-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E006-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E006-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E007-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E007-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E008-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E008-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E017-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E017-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E021-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E021-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E022-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E022-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E029-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E029-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E032-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E032-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E034-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E034-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E046-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E046-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E050-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E050-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E055-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E055-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E056-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E056-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E059-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E059-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E080-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E080-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E084-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E084-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E085-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E085-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E089-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E089-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E090-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E090-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E091-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E091-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E092-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E092-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E093-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E093-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E094-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E094-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E097-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E097-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E098-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E098-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E100-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E100-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E109-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E109-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E114-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E114-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E116-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E116-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E117-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E117-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E118-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E118-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E119-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E119-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E120-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E120-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E121-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E121-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E122-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E122-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E123-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E123-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E124-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E124-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E125-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E125-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E126-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E126-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E127-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E127-DNase.enh.u.bedgraph
cat ../processed_bedgraphs/E128-DNase.enh.bedgraph | uniq > ../processed_bedgraphs/E128-DNase.enh.u.bedgraph

rm ../processed_bedgraphs/*-DNase.enh.bedgraph
