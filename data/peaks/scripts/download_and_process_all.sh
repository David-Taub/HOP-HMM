# Downloaded from roadmap site (2013):
# from page http://egg2.wustl.edu/roadmap/web_portal/processed_data.html
# http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
# downloaded E109 small intentine and E034 T cells
# H3K27ac, H3K27me3, H3K4me3, H3K4me1


# bedtools installation
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make

# in cygwin on windows laptop:
# cd  ../../cygdrive/c/Users/booga/Dropbox/projects/CompGenetics/data/peaks/scripts

mkdir ../raw_bed
cd ../raw_bed
#HG19 genome size:
wget -O hg19.chrom.sizes http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

#HG19 known Genes (2009):
# wget -O hg19.KnownGenes.bed "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=597259395_zm48O4B23IAG8yLNYr5ekoweAIGe&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED"
wget -O hg19.KnownGenes.bed "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=652268259_VIkdFtkCyWYPpumVNPdYEZxkg5MP&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED"
bedtools slop -i ./hg19.KnownGenes.bed -g ./hg19.chrom.sizes -b 5000 > ./hg19.KnownGenes.slopped.bed
ls -la
rm ./hg19.KnownGenes.bed

mkdir ../processed

../scripts/download_and_process_one.sh E003
../scripts/download_and_process_one.sh E004
../scripts/download_and_process_one.sh E005
../scripts/download_and_process_one.sh E006
../scripts/download_and_process_one.sh E007
../scripts/download_and_process_one.sh E008
../scripts/download_and_process_one.sh E017
../scripts/download_and_process_one.sh E021
../scripts/download_and_process_one.sh E022
../scripts/download_and_process_one.sh E029
../scripts/download_and_process_one.sh E032
../scripts/download_and_process_one.sh E034
../scripts/download_and_process_one.sh E046
../scripts/download_and_process_one.sh E050
../scripts/download_and_process_one.sh E055
../scripts/download_and_process_one.sh E056
../scripts/download_and_process_one.sh E059
../scripts/download_and_process_one.sh E080
../scripts/download_and_process_one.sh E084
../scripts/download_and_process_one.sh E085
../scripts/download_and_process_one.sh E089
../scripts/download_and_process_one.sh E090
../scripts/download_and_process_one.sh E091
../scripts/download_and_process_one.sh E092
../scripts/download_and_process_one.sh E093
../scripts/download_and_process_one.sh E094
../scripts/download_and_process_one.sh E097
../scripts/download_and_process_one.sh E098
../scripts/download_and_process_one.sh E100
../scripts/download_and_process_one.sh E109
../scripts/download_and_process_one.sh E114
../scripts/download_and_process_one.sh E116
../scripts/download_and_process_one.sh E117
../scripts/download_and_process_one.sh E118
../scripts/download_and_process_one.sh E119
../scripts/download_and_process_one.sh E120
../scripts/download_and_process_one.sh E121
../scripts/download_and_process_one.sh E122
../scripts/download_and_process_one.sh E123
../scripts/download_and_process_one.sh E124
../scripts/download_and_process_one.sh E125
../scripts/download_and_process_one.sh E126
../scripts/download_and_process_one.sh E127
../scripts/download_and_process_one.sh E128

#download fasta of genome
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

gunzip -f *.gz
rm -f *.gz

#then in matlab

# cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
# peaks.beds2matsNoSeq()
# peaks.mergePeakFiles()


# load('data/peaks/raw/roadmap/merged/mergedPeaks.mat');
# peaks.minimizeMergePeak(mergedPeaks);
# mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
# %mergedPeaksMin = load('data/peaks/randGenerated/mergedPeaksMinimized.mat');
# mainPWMCor(mergedPeaksMin);
