# Downloaded from roadmap site (2013):
# from page http://egg2.wustl.edu/roadmap/web_portal/processed_data.html
# http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
# downloaded E109 small intentine and E034 T cells
# H3K27ac, H3K27me3, H3K4me3, H3K4me1


# HG19 known Genes (2009):
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=597259395_zm48O4B23IAG8yLNYr5ekoweAIGe&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED



cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/raw/roadmap
#HG19 genome size:
# wget -O hg19.chrom.sizes http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

#HG19 known Genes (2009):
# wget -O hg9.KnownGenes.bed "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=597259395_zm48O4B23IAG8yLNYr5ekoweAIGe&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED"
# bedtools slop -i ../hg19.KnownGenes.bed -g ../hg19.chrom.sizes.txt -b 5000 > ../hg19.KnownGenes.slopped.bed
ls -la ..


mkdir bed
mkdir processed

/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E003
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E004
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E005
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E006
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E007
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E008
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E017
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E021
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E022
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E029
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E032
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E034
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E046
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E050
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E055
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E056
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E059
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E080
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E084
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E085
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E089
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E090
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E091
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E092
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E093
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E094
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E097
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E098
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E100
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E109
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E114
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E116
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E117
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E118
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E119
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E120
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E121
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E122
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E123
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E124
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E125
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E126
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E127
/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks/download_and_process_one.sh E128

#then in matlab

# cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
# peaks.beds2mats()
# peaks.mergePeakFiles()


# load('data/peaks/raw/roadmap/merged/mergedPeaks.mat');
# peaks.minimizeMergePeak(mergedPeaks);
# mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
# %mergedPeaksMin = load('data/peaks/randGenerated/mergedPeaksMinimized.mat');
# mainPWMCor(mergedPeaksMin);

