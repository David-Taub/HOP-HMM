
cd ../raw

ls -la
wget -O $1-H3K27ac.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27ac.narrowPeak.gz"
wget -O $1-H3K4me3.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me3.narrowPeak.gz"
wget -O $1-H3K4me1.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me1.narrowPeak.gz"
wget -O $1-H3K27me3.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27me3.narrowPeak.gz"

gunzip -f *.gz
 rm -f *.gz

#take only H3k27ac peaks with H3k4me1 peaks
bedtools intersect -a $1-H3K27ac.narrowPeak -b $1-H3K4me1.narrowPeak -wa >$1-H3K27ac.cleaned.narrowPeak
ls -la

#remove duplicates
awk '!x[$0]++' $1-H3K27ac.cleaned.narrowPeak > $1-H3K27ac.tmp.narrowPeak
mv -f $1-H3K27ac.tmp.narrowPeak $1-H3K27ac.cleaned.narrowPeak
ls -la


#remove H3k27me3 peaks
bedtools subtract -a $1-H3K27ac.cleaned.narrowPeak -b $1-H3K27me3.narrowPeak -A > $1-H3K27ac.tmp.narrowPeak
mv -f $1-H3K27ac.tmp.narrowPeak $1-H3K27ac.cleaned.narrowPeak
ls -la

#remove H3k4me3 peaks
bedtools subtract -a $1-H3K27ac.cleaned.narrowPeak -b $1-H3K4me3.narrowPeak -A > $1-H3K27ac.tmp.narrowPeak
mv -f $1-H3K27ac.tmp.narrowPeak $1-H3K27ac.cleaned.narrowPeak
ls -la

#remove padded known genes
bedtools subtract -a $1-H3K27ac.cleaned.narrowPeak -b ./hg19.KnownGenes.slopped.bed -A > $1-H3K27ac.tmp.narrowPeak
mv -f $1-H3K27ac.tmp.narrowPeak $1-H3K27ac.cleaned.narrowPeak
ls -la
mv -f $1-H3K27ac.cleaned.narrowPeak ../processed/$1-H3K27ac.cleaned.narrowPeak

