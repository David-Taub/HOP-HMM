
cd ../raw_bed

ls -la
wget -O $1-H3K27ac.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27ac.narrowPeak.gz"
wget -O $1-H3K4me3.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me3.narrowPeak.gz"
wget -O $1-H3K4me1.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me1.narrowPeak.gz"
wget -O $1-H3K27me3.narrowPeak.gz "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27me3.narrowPeak.gz"

gunzip -f ./*.gz
 rm -f ./*.gz

ls -la
# refine with background
bedtools sort -faidx hg19.chrom.sizes -i $1-H3K27me3.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
bedtools sort -faidx hg19.chrom.sizes -i $1-H3K27ac.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
bedtools sort -faidx hg19.chrom.sizes -i $1-H3K4me1.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
bedtools sort -faidx hg19.chrom.sizes -i $1-H3K4me3.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
rm -f tmp.*

#take only H3k27ac peaks with H3k4me1 peaks
bedtools intersect -a $1-H3K27ac.narrowPeak -b $1-H3K4me1.narrowPeak -wa >$1-H3K27ac.cleaned.narrowPeak
ls -la

# sort -V -k 1,3 "total_merged.narrowPeak" | sortBed | tee total_merged.tmp.narrowPeak | sort -c -k1,1 -k2,2n || true
# # bedtools sort -i total_merged.narrowPeak > total_merged.tmp.narrowPeak
# bedtools merge -i total_merged.tmp.narrowPeak -d 10000 > total_merged.narrowPeak
# rm -f total_merged.tmp.narrowPeak

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


