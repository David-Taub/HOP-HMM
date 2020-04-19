cd ../raw_bed

###############################
# Downloads
###############################

wget -O $1-H3K27ac.narrowPeak.gz "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27ac.narrowPeak.gz"
if [ ! -f $1-H3K27ac.narrowPeak.gz ]; then
    echo "$1-H3K27ac.narrowPeak.gz not found!"
    exit 1
fi
wget -O $1-DNase.macs2.narrowPeak.gz "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-DNase.macs2.narrowPeak.gz"
if [ ! -f $1-DNase.macs2.narrowPeak.gz ]; then
    echo "$1-DNase.macs2.narrowPeak.gz not found!"
    exit 1
fi
wget -O $1-H3K4me3.narrowPeak.gz "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me3.narrowPeak.gz"
if [ ! -f $1-H3K4me3.narrowPeak.gz ]; then
    echo "$1-H3K4me3.narrowPeak.gz not found!"
    exit 1
fi
wget -O $1-H3K4me1.narrowPeak.gz "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K4me1.narrowPeak.gz"
if [ ! -f $1-H3K4me1.narrowPeak.gz ]; then
    echo "$1-H3K4me1.narrowPeak.gz not found!"
    exit 1
fi
wget -O $1-H3K27me3.narrowPeak.gz "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$1-H3K27me3.narrowPeak.gz"
if [ ! -f $1-H3K27me3.narrowPeak.gz ]; then
    echo "$1-H3K27me3.narrowPeak.gz not found!"
    exit 1
fi

echo "done downloading"
gunzip -f ./*.gz
rm -f ./*.gz
echo "done extracting"

find . -empty -delete

if [ ! -f $1-H3K27ac.narrowPeak ]; then
    echo "$1-H3K27ac.narrowPeak not found!"
    exit 1
fi
if [ ! -f $1-DNase.macs2.narrowPeak ]; then
    echo "$1-H3K27me3.narrowPeak not found!"
    exit 1
fi
if [ ! -f $1-H3K4me1.narrowPeak ]; then
    echo "$1-H3K4me1.narrowPeak not found!"
    exit 1
fi
if [ ! -f $1-H3K4me3.narrowPeak ]; then
    echo "$1-H3K4me3.narrowPeak not found!"
    exit 1
fi
if [ ! -f $1-H3K27me3.narrowPeak ]; then
    echo "$1-H3K27me3.narrowPeak not found!"
    exit 1
fi
echo "Asserted files existing"

###############################
# Background
###############################
echo "Refining background"

bedtools sort -faidx hg19.chrom.sizes -i $1-H3K27ac.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
bedtools sort -faidx hg19.chrom.sizes -i $1-DNase.macs2.narrowPeak > tmp.1.narrowPeak
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
bedtools sort -faidx hg19.chrom.sizes -i $1-H3K27me3.narrowPeak > tmp.1.narrowPeak
bedtools complement  -i tmp.1.narrowPeak -g hg19.chrom.sizes > tmp.2.narrowPeak
bedtools intersect -a tmp.2.narrowPeak -b background.narrowPeak > tmp.3.narrowPeak
mv -f tmp.3.narrowPeak background.narrowPeak
rm -f tmp.*

find . -empty -delete

ls -la background.narrowPeak

if [ ! -f background.narrowPeak ]; then
    echo "background.narrowPeak not found!"
    exit 1
fi
echo "bedtools run on background completed"


###############################
# Enhancers
###############################
#take only peaks with DNASE-I peaks
bedtools intersect -a $1-DNase.macs2.narrowPeak -b $1-H3K27ac.narrowPeak -wa > tmp.1.$1-H3K27ac.cleaned.narrowPeak
#take only H3k27ac peaks with H3k4me1 peaks
bedtools intersect -a tmp.1.$1-H3K27ac.cleaned.narrowPeak -b $1-H3K4me1.narrowPeak -wa > tmp.2.$1-H3K27ac.cleaned.narrowPeak

#remove duplicates
awk '!x[$0]++' tmp.2.$1-H3K27ac.cleaned.narrowPeak > tmp.3.$1-H3K27ac.cleaned.narrowPeak

#remove H3k27me3 peaks
bedtools subtract -a tmp.3.$1-H3K27ac.cleaned.narrowPeak -b $1-H3K27me3.narrowPeak -A > tmp.4.$1-H3K27ac.cleaned.narrowPeak

#remove H3k4me3 peaks
bedtools subtract -a tmp.4.$1-H3K27ac.cleaned.narrowPeak -b $1-H3K4me3.narrowPeak -A > tmp.5.$1-H3K27ac.cleaned.narrowPeak

#remove padded known genes
bedtools subtract -a tmp.5.$1-H3K27ac.cleaned.narrowPeak -b ./hg19.KnownGenes.slopped.bed -A > tmp.6.$1-H3K27ac.cleaned.narrowPeak
ls -al tmp.*
mv -f tmp.6.$1-H3K27ac.cleaned.narrowPeak ../processed/$1-H3K27ac.cleaned.narrowPeak
rm -f tmp.*

#############################
# Verify Output
#############################
ls -al ../processed/$1-H3K27ac.cleaned.narrowPeak
echo "bedtools run completed"

if [ ! -f ../processed/$1-H3K27ac.cleaned.narrowPeak ]; then
    echo "../processed/$1-H3K27ac.cleaned.narrowPeak not found!"
    exit 1
fi
find . -empty -delete