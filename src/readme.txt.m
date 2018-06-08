
% first we choose the PWMs that are most different between the following tissues,
% including background and genes:
cd /cygdrive/c/users/booga/Dropbox/projects/HOP-HMM/src
cd ~/projects/HopHMM/src
k = 40
peaks.beds2mats(500)
peaks.mergePeakFiles(true, true)
load('../data/peaks/mergedPeaks.mat');
peaks.minimizeMergePeak(mergedPeaks, 500, tissueNames)
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat');
% delete('../data/precomputation/SelectedPWMs.mat');
% mainGenSequences(10000, 500, 5, 1000, false);
% mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized_fake.mat');
% chooseBestPWMs(mergedPeaksMin, [1,2,3,5], 1000)
JasparDataProcessing.mainPreprocessPWMs();
chooseBestPWMs(mergedPeaksMin, [10, 20, 30, backgroundIndex], 1000)

geneIndex = find(strcmp('genes', tissueNames))
backgroundIndex = find(strcmp('background', tissueNames))

pretrain(mergedPeaksMin, k, [10, 20, 30], backgroundIndex)



% then we find multiEnhancers, and try to learn them
peaks.beds2matsNoSeq()
peaks.mergePeakFiles(false, false)
load('../data/peaks/mergedPeaksNoBackground.mat');
multiEnhancers = peaks.multiEnhancerCaller(mergedPeaks, 10000, [10, 20, 30]);
load('../data/multiEnhancers.mat');
k = 40
mainRealData(multiEnhancers, 4, k, false, true);

mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat');


mergedPeaksMin =
        lengths: [121279×1 int32]
       overlaps: [121279×46 double]
           seqs: [121279×500 double]
    tissueNames: {1×46 cell}
