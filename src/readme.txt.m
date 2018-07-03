
% first we choose the PWMs that are most different between the following tissues,
% including background and genes:
cd C:\Users\booga\Dropbox\projects\HOP-HMM\src
cd /cygdrive/c/users/booga/Dropbox/projects/HOP-HMM/src
cd ~/projects/HopHMM/src
cd ~/projects/HOP-HMM/src

% TEST

k = 40
JasparDataProcessing.JasparTxtToMat();
delete('../data/precomputation/SelectedPWMs.mat');
delete('../data/precomputation/pcPWMp.mat');
delete('../data/precomputation/pretrainedTheta.mat')
mainGenSequences(3000, 500, 4, 1000, false);
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized_fake.mat');
backgroundIndex = 4;
selectedPWMs = JasparDataProcessing.mainPreprocessPWMs(0.25, 0.55, mergedPeaksMin, [1,2,3,backgroundIndex], k)
pretrain(mergedPeaksMin, k, [1,2,3], backgroundIndex)
mainGenSequences(300, 4000, 4, 1000, true);
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized_fake.mat');
mainRealData(mergedPeaksMin, 4, k, false, true, false);

% REAL
delete('../data/precomputation/SelectedPWMs.mat');
delete('../data/precomputation/pcPWMp.mat');
delete('../data/precomputation/pretrainedTheta.mat')
peaks.beds2mats(500)
peaks.mergePeakFiles(true, true)
load('../data/peaks/mergedPeaks.mat');
peaks.minimizeMergePeak(mergedPeaks, 500, tissueNames)
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat');
geneIndex = find(strcmp('genes', tissueNames))
backgroundIndex = find(strcmp('background', tissueNames))
pretrain(mergedPeaksMin, k, [10, 20, 30], backgroundIndex)

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
