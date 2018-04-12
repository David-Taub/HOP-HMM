
% first we choose the PWMs that are most different between the following tissues,
% including background and genes:

peaks.beds2mats(500)
peaks.mergePeakFiles(true, true)
load('../data/peaks/mergedPeaks.mat');
peaks.minimizeMergePeak(mergedPeaks, 500, tissueNames)
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat');
JasparDataProcessing.mainPreprocessPWMs()
chooseBestPWMs(mergedPeaksMin, [10, 20, 30, 45, 46])

% then we find multiEnhancers, and try to learn them
peaks.beds2matsNoSeq()
peaks.mergePeakFiles(false, false)
load('../data/peaks/mergedPeaksNoBackground.mat');
multiEnhancers = peaks.multiEnhancerCaller(mergedPeaks, 10000, [10, 20, 30]);
load('../data/multiEnhancers.mat');
mainRealData(multiEnhancers, 4, 40, false, true);