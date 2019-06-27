1. Download data + preprocess: D:\projects\HOP-HMM\data\peaks\scripts\download_and_process_all.sh
2. In Matlab - JASPAR preprocess JasparDataProcessing.JasparTxtToMat();
3. In Matlab - Gen synthetic data + plots via the main*.m files in src
4. In Matlab - Real data:
delete('../data/precomputation/SelectedPWMs.mat');
delete('../data/precomputation/pcPWMp.mat');
delete('../data/precomputation/pretrainedTheta.mat');
k = 40
peaks.beds2mats(500)
peaks.mergePeakFiles(true, true)
load('../data/peaks/mergedPeaks.mat');
peaks.minimizeMergePeak(mergedPeaks, 500, tissueNames)
mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat');
geneIndex = find(strcmp('genes', tissueNames))
backgroundIndex = find(strcmp('background', tissueNames))
m = 3;
tissueList = [14, 75, backgroundIndex];
tissueList = mainSelectTissues(mergedPeaksMin, backgroundIndex, m);
selectedPWMs = JasparDataProcessing.mainPreprocessPWMs(0.25, 0.55, mergedPeaksMin, tissueList, k);
pretrain(mergedPeaksMin, k, tissueList(1:end-1), backgroundIndex)

peaks.beds2matsNoSeq()
peaks.mergePeakFiles(false, false)
load('../data/peaks/mergedPeaksNoBackground.mat');
multiEnhancers = peaks.multiEnhancerCaller(mergedPeaks, 10000, tissueList(1:end-1));
load('../data/multiEnhancers.mat');
k = 40
mainRealData(multiEnhancers, 4, k, false, true, false);


mergedPeaksMin =
        lengths: [121279×1 int32]
       overlaps: [121279×46 double]
           seqs: [121279×500 double]
    tissueNames: {1×46 cell}
