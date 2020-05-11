
function mainViterbiViolin()
    conf.startWithBackground = false;
    conf.doEnhSpecific = true;
    conf.seqsPerTissue = 1000;
    conf.maxIters = 25;
    conf.repeat = 1;
    conf.canCrossLayer = true;
    conf.patience = 4;
    conf.L = 2000;
    conf.peakMaxLength = 1000;
    conf.peakMinL = 200;
    conf.peakMaxL = 1500;
    conf.withExponent = false;
    conf.order = 3;
    conf.m = 3;
    conf.k = 20;
    conf.withBackground = true;
    conf.withGenes = false;
    conf.minSamplesCount = 10;
    conf.sequencesToShow = 5;
    conf.backgroundAmount = 1;
    conf.doESharing = false;
    conf.doGTBound = true;
    conf.doResampling = false;
    conf.topPercent = 0.5;
    % conf.tissueList = [2, 37];
    % conf.tissueList = [3, 23];
    conf.tissueList = [2, 18];
    conf.startTUniform = false;
    main(conf);
end

function main(conf)
    dbstop if error
    close all;
    mergedPeaksMin = peaks.minimizeMergePeak(conf.topPercent, conf.doEnhSpecific, conf.withBackground, conf.withGenes,...
                                             conf.seqsPerTissue, conf.L, conf.peakMinL, conf.peakMaxL, conf.tissueList,...
                                             conf.minSamplesCount);
    N = size(mergedPeaksMin.seqs, 1);
    testTrainRatio = 0.01;
    selectedPWMs = misc.PWMsFeatureSelect(mergedPeaksMin, conf.k);
    params = misc.genParams(conf.m, selectedPWMs, conf.backgroundAmount, conf.L, conf.order, ...
                            conf.doESharing, conf.doGTBound);
    [test, train] = misc.crossValidationSplit(params, mergedPeaksMin, testTrainRatio);
    [thetaEst, ~] = EM.EM(train, params, conf.maxIters, conf.patience, conf.repeat);
    % N x L x 2
    plotViolins(params, thetaEst, train, mergedPeaksMin.tissueEIDs)
end

function plotViolins(params, thetaEst, dataset, tissueEIDs)
    YEstViterbi = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);
    % estMask - N x L
    estMask1 = YEstViterbi(:, :, 1) == 1;
    estMask2 = YEstViterbi(:, :, 1) == 2;
    estMask3 = YEstViterbi(:, :, 1) == 3;
    % bedGraphs - tissues x tracks
    bedGraphs = misc.readAllBedGraphs(tissueEIDs, {'DNase'});
    % N x L x tissues
    tracks = getTracks(dataset, bedGraphs, 1);
    % estVals - N x L
    tracks1 = tracks(:, :, 1);
    estVals1 = tracks1(estMask1);
    estVals2 = tracks1(estMask2);
    estVals3 = tracks1(estMask3);
    outpath = misc.pathMaker(params, size(estMask1, 1), size(estMask1, 2), 0, '3violine_1', '.jpg');
    violinViterbi3(estVals1, estVals2, estVals3, outpath);

    tracks2 = tracks(:, :, 2);
    estVals1 = tracks2(estMask1);
    estVals2 = tracks2(estMask2);
    estVals3 = tracks2(estMask3);
    outpath = misc.pathMaker(params, size(estMask1, 1), size(estMask1, 2), 0, '3violine_2', '.jpg');
    violinViterbi3(estVals1, estVals2, estVals3, outpath);
end


% TODO
% estVals1 - N1, 1
% estVals2 - N2, 1
% estVals3 - N3, 1

function violinViterbi3(estVals1, estVals2, estVals3, outpath)
    % outpath = 'c:\tmp\tmp.jpg'
    % estVals1 = r1;
    % estVals2 = r2;
    % estVals3 = r3;
    fig = figure('units', 'pixels', 'Position', [0 0 1500 1000]);
    maxRange = max([estVals1; estVals2; estVals3], [], 1);
    minRange = min([estVals1; estVals2; estVals3], [], 1);
    % show.distributionPlot({estVals1, estVals2, estVals3}, 'showMM', 0);
    show.distributionPlot(estVals1, 'showMM', 0, 'color', 'r', 'xValues', 1);
    show.distributionPlot(estVals2, 'showMM', 0, 'color', 'b', 'xValues', 2);
    show.distributionPlot(estVals3, 'showMM', 0, 'color', 'k', 'xValues', 3);
    jf = java.text.DecimalFormat;
    set(gca,'xtick', [1:3], 'xticklabel', {'State 1', 'State 2', 'State 3'});
    ylabel('DNase -log(p-value)');
    ylh = get(gca,'ylabel');
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation', 0, 'Position', ylp-[ext(3) / 2 0 0], 'VerticalAlignment','middle', 'HorizontalAlignment','center')
    legend({'True TFBS positions', 'Viterbi estimated TFBS positions'})
    title('Posterior Probability of TFBS and non-TFBS Locations');
    saveas(fig, outpath);
end

% bedGraphs - tissues x tracks
% tracks - N x L x tissues
function tracks = getTracks(dataset, bedGraphs, trackInd)
    [N, L] = size(dataset.X);
    tissues = size(bedGraphs, 1);
    tracks = zeros(N, L, tissues);
    for tissueInd = 1:size(bedGraphs, 1)
        fprintf('getting tracks for tissue index %d\n', tissueInd)
        for seqInd = 1:N
            chr = dataset.chrs{seqInd};
            to = dataset.starts(seqInd) + L;
            from = dataset.starts(seqInd);
            tracks(seqInd, :, tissueInd) = getTrack(bedGraphs{tissueInd, trackInd}, chr, from, to);
            % if all(tracks(i, :, trackInd) == 0)
            %     fprintf('Seq has not enough samples\n', trackInd)
            %     return;
            % end
        end
    end
end

% ret - L x 1
function ret = getTrack(bedGraph, trackChr, trackFrom, trackTo)
    MIN_SAMPLES_IN_SEQ = 100;
    mask = strcmp(bedGraph.chrs, trackChr) & (bedGraph.tos >= trackFrom) & (bedGraph.froms <= trackTo);
    % dilation of bin vector
    mask = [mask(2:end); false] | mask | [false; mask(1:end-1)];
    % linear interp
    foundSamples = sum(mask);
    if foundSamples == 0
        ret = zeros(trackTo - trackFrom, 1);
        return
    end
    knownPoints = double([(bedGraph.tos(mask) + bedGraph.froms(mask)) / 2]');
    knownVals = bedGraph.vals(mask)';
    wantedPoints = double([trackFrom : trackTo - 1]);
    ret = interp1(knownPoints, knownVals, wantedPoints);
end

