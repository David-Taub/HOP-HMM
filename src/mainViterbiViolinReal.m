
function mainViterbiViolin()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.doResampling = false;
    conf.maxIters = 1000;
    conf.canCrossLayer = true;
    conf.patience = 4;
    conf.L = 1000;
    conf.N = 100;
    conf.withExponent = false;
    conf.repeat = 1;
    conf.order = 2;
    conf.m = 5;
    conf.k = 10;
    conf.backgroundAmount = 1;
    conf.doGTBound = false;
    main(conf);

end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing, conf.doGTBound, conf.doResampling);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    outpath = misc.pathMaker(params, conf.N, conf.L, 'viterbiViolin', '.jpg');
    subtitle = sprintf('m=%d, k=%d, %d%% of data', conf.m, conf.k);
    dataset.title = subtitle;
    dataset.X = mergedPeaksMin.seqs;
    dataset.theta = mergedPeaksMin.theta;
    dataset.pcPWMp = misc.preComputePWMp(mergedPeaksMin.seqs, params);
    [thetaEst, ~] = EM.EM(dataset, params, conf.maxIters, conf.patience, conf.repeat);
    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    [~, ~, ~, ~, ~, psi] = EM.EStep(params, thetaEst, dataset.X, dataset.pcPWMp);
    % N x L x 2
    YEstViterbi = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);

    estMask1 = YEstViterbi(:, :, 1) == 1;
    estMask2 = YEstViterbi(:, :, 1) == 2;
    estMask3 = YEstViterbi(:, :, 1) == 3;

    bedGraphs = misc.readAllBedGraphs(mergedPeaksMin.tissueEIDs, {'DNase'});
    % N x L x tissues
    % TODO
    tracks = getTracks(bedGraphs)

    estVals1 = tracks(estMask1, 1);
    estVals2 = tracks(estMask2, 1);
    estVals3 = tracks(estMask3, 1);
    violinViterbi3(estVals1, estVals2, estVals3);

    estVals1 = tracks(estMask1, 2);
    estVals2 = tracks(estMask2, 2);
    estVals3 = tracks(estMask3, 2);
    % TODO
    violinViterbi3(estVals1, estVals2, estVals3);
end


function violinViterbi3(truePosVals, trueNegVals, estPosVals, estNegVals, outpath)
    fig = figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    maxRange = max([truePosVals; trueNegVals; estPosVals; estNegVals], [], 1);
    minRange = min([truePosVals; trueNegVals; estPosVals; estNegVals], [], 1);
    show.distributionPlot({[truePosVals; minRange; maxRange], [trueNegVals; minRange; maxRange]}, 'histOri', 'right', 'color', 'r', 'widthDiv', [2 2], 'showMM', 0, 'histOpt', 0, 'divFactor', [minRange: maxRange]);
    show.distributionPlot({[estPosVals; minRange; maxRange], [estNegVals; minRange; maxRange]}, 'histOri', 'left', 'color', 'b', 'widthDiv', [2 1], 'showMM', 0, 'histOpt', 0, 'divFactor', [minRange: maxRange]);
    jf = java.text.DecimalFormat;
    tick1 = sprintf('TFBS Positions [n=%s, n=%s]', char(jf.format(length(truePosVals))), char(jf.format(length(estPosVals))));
    tick2 = sprintf('Non-TFBS Positions [n=%s, n=%s]', char(jf.format(length(trueNegVals))), char(jf.format(length(estNegVals))));
    set(gca,'xtick', [1:2], 'xticklabel', {tick1, tick2});
    ylabel({'Log Posterior';'Probability (log\gamma)'});
    ylh = get(gca,'ylabel');
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation', 0, 'Position', ylp-[ext(3) / 2 0 0], 'VerticalAlignment','middle', 'HorizontalAlignment','center')
    legend({'True TFBS positions', 'Viterbi estimated TFBS positions'})
    title('Posterior Probability of TFBS and non-TFBS Locations');
    saveas(fig, outpath);
end