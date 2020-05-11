% generation of figures 20 and 21
function  mainPosterior()
    conf.doESharing = false;
    conf.startWithBackground = true;
    conf.maxIters = 30;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1500;
    conf.N = 500;
    conf.withExponent = false;
    conf.repeat = 2;

    conf.background_g_noise = 1 / 50;

    conf.order = 3;
    conf.m = 5;
    conf.k = 25;
    conf.backgroundAmount = 1;
    conf.doGTBound = 5;
    conf.doResampling = false;
    main(conf);
end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, ...
                            conf.doESharing, conf.doGTBound);
    mergedPeaksMin = misc.genSyntheticMergedPeaksMin(conf.N, conf.L, params, conf.startWithBackground, conf.background_g_noise);
    thetaOrig = mergedPeaksMin.theta;
    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15);
    [thetaEst, ~] = EM.EM(trainDataset, params, conf.maxIters, conf.patience, conf.repeat);

    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    outpath = sprintf('seqSample_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doGTBound, conf.N, conf.L);
    distMat = pdist2(exp(misc.thetaToMat(params, thetaOrig, false)), ...
                     exp(misc.thetaToMat(params, thetaEst, false)))
    perm = misc.munkres(distMat)';
    thetaEstPerm = misc.permTheta(thetaEst, perm);

    show.showTheta(thetaEstPerm);
    show.showTheta(thetaOrig);
    show.seqSampleCertainty(params, thetaEstPerm, testDataset, 3, outpath);
end
