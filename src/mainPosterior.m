% generation of figures 20 and 21
function  mainPosterior()
    conf.doESharing = false;
    conf.startWithBackground = true;
    conf.maxIters = 220;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1500;
    conf.N = 1000;
    conf.withExponent = false;
    conf.repeat = 6;

    conf.background_g_noise = 0.007;

    conf.order = 3;
    conf.m = 5;
    conf.k = 10;
    conf.backgroundAmount = 1;
    conf.doGTBound = true;
    conf.doResampling = false;
    main(conf);
end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, ...
                            conf.doESharing, conf.doGTBound, conf.doResampling);
    mergedPeaksMin = misc.genSyntheticMergedPeaksMin(conf.N, conf.L, params, conf.startWithBackground, conf.background_g_noise);
    thetaOrig = mergedPeaksMin.theta;
    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15);
    [thetaEst, ~] = EM.EM(trainDataset, params, conf.maxIters, conf.patience, conf.repeat, false);

    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    outpath = sprintf('seqSample_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doGTBound, conf.N, conf.L);
    show.seqSampleCertainty(params, thetaEst, testDataset, 3, outpath);
    show.showTheta(thetaEst);
    show.showTheta(thetaOrig);
end
