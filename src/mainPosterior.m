
function  mainPosterior()
    conf.doESharing = false;
    conf.startWithBackground = true;
    conf.maxIters = 100;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1500;
    conf.N = 1000;
    conf.withExponent = false;
    conf.repeat = 3;

    conf.order = 3;
    conf.m = 5;
    conf.k = 10;
    conf.backgroundAmount = 1;
    conf.doBound = false;
    conf.doResample = false;
    main(conf);
end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15);
    [thetaEst, ~] = EM.EM(trainDataset, params, conf.maxIters, conf.doBound, conf.doResample, conf.patience, conf.repeat);

    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    outpath = sprintf('Posterior_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.N, conf.L);
    show.seqSampleCertainty(params, thetaEst, testDataset, 3, outpath);
    outpath = sprintf('Posterior2_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.N, conf.L);
    show.seqSampleCertainty(params, thetaEst, testDataset, 3, outpath);
    outpath = sprintf('Posterior3_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.N, conf.L);
    show.seqSampleCertainty(params, thetaEst, testDataset, 3, outpath);
    show.showTheta(thetaEst);
    show.showTheta(thetaOrig);
end
