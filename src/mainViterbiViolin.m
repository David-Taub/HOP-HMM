
function mainViterbiViolin()
    conf.doESharing = false;
    conf.startWithBackground = false;
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
    conf.doBound = false;
    main(conf);

end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    outpath = sprintf('viterbi_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.N, conf.L);
    subtitle = sprintf('m=%d, k=%d, %d%% of data', conf.m, conf.k);
    dataset.title = subtitle;
    dataset.X = mergedPeaksMin.seqs;
    dataset.theta = mergedPeaksMin.theta;
    dataset.Y = mergedPeaksMin.Y;
    dataset.Y2 = mergedPeaksMin.Y2;
    dataset.pcPWMp = misc.preComputePWMp(mergedPeaksMin.seqs, params);
    [thetaEst, ~] = EM.EM(dataset, params, conf.maxIters, conf.doBound, conf.patience, conf.repeat);
    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    [~, ~, ~, ~, ~, psi] = EM.EStep(params, thetaEst, dataset.X, dataset.pcPWMp);
    % N x L x 2
    YEstViterbi = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);
    % N x L x 2
    Ymerged = cat(3, dataset.Y, dataset.Y2);

    % N x m x k x L
    estMask = genPWMMask(params, YEstViterbi, conf.N, conf.L) == 1;
    % N x m x k x L
    trueMask = genPWMMask(params, Ymerged, conf.N, conf.L) == 1;

    estPosVals = psi((estMask) & (psi > -inf));
    estNegVals = psi((~estMask) & (psi > -inf));

    truePosVals = psi((trueMask) & (psi > -inf));
    trueNegVals = psi((~trueMask) & (psi > -inf));
    show.violinViterbi(truePosVals, trueNegVals, estPosVals, estNegVals, outpath);
end


% Y - N x L x 2
% mask - N x m x k x L
function mask = genPWMMask(params, Y, N, L)
    mask = [];
    for l = 1:params.k
        layerMask = [];
        for i = 1:params.m
            PWMStateMask = (Y(:, 2:end, 2) == l) & (Y(:, 1:end - 1, 2) == 0) & (Y(:, 1:end - 1, 1) == i);
            % N x L
            PWMStateMask = cat(2, PWMStateMask, false(N, 1, 1));
            layerMask = cat(3, layerMask, PWMStateMask);
        end
        assert(all(size(layerMask) == [N, L, params.m]));
        layerMask = permute(layerMask, [1, 3, 4, 2]);
        mask = cat(3, mask, layerMask);
    end
    assert(all(size(mask) == [N, params.m, params.k, L]));
end
