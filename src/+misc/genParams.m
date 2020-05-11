% m - number of base states in the model
% k - number of sub states each base state has in the model
% backgroundAmount - number of background states
function params = genParams(m, k, backgroundAmount, L, order, doESharing, doGTBound)
    params.doESharing = doESharing;
    params.m = m;
    if length(k) > 1
        selectedPWMs = k;
        k = length(selectedPWMs);
    else
        selectedPWMs = 1:k;
    end
    params.k = k;
    params.order = order;
    params.backgroundAmount = backgroundAmount;
    params.enhancerAmount = m - backgroundAmount;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();
    params.PWMs = params.PWMs(selectedPWMs, :, :);
    params.lengths = params.lengths(selectedPWMs);
    params.names = {params.names{selectedPWMs}};
    params.learningRate = 1;

    params.doGTBound = doGTBound;
    params.doResampling = false;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);
    params.batchSize = 100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Real Sequences Params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % params.topPercent = 0.8;
    % params.doEnhSpecific = true;
    % params.withBackground = true;
    % params.withGenes = false;
    % params.seqsPerTissue = 1000;
    % params.peakMinL = 100;
    % params.peakMaxL = 1500;
    % params.tissueList = [2, 18];
    params.minSamplesCount = 5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regularization Params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.EPS = 10 ^ -7;
    maxBbMotifRatio = 1 / 10000; % maximal ratio of all motif letters in BG section
    maxEnhMotifRatio = 1 / 1000; % maximal ratio of a single motif letters in enhancer section
    % minEnhMotifRatio = 1 / 15; % minimal ratio of a single motif letters in enhancer section
    maxCrossEnhRatio = 1 / 100000; % maximal ratio of transition to any another enhancer
    maxCrossBgRatio = 1 / 3000; % maximal ratio of transition to any another enhancer
    maxEnhLen = 1000;
    minEnhLen = 200;
    maxBgLen = 2000;
    minBgLen = 500;

    params.ACTIVE_TFS = 3;
    params.maxEnhMotif = maxEnhMotifRatio;
    % params.minEnhMotifTotal = minEnhMotifRatio / mean(params.lengths);
    params.minEnhMotif = params.EPS;
    params.maxCrossEnh = maxCrossEnhRatio;
    params.minCrossEnh = params.EPS;
    params.maxBgMotif = maxBbMotifRatio / (mean(params.lengths) * params.k);
    params.minBgMotif = params.EPS;
    params.maxStayEnh = 1 - (1 / maxEnhLen);
    params.minStayEnh = 1 - (1 / 100);
    params.maxStayBg = 1 - (1 / minBgLen);
    params.minStayBg = 1 - (1 / maxBgLen);
    params.maxCrossBg = 1 / minBgLen;
    params.minCrossBg = params.EPS;
    params.maxBgToEnh = 1 / (minBgLen * params.m);
    params.minBgToEnh = 1 / (maxBgLen * params.m);
    params.maxEnhToBg = 1 / minEnhLen;
    params.minEnhToBg = 1 / maxEnhLen;


    [params.maxT, params.maxG] = genMaxGT(params);
    [params.minT, params.minG] = genMinGT(params);
end


function [minT, minG] = genMinGT(params)
    minT = ones(params.m) .* params.minCrossEnh;
    minT(end - params.backgroundAmount + 1: end, end - params.backgroundAmount + 1: end) = params.minCrossBg;
    minT(eye(params.m) == 1) = params.minStayEnh;
    minT(:, end - params.backgroundAmount + 1: end) = params.minEnhToBg;
    minT(end - params.backgroundAmount + 1:end, :) = params.minBgToEnh;
    for t = params.enhancerAmount + 1:params.m
        minT(t, t) = params.minStayBg;
    end
    minG = ones(params.m, params.k) .* params.minEnhMotif;
    % minG(end - params.backgroundAmount + 1:end, :) = params.minBgMotif;
    % assert(all(sum(minG, 2) + sum(minT, 2) < 1))
end


function [maxT, maxG] = genMaxGT(params)
    maxT = ones(params.m) .* params.maxCrossEnh;
    maxT(end - params.backgroundAmount + 1: end, end - params.backgroundAmount + 1: end) = params.maxCrossBg;
    maxT(eye(params.m) == 1) = params.maxStayEnh;
    maxT(:, end - params.backgroundAmount + 1: end) = params.maxEnhToBg;
    maxT(end - params.backgroundAmount + 1:end, :) = params.maxBgToEnh;
    for t = params.enhancerAmount + 1:params.m
        maxT(t, t) = params.maxStayBg;
    end
    maxG = ones(params.m, params.k) .* params.maxEnhMotif;
end