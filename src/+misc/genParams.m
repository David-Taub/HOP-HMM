% m - number of base states in the model
% k - number of sub states each base state has in the model
% backgroundAmount - number of background states
function params = genParams(m, k, backgroundAmount, L, order, doESharing, doGTBound,...
                            doResampling)
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
    params.doResampling = doResampling;

    params.batchSize = 100;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);
    params.EPS = 10 ^ -7;
    maxBbMotifRatio = 1 / 90; % maximal ratio of all motif letters in BG section
    maxEnhMotifRatio = 1 / 20; % maximal ratio of a single motif letters in enhancer section
    % minEnhMotifRatio = 1 / 15; % minimal ratio of a single motif letters in enhancer section
    maxCrossEnhRatio = 1 / 100; % maximal ratio of transition to any another enhancer
    maxCrossBgRatio = 1 / 30; % maximal ratio of transition to any another enhancer
    maxEnhLen = 1000;
    minEnhLen = 400;
    maxBgLen = 4000;
    minBgLen = 500;

    params.maxEnhMotif = maxEnhMotifRatio / mean(params.lengths);
    % params.minEnhMotifTotal = minEnhMotifRatio / mean(params.lengths);
    params.minEnhMotif = params.EPS;
    params.maxCrossEnh = maxCrossEnhRatio / minEnhLen;
    params.minCrossEnh = params.EPS;
    params.maxBgMotif = maxBbMotifRatio / ( mean(params.lengths) * params.k);
    params.minBgMotif = params.EPS;
    params.maxStayEnh = 1;
    params.minStayEnh = 1 - (1 / minEnhLen);
    params.maxStayBg = 1 - (1 / maxBgLen);
    params.minStayBg = 1 - (1 / minBgLen);
    params.maxCrossBg = 1 / minBgLen;
    params.minCrossBg = params.EPS;
    params.maxBgToEnh = 1 / minBgLen;
    params.minBgToEnh = 1 / maxBgLen
    params.maxEnhToBg = 1 / minEnhLen;
    params.minEnhToBg = 1 / (maxEnhLen * params.backgroundAmount);


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
    minG(end - params.backgroundAmount + 1:end, :) = params.minBgMotif;
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
    maxG(end - params.backgroundAmount + 1:end, :) = params.maxBgMotif;
end