% m - number of base states in the model
% k - number of sub states each base state has in the model
function params = genParams(m, k)

    params.m = m;
    params.order = 3;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();

    params.backgroundRatio = 0.90;
    params.crossEnhancer = 0.02; % probability that for a cross enhancer transition at end of enhancer
    params.enhancerLength = 500;
    params.enhancerMotifsRatio = 0.10; % max limit

    meanEnhancerLength = (params.enhancerLength * (1 + params.crossEnhancer));
    PTotalBaseToOthers = 1 / ((1 - params.enhancerMotifsRatio) * params.enhancerLength);

    params.PTotalBaseToSub = misc.ratio2TransitionProb(mean(params.lengths, 2), params.enhancerMotifsRatio);
    params.PBackgroundToEnhancer = misc.ratio2TransitionProb(meanEnhancerLength, 1 - params.backgroundRatio) / (params.m - 1);
    params.PEnhancerToBackground = PTotalBaseToOthers * (1 - params.crossEnhancer);
    params.PCrossEnhancers = PTotalBaseToOthers * params.crossEnhancer / (params.m - 2);
    if params.m == 2
        params.PCrossEnhancers = 0;
    end
    params.maxPRatio = 0.5;




    params.batchSize = 2;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);

    params.minEnhLen = 300;
    params.minBgLen = 2000;
    params.maxEnhLen = 1000;
    params.maxBgLen = 10000;
    params.maxMotif = (1/10) ./ params.k;
    params.maxBgMotif = (1/2000) ./ params.k;
    params.minCrossEnh = eps;
    params.minEnhMotif = eps;
    params.minBgMotif = eps;



    [params.maxT, params.maxG] = genMaxGT(params);
    [params.minT, params.minG] = genMinGT(params);
    params
    try
        loaded = load('../data/precomputation/SelectedPWMs.mat');;
        inds = loaded.selectedPWMs;
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};
        fprintf('loaded optimized PWMs for tissue list', loaded.tissueList)
    catch
        params.PWMs = params.PWMs(1:params.k, :, :);
        params.lengths = params.lengths(1:params.k);
        params.names = {params.names{1:params.k}};
        fprintf('loaded non-optimized PWMs')
    end
end



function [minT, minG] = genMinGT(params)

    minStayEnh = 1 - (1 / params.minEnhLen); %denominator ~ min Enh length
    minStayBg = 1 - (1 / params.minBgLen); %denominator ~ min bg length
    minBgToEnh = 1 / params.maxBgLen; %denominator ~ max bg length
    minEnhToBg = 1 / params.maxEnhLen; %denominator ~ max Enh length
    params.minCrossEnh = eps;

    minT = ones(params.m) .* params.minCrossEnh;
    minT(eye(params.m) == 1) = minStayEnh;
    minT(:, params.m) = minEnhToBg;
    minT(params.m, :) = minBgToEnh;
    minT(params.m, params.m) = minStayBg;
    minG = ones(params.m, params.k) .* params.minEnhMotif;
    minG(params.m, :) = params.minBgMotif;
end

function [maxT, maxG] = genMaxGT(params)

    maxStayEnh = 1 - (1 / params.maxEnhLen); %numenator ~ max Enh length
    maxStayBg = 1 - (1 / params.maxBgLen); %numenator ~ max bg length
    maxBgToEnh = 1 / params.minBgLen; %numenator ~ min bg length
    maxEnhToBg = 1 / params.minEnhLen; %numenator ~ min Enh length
    maxCrossEnh = 1 / params.minEnhLen;

    maxT = ones(params.m) .* maxCrossEnh;
    maxT(eye(params.m) == 1) = maxStayEnh;
    maxT(:, params.m) = maxEnhToBg;
    maxT(params.m, :) = maxBgToEnh;
    maxT(params.m, params.m) = maxStayBg;
    maxG = ones(params.m, params.k) .* params.maxMotif;
    maxG(params.m, :) = params.maxBgMotif;

end