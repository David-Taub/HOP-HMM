% m - number of base states in the model
% k - number of sub states each base state has in the model
% backgroundAmount - number of background states
function params = genParams(m, k, backgroundAmount)

    params.m = m;
    params.order = 3;
    params.backgroundAmount = backgroundAmount;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();

    params.backgroundRatio = 0.90;
    params.crossEnhancer = 0.02; % probability that for a cross enhancer transition at end of enhancer
    params.enhancerLength = 500;
    params.enhancerMotifsRatio = 0.10; % max limit

    meanEnhancerLength = (params.enhancerLength * (1 + params.crossEnhancer));
    PTotalBaseToOthers = 1 / ((1 - params.enhancerMotifsRatio) * params.enhancerLength);

    params.PTotalBaseToSub = misc.ratio2TransitionProb(mean(params.lengths, 2), params.enhancerMotifsRatio);
    params.PBackgroundToEnhancer = misc.ratio2TransitionProb(meanEnhancerLength, 1 - params.backgroundRatio) / ((params.m - 1) * params.backgroundAmount);
    params.PEnhancerToBackground = PTotalBaseToOthers * (1 - params.crossEnhancer) / params.backgroundAmount;
    params.PCrossEnhancers = PTotalBaseToOthers * params.crossEnhancer / (params.m - 1 - params.backgroundAmount);

    if params.m - params.backgroundAmount == 1
        params.PCrossEnhancers = 0;
    end

    params.maxPRatio = 3;
    params.batchSize = 50;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);


    params.minEnhLen = 300;
    params.minBgLen = 2000;
    params.maxEnhLen = 1000;
    params.maxBgLen = 10000;
    params.maxMotif = (1/10) ./ params.k;
    params.maxBgMotif = (1/2000) ./ params.k;
    params.minCrossEnh = eps;
    params.maxCrossEnh = 1 / (params.minEnhLen * (params.m - params.backgroundAmount) * 10)
    params.minEnhMotif = eps;
    params.minBgMotif = eps;



    [params.maxT, params.maxG] = genMaxGT(params);
    [params.minT, params.minG] = genMinGT(params);
    params
    try
        selectedPWMsFilepath = '../data/precomputation/SelectedPWMs.mat';
        loaded = load(selectedPWMsFilepath);
        inds = loaded.selectedPWMs;
        params.k = length(inds);
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};

        % k x n x J -> J x n x k
        PWMImage = permute(params.PWMs, [3,2,1]);
        PWMImage = cat(2, PWMImage, zeros(params.J, 4, params.k ));
        PWMImage = PWMImage(:, :);
        figure
        imagesc(PWMImage);
        title('PWMs');
        drawnow;
        fprintf('loaded feature selected PWMs from %s\n', selectedPWMsFilepath)
    catch
        params.PWMs = params.PWMs(1:params.k, :, :);
        params.lengths = params.lengths(1:params.k);
        params.names = {params.names{1:params.k}};
        fprintf('loaded non-optimized PWMs')
    end
end



function [minT, minG] = genMinGT(params)
    minStayEnh = 1 - (1 / params.minEnhLen);
    minStayBg = 1 - (params.backgroundAmount / params.minBgLen);
    minCrossBg = params.backgroundAmount / params.maxBgLen;
    minBgToEnh = 1 / params.maxBgLen;
    minEnhToBg = 1 / (params.maxEnhLen * params.backgroundAmount);

    minT = ones(params.m) .* params.minCrossEnh;
    minT(params.m - params.backgroundAmount + 1: end, params.m - params.backgroundAmount + 1: end) = minCrossBg;
    minT(eye(params.m) == 1) = minStayEnh;
    minT(:, params.m - params.backgroundAmount + 1: end) = minEnhToBg;
    minT(params.m - params.backgroundAmount + 1:end, :) = minBgToEnh;
    for t = params.m - params.backgroundAmount + 1:params.m
        minT(t, t) = minStayBg;
    end
    minG = ones(params.m, params.k) .* params.minEnhMotif;
    minG(params.m - params.backgroundAmount + 1:end, :) = params.minBgMotif;
end

function [maxT, maxG] = genMaxGT(params)
    maxStayEnh = 1 - (1 / params.maxEnhLen);
    maxStayBg = 1 - (params.backgroundAmount / params.maxBgLen);
    maxCrossBg = params.backgroundAmount / params.minBgLen;
    maxBgToEnh = 1 / params.minBgLen;
    maxEnhToBg = 1 / (params.minEnhLen * params.backgroundAmount);

    maxT = ones(params.m) .* params.maxCrossEnh;
    maxT(params.m - params.backgroundAmount + 1: end, params.m - params.backgroundAmount + 1: end) = maxCrossBg;
    maxT(eye(params.m) == 1) = maxStayEnh;
    maxT(:, params.m - params.backgroundAmount + 1: end) = maxEnhToBg;
    maxT(params.m - params.backgroundAmount + 1:end, :) = maxBgToEnh;
    for t = params.m - params.backgroundAmount + 1:params.m
        maxT(t, t) = maxStayBg;
    end
    maxG = ones(params.m, params.k) .* params.maxMotif;
    maxG(params.m - params.backgroundAmount + 1:end, :) = params.maxBgMotif;
end