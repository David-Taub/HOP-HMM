% m - number of base states in the model
% k - number of sub states each base state has in the model
% backgroundAmount - number of background states
function params = genParams(m, k, backgroundAmount, L, order, doESharing)
    params.doESharing = doESharing;
    params.m = m;
    params.k = k;
    params.order = order;
    params.backgroundAmount = backgroundAmount;
    params = loadPWMs(params);

    % params.enhancerLength = 500;
    % meanEnhancerCountInSeq = 1.5;
    % % percent of background out of the sequence
    % params.backgroundRatio = 1 - ((params.enhancerLength * meanEnhancerCountInSeq) / L);
    % % probability that for a cross enhancer transition at end of enhancer
    % params.crossEnhancer = 0.02;
    % params.enhancerMotifsRatio = 0.10; % max limit

    % PTotalBaseToOthers = 1 / ((1 - params.enhancerMotifsRatio) * params.enhancerLength);

    % params.PTotalBaseToSub = misc.ratio2TransitionProb(mean(params.lengths, 2), params.enhancerMotifsRatio);
    % params.PBackgroundToEnhancer = misc.ratio2TransitionProb(params.enhancerLength, 1 - params.backgroundRatio) / ((params.m - 1) * params.backgroundAmount);
    % params.PEnhancerToBackground = PTotalBaseToOthers * (1 - params.crossEnhancer) / params.backgroundAmount;
    % params.PCrossEnhancers = PTotalBaseToOthers * params.crossEnhancer / (params.m - 1 - params.backgroundAmount);

    % if params.m - params.backgroundAmount == 1
    %     params.PCrossEnhancers = 0;
    % end

    % params.maxPRatio = 3;
    params.batchSize = 1000;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);
    params.EPS = 10 ^ -7;
    maxBbMotifRatio = 1 / 90; % maximal ratio of all motif letters in BG section
    maxEnhMotifRatio = 1 / 10; % maximal ratio of a single motif letters in enhancer section
    maxCrossEnhRatio = 1 / 15; % maximal ratio of transition to any another enhancer
    maxCrossBgRatio = 1 / 20; % maximal ratio of transition to any another enhancer
    maxEnhLen = min(700, floor(0.7 * L));
    minEnhLen = min(100, floor(0.3 * L));
    maxBgLen = L - minEnhLen;
    minBgLen = L - maxEnhLen;

    params.maxEnhMotif = maxEnhMotifRatio / mean(params.lengths);
    params.minEnhMotif = params.EPS;
    params.maxCrossEnh = maxCrossEnhRatio / (minEnhLen * (params.m - params.backgroundAmount))
    params.minCrossEnh = params.EPS;
    params.maxBgMotif = maxBbMotifRatio / ( mean(params.lengths) * params.k);
    params.minBgMotif = params.EPS;
    params.maxStayEnh = 1 - (1 / maxEnhLen);
    params.minStayEnh = 1 - (1 / minEnhLen);
    params.maxStayBg = 1 - (1 / maxBgLen);
    params.minStayBg = 1 - (1 / minBgLen);
    params.maxCrossBg = 1 / minBgLen;
    params.minCrossBg = params.EPS;
    params.maxBgToEnh = 1 / minBgLen;
    params.minBgToEnh = 1 / maxBgLen;
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
    for t = params.m - params.backgroundAmount + 1:params.m
        minT(t, t) = params.minStayBg;
    end
    minG = ones(params.m, params.k) .* params.minEnhMotif;
    minG(end - params.backgroundAmount + 1:end, :) = params.minBgMotif;
end


function [maxT, maxG] = genMaxGT(params)
    maxT = ones(params.m) .* params.maxCrossEnh;
    maxT(end - params.backgroundAmount + 1: end, end - params.backgroundAmount + 1: end) = params.maxCrossBg;
    maxT(eye(params.m) == 1) = params.maxStayEnh;
    maxT(:, end - params.backgroundAmount + 1: end) = params.maxEnhToBg;
    maxT(end - params.backgroundAmount + 1:end, :) = params.maxBgToEnh;
    for t = params.m - params.backgroundAmount + 1:params.m
        maxT(t, t) = params.maxStayBg;
    end
    maxG = ones(params.m, params.k) .* params.maxEnhMotif;
    maxG(end - params.backgroundAmount + 1:end, :) = params.maxBgMotif;
end


function params = loadPWMs(params)
    [params.PWMs, params.lengths, params.names] = misc.PWMs();
    try
        SELECTED_PWMS_FILEPATH = '../data/precomputation/SelectedPWMs.mat';
        loaded = load(SELECTED_PWMS_FILEPATH);
        inds = loaded.selectedPWMs;
        params.k = length(inds);
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};

        % k x n x J -> J x n x k
        PWMImage = permute(params.PWMs, [3,2,1]);
        PWMImage = cat(2, PWMImage, zeros(params.J, 4, params.k ));
        PWMImage = PWMImage(:, :);
        figure;
        imagesc(PWMImage);
        title('PWMs');
        drawnow;
        fprintf('loaded feature selected PWMs from %s\n', SELECTED_PWMS_FILEPATH)
    catch
        params.PWMs = params.PWMs(1:params.k, :, :);
        params.lengths = params.lengths(1:params.k);
        params.names = {params.names{1:params.k}};
        fprintf('loaded non-optimized PWMs');
    end
end
