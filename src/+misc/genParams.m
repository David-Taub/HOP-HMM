% m - number of base states in the model
% k - number of sub states each base state has in the model
function params = genParams(m, k)

    params.m = m;
    params.order = 3;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();

    params.backgroundRatio = 0.7;
    params.crossEnhancer = 0.1; % probability that for a cross enhancer transition at end of enhancer
    params.enhancerLength = 500;
    params.enhancerMotifsRatio = 0.1;

    meanEnhancerLength = (params.enhancerLength * (1 + params.crossEnhancer));
    PTotalBaseToOthers = 1/((1-params.enhancerMotifsRatio)*params.enhancerLength);

    params.PTotalBaseToSub = misc.ratio2TransitionProb(mean(params.lengths, 2), params.enhancerMotifsRatio);
    params.PBackgroundToEnhancer = misc.ratio2TransitionProb(meanEnhancerLength, 1-params.backgroundRatio) / (params.m-1);
    params.PEnhancerToBackground = PTotalBaseToOthers * (1-params.crossEnhancer);
    params.PCrossEnhancers = PTotalBaseToOthers * params.crossEnhancer / (params.m-2);
    params.maxPRatio = 5;
    params


    params.batchSize = 2;
    [kMax, params.n, params.J] = size(params.PWMs);
    params.k = min(k, kMax);
    try
        loaded = load('../data/temp/G.mat');
        inds = loaded.inds;
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};
    catch
        params.PWMs = params.PWMs(1:params.k, :, :);
        params.lengths = params.lengths(1:params.k);
        params.names = {params.names{1:params.k}};
    end
end