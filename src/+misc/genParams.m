% m - number of base states in the model
% k - number of sub states each base state has in the model
function params = genParams(m, k)

    params.m = m;
    params.order = 3;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();
    params.tEpsilon = 1/1000;
    params.backgroundRatio = 0.9;
    params.crossEnhancer = 0.1;
    params.enhancerLength = 500;
    params.enhancerMotifsRatio = 0.2;
    params.batchSize = 2;
    [~, params.n, params.J] = size(params.PWMs);
    params.k = k;
    try
        loaded = load('data/temp/G.mat');
        inds = loaded.inds;
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};
    catch
        inds = 1:params.k;
        params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
        params.lengths = params.lengths(inds(end-params.k+1:end));
        params.names = {params.names{inds(end-params.k+1:end)}};
    end
end