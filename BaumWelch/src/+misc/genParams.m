
function params = genParams()
    params.m = 1;
    params.order = 3;
    [params.PWMs, params.lengths, params.names] = misc.PWMs();
    params.tEpsilon = 1/500;
    params.batchSize = 2;
    [~, params.n, params.J] = size(params.PWMs);
    loaded = load('data/temp/G.mat');
    inds = loaded.inds;
    params.k = 50;
    params.PWMs = params.PWMs(inds(end - params.k + 1 : end), :, :);
    params.lengths = params.lengths(inds(end-params.k+1:end));
    params.names = {params.names{inds(end-params.k+1:end)}};
end