% mainGenSequences(20, 50);
function mergedPeaksMin = mainGenSequences(N, L, m)
    params.m = m;
    params.order = 3;
    params.N = N;
    params.L = L;
    params.tEpsilon = 0.0001;
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [params.k, params.n, params.J] = size(PWMs);

    [originalTheta] = BaumWelchPWM.genThetaJ(params);
    for i = 1:m
        originalTheta.G(i, :) = matUtils.logMakeDistribution(log([ones(1, 3)/3, eps * ones(1, params.k-3)]));
        % originalTheta.G(i, :) = originalTheta.G(i, randperm(params.k));
    end
    [seqs, Y] = BaumWelchPWM.genSequencesJ(originalTheta, params);
    overlaps = matUtils.vec2mat(Y(:, 1)', params.m)';
    lengths = ones(params.N, 1) * params.L;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'originalTheta');
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.originalTheta = originalTheta;

end