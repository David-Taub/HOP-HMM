% mainGenSequences(20, 50);
function mergedPeaksMin = mainGenSequences(N, L)
    params.m = 1;
    params.order = 3;
    params.N = N;
    params.L = L;
    params.tEpsilon = 0.0001;
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [params.k, params.n, params.J] = size(PWMs);

    [originalTheta] = BaumWelchPWM.genThetaJ(params);
    originalTheta.G(1,:) = matUtils.logMakeDistribution(log([ones(1,3)/3, eps * ones(1, params.k-3)]));
    [seqs, ~] = BaumWelchPWM.genSequencesJ(originalTheta, params);
    overlaps = ones(params.N, 1);
    lengths = ones(params.N, 1) * params.L;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'originalTheta');
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.originalTheta = originalTheta;

end