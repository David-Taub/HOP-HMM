% mainGenSequences(20, 50);
function [X, Y] = mainGenSequences(N, L)
    params.m = 1;
    params.order = 3;
    params.N = N;
    params.L = L;
    params.tEpsilon = 0.0001;
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [params.k, params.n, params.J] = size(PWMs);

    [originalTheta] = BaumWelchPWM.genThetaJ(params);
    originalTheta.G(1,:) = log([ones(1,3)/3, zeros(1, params.k-3)]);
    [X, Y] = BaumWelchPWM.genSequencesJ(originalTheta, params);

    overlaps = ones(params.N, 1);
    lengths = ones(params.N, 1) * params.L;
    seqs = X;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'originalTheta');
end