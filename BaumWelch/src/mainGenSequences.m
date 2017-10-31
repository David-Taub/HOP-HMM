
% mainGenSequences(3000, 500, 5, false);
function mergedPeaksMin = mainGenSequences(N, L, m, isMixed)
    clear pcPWMp
    delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
    params = misc.genParams();
    PWM_NOISE_FACTOR = 0.0;
    params.m = m;
    params.N = N;
    params.L = L;
    params.tEpsilon = isMixed * (1 / L);

    originalTheta = misc.genThetaJ(params);
    for i=1:params.m
        originalTheta.E(i, :) = [ 9298, 6036, 7032, 4862, 7625, 7735, 6107, 6381, 8470,...
                                  1103, 7053, 6677, 4151, 4202, 3203, 4895, 4996, 6551,...
                                  4387, 3150, 4592, 7331, 5951, 6900, 6870, 1289, 5710,...
                                  6244, 4139, 7391, 4423, 6990, 7396, 9784, 7709, 4184,...
                                  1317, 1586, 1363, 1117, 8430, 1545, 7243, 7822, 5593,...
                                  9804, 6587, 6132, 5651, 5732, 4213, 3952, 5459, 8393,...
                                  7066, 8345, 5533, 1235, 4689, 7744, 5493, 7510, 4973,...
                                  9313]';
        % originalTheta.E(i, :) = originalTheta.E(i, :) + (rand(1, params.n ^ params.order) * 3000);
    end
    originalTheta.E = log(bsxfun(@times, originalTheta.E, 1./sum(originalTheta.E, length(size(originalTheta.E)))));
    F = sum(exp(originalTheta.G), 2);
    for i = 1:m
        originalTheta.G(i, :) = matUtils.logMakeDistribution(log(((mod(1:params.k, m) == (i-1)) + eps) .* rand(1, params.k)));
        originalTheta.G(i, :) = matUtils.logMakeDistribution(log(exp(originalTheta.G(i, :)) + (rand(1, params.k) .* PWM_NOISE_FACTOR)));
    end
    originalTheta.G = originalTheta.G + repmat(log(F), [1, params.k]);
    [seqs, Y] = misc.genSequencesJ(originalTheta, params);
    Y2 = Y;
    Y = mode(Y(:,:,1), 2);
    overlaps = matUtils.vec2mat(Y(:, 1)', params.m)';
    lengths = ones(params.N, 1) * params.L;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'originalTheta', 'Y', 'originalTheta');
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.originalTheta = originalTheta;
    mergedPeaksMin.Y = Y
    mergedPeaksMin.Y2 = Y2;
end

