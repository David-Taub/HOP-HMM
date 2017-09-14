
% mainGenSequences(20, 50);
function mergedPeaksMin = mainGenSequences(N, L, m)
    clear pcPWMp
    delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
    params.m = m;
    params.order = 3;
    params.N = N;
    params.L = L;
    params.tEpsilon = 1 / L;
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [params.k, params.n, params.J] = size(PWMs);

    [theta] = BaumWelchPWM.genThetaJ(params);
    theta.F = log(ones(params.m, 1) * 0.02);
    for i=1:params.m
        theta.E(i, :) = [ 9298, 6036, 7032, 4862, 7625, 7735, 6107, 6381, 8470,...
                          1103, 7053, 6677, 4151, 4202, 3203, 4895, 4996, 6551,...
                          4387, 3150, 4592, 7331, 5951, 6900, 6870, 1289, 5710,...
                          6244, 4139, 7391, 4423, 6990, 7396, 9784, 7709, 4184,...
                          1317, 1586, 1363, 1117, 8430, 1545, 7243, 7822, 5593,...
                          9804, 6587, 6132, 5651, 5732, 4213, 3952, 5459, 8393,...
                          7066, 8345, 5533, 1235, 4689, 7744, 5493, 7510, 4973,...
                          9313]';
        % theta.E(i, :) = theta.E(i, :) + (rand(1, params.n ^ params.order) * 3000);
    end
    theta.E = log(bsxfun(@times, theta.E, 1./sum(theta.E, length(size(theta.E)))));

    for i = 1:m
        theta.G(i, :) = matUtils.logMakeDistribution(log(((mod(1:params.k, m) == (i-1)) + eps) .* rand(1, params.k)));
    end
    [seqs, Y] = BaumWelchPWM.genSequencesJ(theta, params);
    overlaps = matUtils.vec2mat(Y(:, 1)', params.m)';
    lengths = ones(params.N, 1) * params.L;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'theta', 'Y');
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.originalTheta = theta;
    mergedPeaksMin.Y = Y;
end