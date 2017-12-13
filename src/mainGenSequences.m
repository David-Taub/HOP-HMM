
% mainGenSequences(100, 10000, 5, 25, true, true);
function mergedPeaksMin = mainGenSequences(N, L, m, k, isMixed, withBackground)
    dbstop if error
    clear pcPWMp
    delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
    params = misc.genParams(k);
    params.m = m;
    params.tEpsilon = isMixed * params.tEpsilon;
    theta = genHumanTheta(params, withBackground);

    [seqs, Y] = misc.genSequences(theta, params, N, L);
    Y2 = Y;
    Y = mode(Y(:,:,1), 2);
    overlaps = matUtils.vec2mat(Y(:, 1)', params.m)';
    lengths = ones(N, 1) * L;
    save(fullfile('data', 'dummyDNA.mat'), 'seqs', 'lengths', 'overlaps', 'theta', 'Y', 'Y2', 'theta');
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.theta = theta;
    mergedPeaksMin.Y = Y
    mergedPeaksMin.Y2 = Y2;
end

function theta = genHumanTheta(params, withBackground)
    PWM_NOISE_FACTOR = 0.0;
    theta = misc.genTheta(params);
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
    PTotalBaseToSub = params.enhancerMotifsRatio / (mean(params.lengths, 2) * (1-params.enhancerMotifsRatio));
    G = zeros(params.m, params.k);
    for i = 1:params.m
        G(i, :) = ((mod(1:params.k, params.m) == (i-1)) + eps) .* rand(1, params.k);
    end
    G = G ./ repmat(sum(G, 2), [1, params.k]);
    G = G .* PTotalBaseToSub;
    PTotalBaseToOtherBase = 1/((1-params.enhancerMotifsRatio)*params.enhancerLength);
    PStayInBase = 1-PTotalBaseToSub-PTotalBaseToOtherBase;
    PCrossBaseState = (1-PStayInBase-PTotalBaseToSub) ./ (params.m-1);
    T = ones(params.m) * PCrossBaseState;
    T(eye(params.m) == 1) = PStayInBase;
    if withBackground
        theta.G(end, :) = log(eps);
        PBackgroundToEnhancer = (1-params.backgroundRatio)/(params.backgroundRatio*params.enhancerLength*(params.m-1));
        PCrossBaseState = params.crossEnhancer*(1-PStayInBase-PTotalBaseToSub)/(params.m-2);
        PBaseToBackground = (1-params.crossEnhancer) * (1-PStayInBase-PTotalBaseToSub);
        T = ones(params.m) * PCrossBaseState;
        T(eye(params.m) == 1) = PStayInBase;
        T(:, end) = PBaseToBackground;
        T(end, :) = PBackgroundToEnhancer;
        T(end, end) = 1 - (sum(T(end, 1:end-1), 2) + sum(G(end, :), 2));
    end
    theta.T = log(T);
    theta.G = log(G);

end
