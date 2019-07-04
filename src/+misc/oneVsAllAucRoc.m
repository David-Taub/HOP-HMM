% aucRocs  - m x k
function aucRocs = oneVsAllAucRoc(X, Y, PWMs, lengths)
    EXPECTED_NUM_OF_PEAKS_IN_SEQ = 2;
    n = 4;
    m = max(Y(:));
    k = size(PWMs, 1);
    aucRocs = zeros(m, k);
    Xs1H = matUtils.mat23Dmat(X, n);
    for pwmId = 1:k
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        for tissueID = 1:m
            tissueMask = Y(:, 1) == tissueID;
            % N x L
            pos = misc.maxN(PWMLogLike(tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
            neg = misc.maxN(PWMLogLike(~tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
            aucRocs(tissueID, pwmId) = matUtils.getAucRoc(pos(:), neg(:), false, true);
            [v, i] = max(aucRocs(:));
            fprintf('Best 1vsAll AucRocs for tissue %d / %d of PWM %d / %d is %.2f / %.2f (%d)\n', tissueID, m, ...
                    pwmId, k, aucRocs(tissueID, pwmId), v, i);
        end
    end
end