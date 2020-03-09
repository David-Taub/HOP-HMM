% Xs1H - N x L x n
% Y - N x L
% PWM - k x n x J
% PWM - 1 x J
% aucRocs  - m x k
function aucRocs = oneVsAllAucRoc(Xs1H, Y, PWMs, lengths)
    EXPECTED_NUM_OF_PEAKS_IN_SEQ = 3;
    [k, n, J] = size(PWMs);
    m = max(Y(:));
    aucRocs = zeros(m, k);
    % N x L x n
    for pwmId = 1:k
        % N x L
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        for tissueID = 1:m
            tissueMask = Y(:, 1) == tissueID;
            % mask x EXPECTED_NUM_OF_PEAKS_IN_SEQ
            pos = misc.maxN(PWMLogLike(tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
            neg = misc.maxN(PWMLogLike(~tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
            aucRocs(tissueID, pwmId) = matUtils.getAucRoc(pos(:), neg(:), false, true);
            [v, i] = max(aucRocs(:));
            fprintf('Best 1vsAll AucRocs for tissue %d / %d of PWM %d / %d is %.2f / %.2f (%d)\n', tissueID, m, ...
                    pwmId, k, aucRocs(tissueID, pwmId), v, i);
        end
    end
end