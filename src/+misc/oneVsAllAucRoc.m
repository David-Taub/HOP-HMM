% aucRocs  - m x k
function aucRocs = oneVsAllAucRoc(X, Y, m, k, PWMs, lengths)
    expected_num_of_peaks_in_seq = 2;
    n = 4;
    aucRocs = zeros(m, k);
    Xs1H = matUtils.mat23Dmat(X, n);
    for pwmId = 1:k
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        for tissueID = 1:m
            tissueMask = Y(:,1) == tissueID;
            % N x L
            pos = misc.maxN(PWMLogLike(tissueMask, :), 2, expected_num_of_peaks_in_seq);
            neg = misc.maxN(PWMLogLike(~tissueMask, :), 2, expected_num_of_peaks_in_seq);
            aucRocs(tissueID, pwmId) = matUtils.getAucRoc(pos(:), neg(:), false, true);
            [v, i] = max(aucRocs(:));
            fprintf('Best 1vsAll AucRocs for tissue %d of PWM %d is %.2f / %.2f (%d)\n', tissueID, pwmId, aucRocs(tissueID, pwmId), v, i);
        end
    end
end