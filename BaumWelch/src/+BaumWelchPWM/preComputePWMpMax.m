
% Xs - N x L
% pcPWMp - N x k x L
function pcPWMp = preComputePWMpMax(Xs)
    [PWMs, lengths] = BaumWelchPWM.PWMs();
    [~, n, J] = size(PWMs);
    [N, ~] = size(Xs);
    %0 padding for handling the edges probability calculation
    Xs1H = cat(2, matUtils.mat23Dmat(Xs, n), zeros(N, J, n));
    % NxJxnxk
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    pcPWMp = preComputePWMpAux(PWMsRep, Xs1H, lengths);
end
% calculating the probability of each position in the sequences
% to be emitted by each PWM
% PWMsRep - NxJxnxk
function out = preComputePWMpAux(PWMsRep, Xs1H, lengths)
    % pcPWMp - N x k x L
    persistent pcPWMp
    [N, J, n, k] = size(PWMsRep);
    [~, L, ~] = size(Xs1H);
    L = L - J;
    % G = 0:5:100;
    % g = length(G);
    PC_PWM_PROBABILITY_FILE = fullfile('data', 'precomputation', 'pcPWMpMax.mat');
    if ~isempty(pcPWMp) && length(size(pcPWMp)) == 2 && all(size(pcPWMp) == [N, k])
        % in memory
        out = pcPWMp;
        return;
    end
    if exist(PC_PWM_PROBABILITY_FILE, 'file') == 2
        fprintf('File cache found. Loading pre-computed PWM from %s...\n', PC_PWM_PROBABILITY_FILE);
        load(PC_PWM_PROBABILITY_FILE, 'pcPWMp');
        if length(size(pcPWMp)) == 2 && all(size(pcPWMp) == [N, k])
            % loaded from file, same size
            out = pcPWMp;
            return;
        else
            fprintf('Data size [%d, %d] mismatch [%d, %d], delete and recalculate?\n', size(pcPWMp,1), size(pcPWMp,2), N, k);
            pause
            delete(PC_PWM_PROBABILITY_FILE);
        end
    end
    pcPWMp = calculate(Xs1H, PWMsRep, lengths, N, k, J, n, L, PC_PWM_PROBABILITY_FILE);
    out = pcPWMp;
end

function pcPWMp = calculate(Xs1H, PWMsRep, lengths, N, k, J, n, L, filePath)
    PSSM_THRESHOLD = 0;
    fprintf('Pre-computing PWM probability on %d sequences\n', size(Xs1H, 1));
    size(repmat(1:J, [N, 1, k]))
    size(repmat(permute(lengths, [1,3,2]), [N,J,1]))
    mask = repmat(1:J, [N, 1, k]) > repmat(permute(lengths, [1,3,2]), [N,J,1]);
    mask = permute(mask, [1, 2, 4, 3]);

    pcPWMp = zeros(N, k);
    imask = repmat(1-mask, [1,1,n,1]);
    % N x J x n x k
    bgPWMs = imask .* repmat(permute([0.2573, 0.2433, 0.2426, 0.2568],[1,3,2]), [N, J, 1, k]);
    for t = 1:L
        H1 = BaumWelchPWM.getPWMp(J, PWMsRep, Xs1H, t, mask);
        H0 = BaumWelchPWM.getPWMp(J, bgPWMs, Xs1H, t, mask);
        pssm = log(H1 ./ (H0 + eps));
        pcPWMp = max(pcPWMp, pssm);

        % pcPWMp(:, :, t) = BaumWelchPWM.getPWMp(J, PWMsRep, Xs1H, t, mask);
        fprintf('\r%d / %d', t, L);
    end
    % pcPWMp = cat(3, zeros(N, k, J), pcPWMp);

    save(filePath, 'pcPWMp', '-v7.3');
    fprintf('\nDone.\n');

end