
% Xs - N x L
% pcPWMp - N x k x L
function pcPWMp = preComputePWMp(Xs)
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
    PC_PWM_PROBABILITY_FILE = fullfile('data', 'precomputation', 'pcPWMp.mat');
    if ~isempty(pcPWMp) && all(size(pcPWMp) == [N, k, L])
        % in memory
        out = pcPWMp;
        return;
    end
    if exist(PC_PWM_PROBABILITY_FILE, 'file') == 2
        fprintf('Loading pre-computed PWM from %s...\n', PC_PWM_PROBABILITY_FILE);
        load(PC_PWM_PROBABILITY_FILE, 'pcPWMp');
        if all(size(pcPWMp) == [N, k, L])
            fprintf('data size mismatch, recalculating\n');
            % loaded from file, same size
            out = pcPWMp;
            return;
        else
            delete(PC_PWM_PROBABILITY_FILE);
        end
    end
    pcPWMp = calculate(Xs1H, PWMsRep, lengths, N, k, J, n, L);
    fprintf('Saving if file...\n');
    save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', '-v7.3');
    fprintf('Done saving.\n');
    out = pcPWMp;
end

function pcPWMp = calculate(Xs1H, PWMsRep, lengths, N, k, J, n, L)
    fprintf('Pre-computing PWM probability on %d sequences\n', size(Xs1H, 1));
    mask = repmat(1:J, [N, 1, k]) > repmat(permute(lengths, [1,3,2]), [N,J,1]);
    mask = permute(mask, [1, 2, 4, 3]);
    pcPWMp = zeros(N, k, L);
    imask = repmat(1-mask, [1,1,n,1]);
    % N x J x n x k
    bgPWMs = imask .* repmat(permute([0.2573, 0.2433, 0.2426, 0.2568],[1,3,2]), [N, J, 1, k]);
    for t = 1:L
        H1 = BaumWelchPWM.getPWMp(J, PWMsRep, Xs1H, t, mask);
        H0 = BaumWelchPWM.getPWMp(J, bgPWMs, Xs1H, t, mask);
        pcPWMp(:, :, t) = log(H1 ./ (H0 + eps));

        % pcPWMp(:, :, t) = BaumWelchPWM.getPWMp(J, PWMsRep, Xs1H, t, mask);
        fprintf('\r%d / %d', t, L);
    end
    % pcPWMp = cat(3, zeros(N, k, J), pcPWMp);
    fprintf('\nDone calculating\n');
end