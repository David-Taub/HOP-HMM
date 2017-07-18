
% Xs - N x L
% pcPWMp - N x k x L
function pcPWMp = preComputePWMpMax(Xs)
    [PWMs, lengths] = BaumWelchPWM.PWMs();
    [~, n, J] = size(PWMs);
    [N, L] = size(Xs);
    %0 padding for handling the edges probability calculation
    % PWMsRep - J x L-J+1 x n x k
    Xs1H = matUtils.mat23Dmat(Xs, n);
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, L-J+1]), [3, 4, 2, 1]);
    pcPWMp = preComputePWMpAux(PWMsRep, Xs1H, lengths);
end
% calculating the probability of each position in the sequences
% to be emitted by each PWM
% PWMsRep - NxJxnxk
function out = preComputePWMpAux(PWMsRep, Xs1H, lengths)
    % pcPWMp - N x k x L
    persistent pcPWMp
    [J, ~, n, k] = size(PWMsRep);
    [N, L, ~] = size(Xs1H);
    % G = 0:5:100;
    % g = length(G);
    PC_PWM_PROBABILITY_FILE = fullfile('data', 'precomputation', 'pcPWMpMax.mat');
    if ~isempty(pcPWMp) && length(size(pcPWMp)) == 2 && all(size(pcPWMp) == [N, k])
        % in memory
        fprintf('Used memory cache\n');
        out = pcPWMp;
        return;
    end
    if exist(PC_PWM_PROBABILITY_FILE, 'file') == 2
        fprintf('File cache found. Loading pre-computed PWM from %s...\n', PC_PWM_PROBABILITY_FILE);
        load(PC_PWM_PROBABILITY_FILE, 'pcPWMp');
        if length(size(pcPWMp)) == 2 && all(size(pcPWMp) == [N, k])
            % loaded from file, same size
            fprintf('Used file cache\n');
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
% Xs1H - N x L x n
% PWMsRep - J x L-J+1 x n x k
% lengths - 1 x k
function pcPWMp = calculate(Xs1H, PWMsRep, lengths, N, k, J, n, L, filePath)
    % PSSM_THRESHOLD = 0;
    fprintf('Pre-computing PWM probability on %d sequences\n', size(Xs1H, 1));
    % mask - J x L-J+1 x 1 x k
    mask = repmat([1:J]', [1, L-J+1, 1, k]) > repmat(permute(lengths, [4,3,1,2]), [J,L-J+1,1,1]);
    pcPWMp = zeros(N, 4 * k);
    imask = repmat(1-mask, [1,1,n,1]);
    % N x J x n x k
    randBaseDist = permute(sum(sum(Xs1H, 2), 1) ./ (N * L), [1,3,2]);
    % randBaseDist = [0.2573, 0.2433, 0.2426, 0.2568]
    randBaseDist
    T1 = 3;
    T2 = 5;
    T3 = 10;
    bgPWMs = imask .* repmat(permute(randBaseDist,[1,3,2]), [J, L-J+1, 1, k]);
    for t = 1:N
        % L-J+1 x k
        H1 = BaumWelchPWM.getPWMpMax(J, PWMsRep, Xs1H(t, :,:), t, mask);
        H0 = BaumWelchPWM.getPWMpMax(J, bgPWMs, Xs1H(t, :,:), t, mask);
        pssm = log(H1 ./ (H0 + eps));
        [pcPWMp(t, 1:k), ~] = max(pssm, [], 1);
        % [pcPWMp(t, 1:k), pcPWMp(t, k+1:2*k)] = max(pssm, [], 1);
        pssm(pssm<0) = 0;
        pssmSorted = sort(pssm, 1);
        pcPWMp(t, 1*k+1:2*k) = sum(pssmSorted(end-T1+1:end, :), 1);
        pcPWMp(t, 2*k+1:3*k) = sum(pssmSorted(end-T2+1:end, :), 1);
        pcPWMp(t, 3*k+1:4*k) = sum(pssmSorted(end-T3+1:end, :), 1);

        % pcPWMp(:, :, t) = BaumWelchPWM.getPWMp(J, PWMsRep, Xs1H, t, mask);
        fprintf('%d / %d ', t, N);
        fprintf('%d', pssm(1:100) > 0);
        fprintf('\n');
    end
    % pcPWMp = cat(3, zeros(N, k, J), pcPWMp);

    save(filePath, 'pcPWMp', '-v7.3');
    fprintf('\nDone.\n');

end