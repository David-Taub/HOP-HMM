
% pcPWMp - N x k x L-1
function pcPWMp = preComputePWMp(Xs)
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [~, n, J] = size(PWMs);
    [N, ~] = size(Xs);
    %0 padding for handling the edges probability calculation
    Xs1H = cat(2, matUtils.mat23Dmat(Xs, n), zeros(N, J, n));
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    pcPWMp = preComputePWMpAux(PWMsRep, Xs1H);
end
% calculating the probability of each position in the sequences
% to be emitted by each PWM
function pcPWMp = preComputePWMpAux(PWMsRep, Xs1H)
    % pcPWMp - N x k x L-1
    [N, J, ~, k] = size(PWMsRep);
    [~, L, ~] = size(Xs1H);
    L = L - J;
    PC_PWM_PROBABILITY_FILE = fullfile('data', 'pcPWMp.mat');
    % persistent p_pcPWMp
    % if isempty(p_pcPWMp)
    try
        fprintf('Loading pre-computed PWM for sequences from %s...\n', PC_PWM_PROBABILITY_FILE);
        load(PC_PWM_PROBABILITY_FILE, 'pcPWMp');
        fprintf('Done.\n');
        % p_pcPWMp = pcPWMp;
    catch
        fprintf('Pre-computing PWM probability on %d sequences', size(Xs1H, 1));
        mask = repmat(1:params.J, [params.N, 1, params.k]) > repmat(permute(params.lengths, [1,3,2]), [N,J,1]);
        pcPWMp = zeros(N, k, L);
        for t = 1:L
            pcPWMp(:, :, t) = BaumWelchPWM.getPWMp(params, PWMsRep, Xs1H, t, mask);
        end
        pcPWMp = cat(3, zeros(N, k, J), pcPWMp);
        fprintf('\n');
        save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', '-v7.3');
        % p_pcPWMp = pcPWMp;
    end
    % end
    % pcPWMp = p_pcPWMp;
end