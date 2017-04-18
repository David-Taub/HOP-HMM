
% pcPWMp - N x k x L-1
function pcPWMp = preComputePWMp(Xs)
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [~, n, J] = size(PWMs);
    [N, ~] = size(Xs);
    Xs1H = cat(2, zeros(N, J, n), matUtils.mat23Dmat(Xs, n));
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    pcPWMp = preComputePWMpAux(PWMsRep, Xs1H);
end

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
        pcPWMp = zeros(N, k, L-1);
        for t = 2:L

            pcPWMp(:, :, t - 1) = BaumWelchPWM.getPWMp(PWMsRep, Xs1H, t);
            sum(pcPWMp(:, :, t - 1))
        end
        pcPWMp = cat(3, zeros(N, k, J), pcPWMp);
        fprintf('\n');
        save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', '-v7.3');
        % p_pcPWMp = pcPWMp;
    end
    % end
    % pcPWMp = p_pcPWMp;
end