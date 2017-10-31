
% Xs - N x L
% pcPWMp - N x k x L
function pcPWMp = preComputePWMp(Xs, params)
    %0 padding for handling the edges probability calculation
    Xs1H = matUtils.mat23Dmat(Xs, params.n);
    % N x J x n x k
    pcPWMp = preComputePWMpAux(Xs1H, params);
end

% calculating the probability of each position in the sequences
% to be emitted by each PWM
% PWMsRep - N x J x n x k
function out = preComputePWMpAux(Xs1H, params)
    % pcPWMp - N x k x L
    persistent pcPWMp
    persistent sample
    [k, ~, J] = size(params.PWMs);
    [N, L, ~] = size(Xs1H);
    % L = L - fJ;
    PC_PWM_PROBABILITY_FILE = fullfile('data', 'precomputation', 'pcPWMp.mat');
    newSample = Xs1H(1:400);
    if ~isempty(pcPWMp) && all(newSample == sample)
        % in memory
        fprintf('Pre-computed PWM probability from memory cache\n');
        out = pcPWMp;
        return;
    end
    try
        if exist(PC_PWM_PROBABILITY_FILE, 'file') == 2
            fprintf('Loading pre-computed PWM from %s...\n', PC_PWM_PROBABILITY_FILE);
            load(PC_PWM_PROBABILITY_FILE, 'pcPWMp');
            if all(newSample == sample)
                fprintf('data match\n');
                % loaded from file, same size
                out = pcPWMp;
                return;
            else
                delete(PC_PWM_PROBABILITY_FILE);
                fprintf('data file mismatch\n');
            end
        end
    catch ME
        fprintf('File not found\n');
        ME
    end
    fprintf('Calculating pre-computed PWM probability\n');
    pcPWMp = calculate(Xs1H, N, k, J, L, params);
    sample = newSample;
    assert(not(any(isnan(pcPWMp(:)))))
    fprintf('Saving data in file %s...\n', PC_PWM_PROBABILITY_FILE);
    save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', 'sample', '-v7.3');
    fprintf('Done saving.\n');

    out = pcPWMp;
end

% PWMsRep - N x J x n x k
% Xs1H - N x L + J x n
% params.lengths - k x 1
% pcPWMp - N x k x L
function pcPWMp = calculate(Xs1H, N, k, J, L, params)
    fprintf('Pre-computing PWM probability on %d sequences\n', size(Xs1H, 1));
    % [params.PWMs, params.lengths] = misc.params.PWMs();
    PWMsRep = permute(repmat(params.PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    mask = repmat(1:J, [N, 1, k]) > repmat(permute(params.lengths, [1,3,2]), [N,J,1]);
    mask = permute(mask, [1, 2, 4, 3]);
    pcPWMp = -inf(N, k, L);
    for t = 1:L
        % see c code https://github.com/yatbear/nlangp/blob/master/Baum-Welch/hmm.c
        pcPWMp(:, :, t) = misc.getPWMp(J, PWMsRep, Xs1H, t, mask);
        if mod(t, 20) == 0
            fprintf('\r%d / %d', t, L);
        end
        if any(isnan(pcPWMp(:)))
            keyboard
        end
    end
    fprintf('\nDone calculating\n');
end