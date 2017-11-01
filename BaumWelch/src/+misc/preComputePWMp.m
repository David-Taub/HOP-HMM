
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
    outOfMotifMask  = repmat(1:J, [N, 1, k]) > repmat(permute(params.lengths, [1,3,2]), [N,J,1]);
    outOfMotifMask = permute(outOfMotifMask, [1, 2, 4, 3]);
    pcPWMp = -inf(N, k, L);
    for startPos = 1:L
        % see c code https://github.com/yatbear/nlangp/blob/master/Baum-Welch/hmm.c
        pcPWMp(:, :, startPos) = PWMLikelihood(J, PWMsRep, Xs1H, startPos, outOfMotifMask);
        if mod(startPos, 20) == 0
            fprintf('\r%d / %d', startPos, L);
        end
    end
    assert not(any(isnan(pcPWMp(:))))
    fprintf('\nDone calculating\n');
end

% res - N x k
% PWMsRep - N x J x n x k
% outOfMotifMask - N x J x 1 x k
% Xs1H - N x L x n
function res = PWMLikelihood(J, PWMsRep, Xs1H, startPos, outOfMotifMask)
    L = size(Xs1H, 2);
    endPos = min(L, startPos+J-1);
    % Note: J2 usually equals J, except when L-startPos is less than J
    J2 = endPos-startPos+1;
    % N x J x n
    lastJXs1H = Xs1H(:, startPos:startPos+J2-1, :);
    k = size(PWMsRep, 4);
    % N x J2 x n x k .* N x J2 x n x k
    res = PWMsRep(:, 1:J2, :, :) .* repmat(lastJXs1H, [1, 1, 1, k]);
    % N x J2 x n x k -> % N x J2 x 1 x k
    res = sum(res, 3);
    % N x J2 x 1 x k
    res = res + outOfMotifMask(:, 1:J2, :, :);
    % N x J2 x 1 x k -> N x 1 x 1 x k
    res = sum(log(res), 2);
    % N x 1 x 1 x k -> N x k
    res = permute(res, [1, 4, 2, 3]);
end