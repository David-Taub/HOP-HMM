
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
function out = preComputePWMpAux(Xs1H, params)
    % pcPWMp - N x k x L
    persistent pcPWMp
    persistent sample
    [N, L, ~] = size(Xs1H);
    % L = L - fJ;
    % TODO: pcPWMp is too big for the memory. Tried matfile and memfile,
    % no good, too slow and unfit for this size. Tall array should work.
    PC_PWM_PROBABILITY_FILE = fullfile('..', 'data', 'precomputation', 'pcPWMp.mat');

    newSample = [Xs1H(1:500), Xs1H(end-499:end), params.k, L, N];
    try
    if ~isempty(pcPWMp) && all(newSample == sample)
        % in memory
        fprintf('Pre-computed PWM probability from memory cache\n');
        out = pcPWMp;
        return;
    end
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
    pcPWMp = calculate(params, Xs1H);
    sample = newSample;
    assert(not(any(isnan(pcPWMp(:)))))
    fprintf('Saving data in file %s...\n', PC_PWM_PROBABILITY_FILE);
    save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', 'sample', '-v7.3');
    fprintf('Done saving.\n');

    out = pcPWMp;
end

% Xs1H - N x L x n
% params.lengths - k x 1
% pcPWMp - N x k x L
function pcPWMp = calculate(params, Xs1H)
    [N, L, ~] = size(Xs1H);
    fprintf('Pre-computing PWM probability on %d sequences\n', size(Xs1H, 1));
    for pwmId = 1:params.k
        % N x L
        newVec = misc.PWMLogLikelihood(params, Xs1H, pwmId);
        fprintf('\r%d / %d', pwmId, params.k);
        pcPWMp(1:N, pwmId, 1:L) = permute(newVec, [1, 3, 2]);
        fprintf('..');
        assert(not(any(isnan(newVec(:)))));
    end
    fprintf('\nDone calculating\n');
end
