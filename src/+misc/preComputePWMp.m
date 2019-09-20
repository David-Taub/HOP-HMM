
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
    filename = misc.pathMaker(params, N, L, 0, 'pcPWMp', '.mat');
    pcPWMProbabilityFile = fullfile('..', 'data', 'precomputation', filename);

    newSample = [Xs1H(1:500), Xs1H(end-499:end), params.k, L, N];
    try
    if ~isempty(pcPWMp) && all(newSample == sample)
        % in memory
        fprintf('Pre-computed PWM probability from memory cache\n');
        out = pcPWMp;
        return;
    end
        if exist(pcPWMProbabilityFile, 'file') == 2
            fprintf('Loading pre-computed PWM from %s...\n', pcPWMProbabilityFile);
            load(pcPWMProbabilityFile, 'pcPWMp');
            if all(newSample == sample)
                fprintf('data match\n');
                % loaded from file, same size
                out = pcPWMp;
                return;
            else
                delete(pcPWMProbabilityFile);
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
    assert(not(any(isnan(pcPWMp(:)))));

    [basedir, ~, ~] = fileparts(pcPWMProbabilityFile);
    if ~exist(basedir, 'dir')
       mkdir(basedir);
    end

    save(pcPWMProbabilityFile, 'pcPWMp', 'sample', '-v7.3');
    fprintf('Saved data in %s\n', pcPWMProbabilityFile);

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
        newVec = misc.PWMLogLikelihood(params.PWMs, params.lengths, Xs1H, pwmId);
        fprintf('\r%d / %d', pwmId, params.k);
        pcPWMp(1:N, pwmId, 1:L) = permute(newVec, [1, 3, 2]);
        fprintf('..');
        assert(not(any(isnan(newVec(:)))));
    end
    fprintf('\nDone calculating\n');
end
