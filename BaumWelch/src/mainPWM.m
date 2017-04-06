% pcPWMp = mainPWM()
function [pcPWMp, XTrain] = mainPWM(pcPWMpIn, XTrainIn)
    % get the n'st most frequent overlap
    L = 500;
    n = [1:1];
    if nargin < 1
        [posSeqs, negSeqs] = loadTommySeqs(L);

        trainLabLength = ceil(size(posSeqs,1) / 2);

        XTrain = [posSeqs(1:trainLabLength, :); negSeqs(1:trainLabLength, :)];
        XTest = [posSeqs(trainLabLength + 1: end, :); negSeqs(trainLabLength + 1:end, :)];
        YTrain = cat(1, ones(trainLabLength, 1), ones(trainLabLength, 1) .* 2);
        YTest = cat(1, ones(size(posSeqs, 1) - trainLabLength,1), ones(size(posSeqs, 1) - trainLabLength,1) .* 2);
        pcPWMp = preComputePWMpAux(XTrain);
    else
        pcPWMp = pcPWMpIn;
        XTrain = XTrainIn;

    end
    learn(XTrain, pcPWMp);


end

% L - sequence lengths
% means unique overlap between tissues, and if n is 1:3 then we take
% the sequences of the three most frequent class
function [posSeqs, negSeqs] = loadTommySeqs(L)
    fprintf('Loading sequences of Tommy...\n');
    negSeqsTrain = matUtils.readSeq(fullfile('data', 'NEnhancers.train.seq'), L);
    posSeqsTrain = matUtils.readSeq(fullfile('data', 'Enhancers.train.seq'), L);
    negSeqsTest = matUtils.readSeq(fullfile('data', 'NEnhancers.test.seq'), L);
    posSeqsTest = matUtils.readSeq(fullfile('data', 'Enhancers.test.seq'), L);

    posSeqs = [posSeqsTest; posSeqsTrain];
    negSeqs = [negSeqsTest; negSeqsTrain];

    N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:N, :);
    negSeqs = negSeqs(1:N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);
end

function [startT, T, Y, E, F, likelihood, gamma] = learn(XTrain, pcPWMp)
    order = 3; m = 10; n = 4; maxIter = 20; tEpsilon = 0.01;
    [PWMs, lengths] = BaumWelchPWM.PWMs();
    [startT, T, Y, E, F, likelihood, gamma] = BaumWelchPWM.EMJ(XTrain, m, maxIter, tEpsilon, order, PWMs, lengths, pcPWMp);
end


function pcPWMp = preComputePWMpAux(Xs)
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [~, n, J] = size(PWMs);
    [N, ~] = size(Xs);
    Xs1H = cat(2, zeros(N, J, n), matUtils.mat23Dmat(Xs, n));
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    pcPWMp = preComputePWMp(PWMsRep, Xs1H);
end

function pcPWMp = preComputePWMp(PWMsRep, Xs1H)
    % pcPWMc - N x k x L-1
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