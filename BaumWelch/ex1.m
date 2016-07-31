% tasks:

function ex1()
    close all;
    m = 2;
    n = 3;
    k = 5000;
    % Create real probaiblities
    [startTReal, TReal, EReal] = genRandEyeParamsExt(m, n);

    % Create X and Y
    [X, Y] = genHmm(k, startTReal, TReal, EReal);
    % [X, Y] = hmmgenerate(k, TReal, EReal)
    
    [~, TGuess, EGuess] = genRandParams(m, n);
    [TMatlab, EMatlab] = hmmtrain(X, TGuess, EGuess);
    [startTMine, TMine, EMine] = EM(X', m, n, 400);
    [TMineP, EMineP] = addPrior(startTMine, TMine, EMine);

    YEst1 = hmmviterbi(X, TReal, EReal);
    YEst2 = hmmviterbi(X, TMatlab, EMatlab);
    YEst3 = hmmviterbi(X, TMineP, EMineP);
    YEst4 = viterbi(X, startTMine, TMine, EMine);
    errors = [  calcError(Y, YEst1);...
                calcError(Y, YEst2);...
                calcError(Y, YEst3);...
                calcError(Y, YEst4)...
             ]
end

function e = calcError(Y, YEst)
    ct = crosstab(Y, YEst);
    correct = max(ct, [], 2);
    e = 1 - sum(correct) / length(Y);
end



function [startT, T, E] = genRandEyeParamsExt(m, n)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    epsilon = 0.03;
    T = eye(m) + epsilon * rand(m);
    % T = rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = rand(m, n);
    E = bsxfun(@times, E, 1 ./ sum(E, 2));
end