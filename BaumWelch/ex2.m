% Try to guess what m was used to create X
% plot the likelihood graph for more and more m
% we should look for the knee in the graph
% also, we want the BIC (Bayesian information criterion) 
% of the esitmation to find the real m
function likelihoodPerM = ex2()
    close all;
    mReal = 6;
    n = 9;
    repetitions = 5;
    ms = 3:10;
    mRealInd = find(ms == mReal);
    epsilon = 0.02;
    k = (1 / epsilon) * mReal * n * 5
    likelihoodPerM = zeros(1, length(ms));
    chi2PerM = zeros(1, length(ms));
    gTestPerM = zeros(1, length(ms));

    for i = 1 : length(ms)

        % errorsReal(end + 1) = calcError(Y, viterbi(X, startTReal, TReal, EReal));
        likelihoodReps = [];
        chi2Reps = [];
        gTestReps = [];
        for t = 1 : repetitions
            [startTReal, TReal, EReal] = genRandEyeParamsExt(mReal, n, epsilon);
            % Create X and Y
            [X, Y] = genHmm(k, startTReal, TReal, EReal);
            [startTMine, TMine, EMine, likelihood] = EM(X, ms(i), n, 400);
            
            % rowsSD(EReal, EMine);

            YEst = viterbi(X, startTMine, TMine, EMine);
            % calcError(Y, YEst)

            % grade from -1 to 1 about how good is YEst
            % likelihood(end + 1) = silhouette(X, YEst, startTMine, TMine, EMine);
            chi2Reps(end + 1) = chi2(YEst, X, n, k, ms(i));
            gTestReps(end + 1) = gTest(YEst, X, n, k, ms(i));
            likelihoodReps(end + 1) = likelihood;

            fprintf('M: %d  repetition: %d Log Likelihood: %f\n', ms(i), t, likelihood)
        end
        likelihoodPerM(i) = mean(likelihoodReps);
        chi2PerM(i) = mean(chi2Reps);
        gTestPerM(i) = mean(gTestReps);
        hold on;
        if i >= mRealInd
            plot(ms(1:i), likelihoodPerM(1:i),'b',[mReal], ...
                likelihoodPerM(mRealInd),'b*', 'erasemode','background');
        else 
            plot(ms(1:i), likelihoodPerM(1:i),'b', 'erasemode','background');
        end

        title('Log(Likelihood) vs #states');
        xlabel('#states');
        ylabel('Log(Likelihood)');
        drawnow;
        hold off;
    end
    paramsAmounts = (ms .^ 2 + ms .* n);
    bic = -2 .* likelihoodPerM + paramsAmounts .* log(k);
    aic = -2 .* likelihoodPerM + 2 .* paramsAmounts;
    figure;
    plot(ms, bic, [mReal], bic(mRealInd), 'b*');
    hold on;
    plot(ms, aic, [mReal], aic(mRealInd), 'b*');
    plot(ms, gTestPerM, [mReal], gTestPerM(mRealInd), 'b*');
    plot(ms, chi2PerM, [mReal], chi2PerM(mRealInd), 'b*');
    title('Models Quality');
    legend('BIC', 'correct model', ...
           'AIC', 'correct model', ...
           'gTest', 'correct model', ...
           'chi2', 'correct model')

    % save('C:\Users\booga\Dropbox\SemesterF\genetics\ex1\x.mat')
end
function [E, O] = chi2Ox(YEst, X, n, k, m)
    O = zeros(n, m);
    Xprob(:, 1) = getProb(X, n);
    E = repmat(Xprob, [1, m]);
    for i = 1 : m
        O(:, i) = getProb(X(YEst == i), n);
    end
end
function ret = gTest(YEst, X, n, k, m)
    [E, O] = chi2Ox(YEst, X, n, k, m);
    relInd = O~=0;
    ret = 2 * sum(O(relInd) .* log(O(relInd) ./ E(relInd)));
end

function ret = chi2(YEst, X, n, k, m)
    [E, O] = chi2Ox(YEst, X, n, k, m);
    ret = ((O - E) .^ 2) ./ E;
    ret = sum(ret(:));
end

function XProb = getProb(X, n)
    k = length(X);
    XProb = sum(vec2mat(X, n), 2) ./ k;
end

function [startT, T, E] = genRandParams(m, n)
    [startT, T, E] = genRandParamsExt(m, n);
    [T, E] = addPrior(startT, T, E);
end

function [T1, E1] = addPrior(startT, T, E)
    [m, n] = size(E);
    E1 = [zeros(1, n); E];
    T1 = [0, startT.'; zeros(m, 1), T];
end



function [startT, T, E] = genRandEyeParamsExt(m, n, epsilon)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + epsilon * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = eye(m, n) + epsilon * rand(m, n);
    E = bsxfun(@times, E, 1 ./ sum(E, 2));
end


function M3d = mat2mat3d(M, n)
    [r,c] = size(M);
    I = repmat(permute([1:n], [1,3,2]), [r, c, 1]);
    M3d = double(I == repmat(M, [1, 1, n]));
end

function sdMat = rowsSD(E1, E2)
    n = size(E1, 1);
    m = size(E2, 1);
    sdMat = zeros(n, m);
    for i=1:n
        for j=1:m
            sdMat(i,j) = norm(E1(i, :) - E2(j, :));
        end
    end
    figure
    subplot(1,3,1);imagesc(sdMat);
    title('Rows Match')
    subplot(1,3,2);imagesc(E1);
    title('Real E')
    subplot(1,3,3);imagesc(E2);
    title('Est. E')
end

