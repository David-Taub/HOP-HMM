% Try to guess what m was used to create X
% plot the likelihood graph for more and more m
% we should look for the knee in the graph
% also, we want the BIC (Bayesian information criterion) 
% of the esitmation to find the real m
function likelihoodPerM = ex3()
    close all;
    mReal = 5;
    n = 6;
    repetitions = 5;
    ms = 3:7;
    mRealInd = find(ms == mReal);
    epsilonE = 0.01;
    epsilonT = 0.05;
    lenFactor = 100; % around 100 yeilds good results
    k = round(1 / mean([epsilonT, epsilonE]) * mReal * lenFactor);
    fprintf('Samples length: %f\n', k);
    dataPerM = [];


    for t = 1 : repetitions
        % errorsReal(end + 1) = calcError(Y, viterbi(X, startTReal, TReal, EReal));
        for i = 1 : length(ms)
            [startTReal, TReal, EReal] = genRandEyeParamsExt(mReal, n, epsilonT, epsilonE);
            % Create X and Y
            [X, Y] = genHmm(k, startTReal, TReal, EReal);
            [startTMine, TMine, EMine, likelihood] = EM(X, ms(i), n, 400);
            if t == 1
                rowsSD(EReal, EMine);
            end

            YEst = viterbi(X, startTMine, TMine, EMine);
            % calcError(Y, YEst)

            % grade from -1 to 1 about how good is YEst
            % likelihood(end + 1) = silhouette(X, YEst, startTMine, TMine, EMine);
            dataPerM(i, t, 1) = likelihood;
            dataPerM(i, t, 2) = chi2(YEst, X, n, k, ms(i));
            dataPerM(i, t, 3) = gTest(YEst, X, n, k, ms(i));

            fprintf('M: %d  repetition: %d Log Likelihood: %f\n', ms(i), t, likelihood)
        end
        
    end

    likelihoodPerM(1,:) = mean(dataPerM(:, :, 1), 2)
    chi2PerM(1,:) = mean(dataPerM(:, :, 2), 2)
    gTestPerM(1,:) = mean(dataPerM(:, :, 3), 2)

    % log likelihood plot
    figure;
    plot(ms(1:i), likelihoodPerM(1:i),'b',[mReal], ...
        likelihoodPerM(mRealInd),'b*');
    title('Log(Likelihood) vs #states');
    xlabel('#states');
    ylabel('Log(Likelihood)');

    % model grades
    BITS_PER_NUMBER = 2^5;
    paramsAmounts = (ms .* (ms - 1) + ms .* (n - 1) + (ms - 1)) .* BITS_PER_NUMBER;
    bic = -2 .* likelihoodPerM + paramsAmounts .* log(k);
    aic = -2 .* likelihoodPerM + 2 .* paramsAmounts;
    figure;
    plot(ms, bic);
    hold on;
    plot(ms, aic);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality I');
    legend('BIC', 'AIC')
    
    figure;
    plot(ms, chi2PerM);
    hold on;
    plot(ms, gTestPerM);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality II');
    legend('chi2', 'gTest')

    LRT = -2.*(likelihoodPerM(1:end-1) - likelihoodPerM(2:end))
    figure;
    plot(ms(2:end), LRT);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality III');

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
    relInd = E~=0;
    ret = ((O(relInd) - E(relInd)) .^ 2) ./ E(relInd);
    ret = sum(ret(:));
end

function XProb = getProb(X, n)
    k = length(X);
    XProb = sum(vec2mat(X, n), 2) ./ k;
    if k==0
        XProb = ones(n, 1) ./ n;
    end
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



function [startT, T, E] = genRandEyeParamsExt(m, n, epsilonT, epsilonE)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + epsilonT * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = eye(m, n) + epsilonE * rand(m, n);
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
    drawnow
end

