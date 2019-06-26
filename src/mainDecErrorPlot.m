
function mainDecErrorPlot()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.maxIters = 100;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1000;
    conf.N = 500;
    conf.withExponent = false;

    conf.order = 2;
    conf.m = 6;
    conf.k = 25;
    conf.backgroundAmount = 1;
    conf.doBound = false;
    conf.doResample = false;
    main(conf);
end

function likelihood = calcLikelihood(params, theta, X, pcPWMp)
    N = size(X, 1);
    [~, ~, pX, ~, ~, ~] = EM.EStep(params, theta, X, pcPWMp);
    likelihood = mean(pX, 1);
end


function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;

    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15)

    thetaOrig = trainDataset.theta;
    thetaEst(1) = misc.genTheta(params, true);
    errorsTrain(1) = rateTheta(params, thetaOrig, thetaEst(1), trainDataset);
    errorsTest(1) = rateTheta(params, thetaOrig, thetaEst(1), testDataset);
    trainLikelihood(1) = calcLikelihood(params, thetaEst(1), trainDataset.X, trainDataset.pcPWMp);
    testLikelihood(1) = calcLikelihood(params, thetaEst(1), testDataset.X, testDataset.pcPWMp);
    trainLikelihoodOrig = calcLikelihood(params, thetaOrig, trainDataset.X, trainDataset.pcPWMp)
    testLikelihoodOrig = calcLikelihood(params, thetaOrig, testDataset.X, testDataset.pcPWMp);

    for i = 1:conf.maxIters
        i
        [thetaEst(i + 1), trainLikelihood(i + 1)] = EM.EMIteration(params, trainDataset, thetaEst(i), conf.doBound, conf.doResample);
        thetaEst(i + 1) = misc.permThetaByAnother(params, thetaOrig, thetaEst(i + 1));
        [~, ~, pX, ~, ~, ~] = EM.EStep(params, thetaEst(i + 1), testDataset.X, testDataset.pcPWMp);
        testLikelihood(i + 1) = matUtils.logMatSum(pX, 1);

        errorsTrain(i + 1) = rateTheta(params, thetaOrig, thetaEst(i + 1) , trainDataset);
        errorsTest(i + 1) = rateTheta(params, thetaOrig, thetaEst(i + 1), testDataset);
        testLikelihood(end)
        convergeSpan = floor(0.1 * conf.maxIters);
        if length(testLikelihood) > floor(0.1 * conf.maxIters) & abs(mean(testLikelihood(end - convergeSpan:end), 2) - testLikelihood(end)) < abs(0.005 * testLikelihood(end))
            break
        end
    end
    close all

    show.showTheta(thetaEst(end));
    show.showTheta(thetaOrig);

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on
    title('Viterbi Classification')
    plot([errorsTrain(:).layerError]);
    plot([errorsTest(:).layerError]);
    ylabel('Misclassification Rate')
    xlabel('EM Iteration')
    legend('Train error', 'Test error');
    outpath = sprintf('DecError_VitErr_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    saveas(gcf, outpath);

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    outpath = sprintf('DecError_ThetaOverIter_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    showTwoThetasOverTime(params, thetaOrig, thetaEst, 'false', outpath);
    % show.showTwoThetas(params, thetaOrig, thetaEst, 'false', subtitle, 'tmp.jpg');

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    outpath = sprintf('DecError_MSLEOverIter_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    title(sprintf('Learned \\theta Error Convergence'))
    hold on;
    plot([errorsTrain(:).msle]);
    ylabel('MSLE of learned \theta (lower is better)');
    xlabel('EM Iteration');
    saveas(gcf, outpath);

    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    outpath = sprintf('DecError_Likelihood_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    title('Log Likelihood Convergence')
    hold on;
    plot(trainLikelihood)
    plot(testLikelihood)
    plot([1, length(trainLikelihood)], [trainLikelihoodOrig, trainLikelihoodOrig])
    plot([1, length(testLikelihood)], [testLikelihoodOrig, testLikelihoodOrig])
    legend('Train likelihood of learned \theta', 'Test likelihood of learned \theta', ...
           'Train likelihood of true \theta', 'Test likelihood of true \theta')
    ylabel('Average log likelihood of sequence (higher is better)');
    xlabel('EM Iteration');
    saveas(gcf, outpath);
    keyboard
end


function showTwoThetasOverTime(params, thetaOrig, thetaEsts, withExponent, outpath)
    DOT_SIZE = 20;
    % colors = ['b', 'r', 'g', 'm'];
    hold on;
    BASE_INTENSITY = 0.6;
    minVal = inf;
    maxVal = -inf;
    for j = 1:length(thetaEsts)
        thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEsts(j));
        % m x length(theta(:)) / m
        thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
        thetaEstMat = misc.thetaToMat(params, thetaEst, true);
        % if isempty(strfind(outpath, 'tmp'))
            % fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        % end
        if withExponent
            thetaOrigMat = exp(thetaOrigMat);
            thetaEstMat = exp(thetaEstMat);
        end
        minVal = floor(min(minVal, min([thetaEstMat(:); thetaOrigMat(:)])));
        maxVal = ceil(max(maxVal, max([thetaEstMat(:); thetaOrigMat(:)])));
        % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
        gradColor = BASE_INTENSITY + (0.9 - BASE_INTENSITY) * (length(thetaEsts) - j) / length(thetaEsts);
        % if j > 1
        %     for t = 1: length(thetaOrigMat(:))
        %         plot([thetaOrigMat(t), thetaOrigMat(t)], [prevTheta(t), thetaEstMat(t)], 'color', [gradColor, gradColor, 1]);
        %     end
        % end
        if j == length(thetaEsts)
            gradColor = 0;
        end
        scatter(thetaOrigMat(:), thetaEstMat(:), DOT_SIZE, [gradColor, gradColor, 1], 'filled');
        prevTheta = thetaEstMat;
    end
    msle = mean((log(exp(thetaOrigMat(:)) + 1) - log(exp(thetaEstMat(:)) + 1)) .^ 2, 1);
    plot([minVal, maxVal], [minVal, maxVal]);
    % legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \\theta vs True \\theta Comparison'));
    xlabel('True \theta probabilities');
    ylabel('Estimated \theta probabilities');
    text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), sprintf('Final EM Iteration MSLE = %.2e', msle), 'FontSize', 14)
    saveas(gcf, outpath);
end



function errors = rateTheta(params, thetaOrig, thetaEst, dataset)
    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    YEst = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);
    [errors.totalError, errors.PWMError, errors.layerError] = calcYEstError(params, dataset, YEst);
    errors.msle = calcThetaError(params, thetaOrig, thetaEst);
end


function msle = calcThetaError(params, thetaOrig, thetaEst)
    thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
    thetaEstMat = misc.thetaToMat(params, thetaEst, true);
    msle = mean((log(exp(thetaOrigMat(:)) + 1) - log(exp(thetaEstMat(:)) + 1)) .^ 2, 1);
end


function [totalError, PWMError, layerError] = calcYEstError(params, dataset, YEst)
    [N, L] = size(dataset.X);
    Ymerged = cat(3, dataset.Y, dataset.Y2);
    PWMMask = dataset.Y2(:, :) > 0;
    PWNCount = sum(PWMMask(:), 1);
    layerAccMask = YEst(:, :, 1) == dataset.Y;
    totalAccMask = all(YEst == Ymerged, 3);
    totalError = 1 - (sum(totalAccMask(:)) / (N * L));
    PWMError =  1 - (sum(totalAccMask(:) & PWMMask(:), 1) / PWNCount);
    layerError =  1 - (sum(layerAccMask(:)) / (N * L));
end

% Y - N x L x 2
% mask - N x m x k x L
function mask = genPWMMask(params, Y, N, L)
    mask = [];
    for l = 1:params.k
        layerMask = [];
        for i = 1:params.m
            PWMStateMask = (Y(:, 2:end, 2) == l) & (Y(:, 1:end - 1, 2) == 0) & (Y(:, 1:end - 1, 1) == i);
            % N x L
            PWMStateMask = cat(2, PWMStateMask, false(N, 1, 1));
            layerMask = cat(3, layerMask, PWMStateMask);
        end
        assert(all(size(layerMask) == [N, L, params.m]));
        layerMask = permute(layerMask, [1, 3, 4, 2]);
        mask = cat(3, mask, layerMask);
    end
    assert(all(size(mask) == [N, params.m, params.k, L]));
end
