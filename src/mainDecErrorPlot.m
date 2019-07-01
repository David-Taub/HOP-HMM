% TODO: BPI bit per instance
% TODO: add tutorial to the post graph
% TODO: add note about areas where the posterior confused between different Enhancers)
% TODO: single permutation
function mainDecErrorPlot()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.maxIters = 150;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1000;
    conf.N = 1000;
    conf.withExponent = false;

    conf.order = 2;
    conf.m = 6;
    conf.k = 25;
    conf.backgroundAmount = 1;
    conf.doGTBound = false;
    conf.doResampling = false;
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
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing, conf.doGTBound, conf.doResampling);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;

    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15)

    thetaOrig = trainDataset.theta;
    thetaEsts(1) = misc.genTheta(params, true);
    % errorsTrain(1) = rateTheta(params, thetaOrig, thetaEsts(1), trainDataset);
    % errorsTest(1) = rateTheta(params, thetaOrig, thetaEsts(1), testDataset);
    trainLikelihood(1) = calcLikelihood(params, thetaEsts(1), trainDataset.X, trainDataset.pcPWMp);
    testLikelihood(1) = calcLikelihood(params, thetaEsts(1), testDataset.X, testDataset.pcPWMp);
    trainLikelihoodOrig = calcLikelihood(params, thetaOrig, trainDataset.X, trainDataset.pcPWMp);
    testLikelihoodOrig = calcLikelihood(params, thetaOrig, testDataset.X, testDataset.pcPWMp);

    for i = 2:conf.maxIters
        i
        [thetaEsts(i), trainLikelihood(i)] = EM.EMIteration(params, trainDataset, thetaEsts(i - 1),...
                                                                   conf.doGTBound);
        % thetaEsts(i) = misc.permThetaByAnother(params, thetaOrig, thetaEsts(i));
        testLikelihood(i) = calcLikelihood(params, thetaEsts(i), testDataset.X, testDataset.pcPWMp);
        trainLikelihood(i) = calcLikelihood(params, thetaEsts(i), trainDataset.X, trainDataset.pcPWMp);

        % errorsTrain(i) = rateTheta(params, thetaOrig, thetaEsts(i) , trainDataset);
        % errorsTest(i) = rateTheta(params, thetaOrig, thetaEsts(i), testDataset);
        testLikelihood(end)
        convergeSpan = floor(0.1 * conf.maxIters);
        if length(trainLikelihood) > floor(0.1 * conf.maxIters) & ...
            abs(mean(trainLikelihood(end - convergeSpan:end), 2) - trainLikelihood(end)) < abs(0.001 * trainLikelihood(end))
            break
        end
    end
    close all

    show.showTheta(thetaEsts(end));
    show.showTheta(thetaOrig);
    origTrainError = rateTheta(params, thetaOrig, thetaOrig, trainDataset);
    origTestError = rateTheta(params, thetaOrig, thetaOrig, testDataset);
    % perm = matUtils.repermuteMat(vectorizedOrig, vectorizedEst);
    perm = misc.munkres(pdist2(misc.thetaToMat(params, thetaOrig, false), misc.thetaToMat(params, thetaEsts(end), false)))';
    for i = 1:length(thetaEsts)
        thetaEst = misc.permTheta(thetaEsts(i), perm);
        errorsTrain(i) = rateTheta(params, thetaOrig, thetaEst , trainDataset);
        errorsTest(i) = rateTheta(params, thetaOrig, thetaEst, testDataset);
    end
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    hold on;
    title('Viterbi Classification Convergence');
    plot([errorsTrain(:).layerError] .* 100, 'LineWidth', 2);
    plot([errorsTest(:).layerError] .* 100, 'LineWidth', 2);
    plot([1, length([errorsTrain(:).layerError])], [origTrainError.layerError, origTrainError.layerError] .* 100)
    plot([1, length([errorsTrain(:).layerError])], [origTestError.layerError, origTestError.layerError] .* 100)
    hold off;
    ylabel('Misclassification Rate (lower is better)');
    ytickformat('percentage');
    xlabel('EM Iteration');
    legend('Viterbi train error with learned \theta', 'Viterbi test error with learned \theta', ...
           'Viterbi train error with true \theta', 'Viterbi test error with true \theta');
    xlim([1, length([errorsTrain(:).layerError])]);
    outpath = misc.pathMaker(params, conf.N, conf.L, 'DecError_VitErr', '.jpg');
    saveas(gcf, outpath);

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    outpath = misc.pathMaker(params, conf.N, conf.L, 'DecError_ThetaOverIter', '.jpg');
    showTwoThetasOverTime(params, thetaOrig, thetaEsts, 'false', outpath);
    % show.showTwoThetas(params, thetaOrig, thetaEst, 'false', subtitle, 'tmp.jpg');

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title(sprintf('Learned \\theta Error Convergence'));
    hold on;
    plot([errorsTrain(:).rmse], 'LineWidth', 2);
    ylabel('RMSE of learned \theta (lower is better)');
    xlabel('EM Iteration');
    ylim([0, max([errorsTrain(:).rmse], [], 2) * 1.1]);
    outpath = misc.pathMaker(params, conf.N, conf.L, 'DecError_RMSEOverIter', '.jpg');
    saveas(gcf, outpath);

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title('Log Likelihood Convergence');
    hold on;
    plot(trainLikelihood, 'LineWidth', 2);
    plot(testLikelihood, 'LineWidth', 2);
    plot([1, length(trainLikelihood)], [trainLikelihoodOrig, trainLikelihoodOrig]);
    plot([1, length(testLikelihood)], [testLikelihoodOrig, testLikelihoodOrig]);
    legend('Train likelihood with learned \theta', 'Test likelihood with learned \theta', ...
           'Train likelihood with true \theta', 'Test likelihood with true \theta', 'Location', 'southeast');
    ylabel('Log likelihood of sequences (higher is better)');
    xlabel('EM Iteration');
    xlim([1, length(trainLikelihood)]);
    outpath = misc.pathMaker(params, conf.N, conf.L, 'DecError_Likelihood', '.jpg');
    saveas(gcf, outpath);
    keyboard
end


function showTwoThetasOverTime(params, thetaOrig, thetaEsts, withExponent, outpath)
    DOT_SIZE = 20;
    % colors = ['b', 'r', 'g', 'm'];
    hold on;
    END_INTENSITY = 0.55;
    START_INTENSITY = 0.9
    minVal = inf;
    maxVal = -inf;
    for j = 1:length(thetaEsts)
        thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEsts(j));
        % m x length(theta(:)) / m
        thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
        thetaEstMat = misc.thetaToMat(params, thetaEst, true);
        % if isempty(strfind(outpath, 'tmp'))
            % fig = figure('units', 'pixels', 'Position', [0 0 1000 1000]);
        % end
        if withExponent
            thetaOrigMat = exp(thetaOrigMat);
            thetaEstMat = exp(thetaEstMat);
        end
        minVal = floor(min(minVal, min([thetaEstMat(:); thetaOrigMat(:)])));
        maxVal = ceil(max(maxVal, max([thetaEstMat(:); thetaOrigMat(:)])));
        % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
        gradColor = END_INTENSITY + (START_INTENSITY - END_INTENSITY) * (length(thetaEsts) - j) / length(thetaEsts);
        % if j > 1
        %     for t = 1: length(thetaOrigMat(:))
        %         plot([thetaOrigMat(t), thetaOrigMat(t)], [prevTheta(t), thetaEstMat(t)],...
        %              'color', [gradColor, gradColor, 1]);
        %     end
        % end
        if j == length(thetaEsts)
            gradColor = 0;
        end
        scatter(thetaOrigMat(:), thetaEstMat(:), DOT_SIZE, [gradColor, gradColor, 1], 'filled');
        prevTheta = thetaEstMat;
    end
    rmse = sqrt(mean((exp(thetaOrigMat(:)) - exp(thetaEstMat(:))) .^ 2, 1));
    plot([minVal, maxVal], [minVal, maxVal]);
    % legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \\theta vs True \\theta Comparison'));
    xlabel('True \theta probabilities');
    ylabel('Estimated \theta probabilities');
    text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), ...
         sprintf('Final EM Iteration RMSE = %.2e', rmse), 'FontSize', 14)
    saveas(gcf, outpath);
end



function errors = rateTheta(params, thetaOrig, thetaEst, dataset)
    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    YEst = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);
    [errors.totalError, errors.PWMError, errors.layerError] = calcYEstError(params, dataset, YEst);
    errors.rmse = calcThetaError(params, thetaOrig, thetaEst);
end


function rmse = calcThetaError(params, thetaOrig, thetaEst)
    thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
    thetaEstMat = misc.thetaToMat(params, thetaEst, true);
    rmse = sqrt(mean((exp(thetaOrigMat(:)) - exp(thetaEstMat(:))) .^ 2, 1));
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
