% TODO: BPI bit per instance
% TODO: add tutorial to the post graph
% TODO: add note about areas where the posterior confused between different Enhancers)
% TODO: single permutation
function mainDecErrorPlot()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.startTUniform = false;
    conf.withExponent = true;
    conf.doResampling = false;
    conf.background_g_noise = 0.050;
    conf.m = 5;
    conf.k = 25;
    conf.backgroundAmount = 1;

    conf.L = 1000;
    conf.maxIters = 40;
    conf.patience = conf.maxIters;
    conf.repeat = 10;
    conf.order = 3;
    conf.N = 1000;


    conf.doGTBound = true;
    main(conf);
    conf.doGTBound = false;
    main(conf);
end


function [likelihood, rmse] = main(conf)
    plotErrors = true;
    ind = 1;
    j = 1;
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing, conf.doGTBound);
    mergedPeaksMin = misc.genSyntheticMergedPeaksMin(conf.N, conf.L, params, conf.startWithBackground, conf.background_g_noise);
    thetaOrig = mergedPeaksMin.theta;
    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15);

    for j = 1:conf.repeat
        [~, ~, thetaEsts] = EM.EM(trainDataset, params, conf.maxIters, conf.patience, 1)
        distMat = pdist2(exp(misc.thetaToMat(params, thetaOrig, false)), ...
                         exp(misc.thetaToMat(params, thetaEsts(end), false)))
        perm = misc.munkres(distMat)';
        for i = 1:conf.maxIters
            thetaEsts(i) = misc.permTheta(thetaEsts(i), perm);
        end

        show.showTheta(thetaOrig);
        show.showTheta(thetaEsts(end));

        close all;
        % raindrop graph
        outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, sprintf('DecError_ThetaOverIter_%d', j), '.jpg');
        showTwoThetasOverTime(conf, params, thetaOrig, thetaEsts, conf.withExponent, outpath);

        if not(plotErrors)
            continue
        end

        for i = 1:conf.maxIters
            testLikelihood(i) = calcLikelihood(params, thetaEsts(i), testDataset.X, testDataset.pcPWMp);
            trainLikelihood(i) = calcLikelihood(params, thetaEsts(i), trainDataset.X, trainDataset.pcPWMp);
            errorsTrain(i) = rateTheta(params, thetaOrig, thetaEsts(i) , trainDataset);
            errorsTest(i) = rateTheta(params, thetaOrig, thetaEsts(i), testDataset);
        end
        runsRepetitions(j).errorsTrain = errorsTrain;
        runsRepetitions(j).trainLikelihood = trainLikelihood;
        runsRepetitions(j).errorsTest = errorsTest;
        runsRepetitions(j).testLikelihood = testLikelihood;
    end


    if not(plotErrors)
        return
    end
    % original theta results
    origTrainError = rateTheta(params, thetaOrig, thetaOrig, trainDataset);
    origTestError = rateTheta(params, thetaOrig, thetaOrig, testDataset);
    origTrainError.likelihood = calcLikelihood(params, thetaOrig, trainDataset.X, trainDataset.pcPWMp);
    origTestError.likelihood = calcLikelihood(params, thetaOrig, testDataset.X, testDataset.pcPWMp);


    viterbiErrorConvergence(params, conf, origTrainError, origTestError, runsRepetitions);
    thetaRMSEConvergence(params, conf, runsRepetitions);
    likelihoodConvergence(params, conf, runsRepetitions, origTrainError, origTestError);
    likelihoodRMSEScatter(params, conf, runsRepetitions);
    likelihoodRMSEScatterConvergence(params, conf, runsRepetitions);
end


function errors = rateTheta(params, thetaOrig, thetaEst, dataset)
    YEst = misc.viterbi(params, thetaEst, dataset.X, dataset.pcPWMp);
    [errors.totalError, errors.PWMError, errors.layerError] = calcYEstError(params, dataset, YEst);
    errors.rmse = calcThetaError(params, thetaOrig, thetaEst);
end


function rmse = calcThetaError(params, thetaOrig, thetaEst)
    thetaOrigMat = misc.thetaToMat(params, thetaOrig, false);
    thetaEstMat = misc.thetaToMat(params, thetaEst, false);
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


function likelihood = calcLikelihood(params, theta, X, pcPWMp)
    N = size(X, 1);
    [~, ~, pX, ~, ~, ~] = EM.EStep(params, theta, X, pcPWMp);
    likelihood = mean(pX, 1);
end


function viterbiErrorConvergence(params, conf, origTrainError, origTestError, runsRepetitions)
    ORIGINAL_TRAIN_PLOT_COLOR =[0.9290, 0.6940, 0.1250];
    ORIGINAL_TEST_PLOT_COLOR = [0.4940, 0.1840, 0.5560];
    TRAIN_PLOT_COLOR = [0, 0.4470, 0.7410];
    TEST_PLOT_COLOR = [0.8500, 0.3250, 0.0980];
    PLOT_TRANSPERENCY = 0.5;
    % error decreasing plot
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title('Viterbi Classification Convergence');
    hold on;

    % original theta error train
    origTrainPlot = plot([1, conf.maxIters], [origTrainError.layerError, origTrainError.layerError] .* 100, ...
                        'LineWidth', 2);
    origTrainPlot.Color = ORIGINAL_TRAIN_PLOT_COLOR;

    % original theta error test
    origTestPlot = plot([1, conf.maxIters], [origTestError.layerError, origTestError.layerError] .* 100, ...
                        'LineWidth', 2);
    origTestPlot.Color = ORIGINAL_TEST_PLOT_COLOR;

    meanTrainError = zeros(1, conf.maxIters);
    meanTestError = zeros(1, conf.maxIters);

    for j = 1:length(runsRepetitions)
        meanTrainError = meanTrainError + ([runsRepetitions(j).errorsTrain(:).layerError] ./ length(runsRepetitions))
        meanTestError = meanTestError + ([runsRepetitions(j).errorsTest(:).layerError] ./ length(runsRepetitions))
        trainPlot = plot([runsRepetitions(j).errorsTrain(:).layerError] .* 100);
        trainPlot.Color = [TRAIN_PLOT_COLOR, PLOT_TRANSPERENCY];
        testPlot = plot([runsRepetitions(j).errorsTest(:).layerError] .* 100);
        testPlot.Color = [TEST_PLOT_COLOR, PLOT_TRANSPERENCY];
    end

    % plot mean error graph
    meanTrainPlot = plot(meanTrainError .* 100, 'LineWidth', 2);
    meanTrainPlot.Color = TRAIN_PLOT_COLOR;
    meanTestPlot = plot(meanTestError .* 100, 'LineWidth', 2);
    meanTestPlot.Color = TEST_PLOT_COLOR;

    hold off;

    ylabel('Misclassification Rate (lower is better)');
    ytickformat('percentage');
    xlabel('EM Iteration');
    legend('Viterbi train error with learned \theta', 'Viterbi test error with learned \theta', ...
           'Viterbi train error with true \theta', 'Viterbi test error with true \theta');
    xlim([1, conf.maxIters]);
    ylim([0,100]);
    outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, 'DecError_VitErr', '.jpg');
    saveas(gcf, outpath);
end


function thetaRMSEConvergence(params, conf, runsRepetitions)
    TRAIN_PLOT_COLOR = [0, 0.4470, 0.7410];
    TEST_PLOT_COLOR = [0.8500, 0.3250, 0.0980];
    PLOT_TRANSPERENCY = 0.5;

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title(sprintf('Learned \\theta Error Convergence'));
    hold on;
    meanTrainError = zeros(1, length([runsRepetitions(1).errorsTrain(:).rmse]));
    for j = 1:length(runsRepetitions)
        % plot([errorsTrain(:).rmse], 'LineWidth', 2);
        meanTrainError = meanTrainError + ([runsRepetitions(j).errorsTrain(:).rmse] ./ length(runsRepetitions))
        trainPlot = plot([runsRepetitions(j).errorsTrain(:).rmse]);
        trainPlot.Color = [TRAIN_PLOT_COLOR, PLOT_TRANSPERENCY];
    end
    meanTrainPlot = plot(meanTrainError, 'LineWidth', 2);
    meanTrainPlot.Color = TRAIN_PLOT_COLOR;

    hold off;
    ylabel('RMSE of learned \theta (lower is better)');
    xlabel('EM Iteration');
    ylim([0, 0.25]);
    outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, 'DecError_RMSEOverIter', '.jpg');
    saveas(gcf, outpath);
end


function likelihoodConvergence(params, conf, runsRepetitions, origTrainError, origTestError)
    ORIGINAL_TRAIN_PLOT_COLOR =[0.9290, 0.6940, 0.1250];
    ORIGINAL_TEST_PLOT_COLOR = [0.4940, 0.1840, 0.5560];
    TRAIN_PLOT_COLOR = [0, 0.4470, 0.7410];
    TEST_PLOT_COLOR = [0.8500, 0.3250, 0.0980];
    PLOT_TRANSPERENCY = 0.5;

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title('Log Likelihood Convergence');
    hold on;
    origTrainPlot = plot([1, conf.maxIters], [origTrainError.likelihood, origTrainError.likelihood], 'LineWidth', 2);
    origTestPlot = plot([1, conf.maxIters], [origTestError.likelihood, origTestError.likelihood], 'LineWidth', 2);
    origTrainPlot.Color = ORIGINAL_TRAIN_PLOT_COLOR;
    origTestPlot.Color = ORIGINAL_TEST_PLOT_COLOR;
    meanTrain = zeros(1, length(runsRepetitions(1).trainLikelihood));
    meanTest = zeros(1, length(runsRepetitions(1).testLikelihood));
    for j = 1:length(runsRepetitions)
        meanTrain = meanTrain + (runsRepetitions(j).trainLikelihood ./ length(runsRepetitions))
        meanTest = meanTest + (runsRepetitions(j).testLikelihood ./ length(runsRepetitions))
        trainPlot = plot(runsRepetitions(j).trainLikelihood);
        testPlot = plot(runsRepetitions(j).testLikelihood);
        trainPlot.Color = [TRAIN_PLOT_COLOR, PLOT_TRANSPERENCY];
        testPlot.Color = [TEST_PLOT_COLOR, PLOT_TRANSPERENCY];
    end
    meanTrainPlot = plot(meanTrain, 'LineWidth', 2);
    meanTestPlot = plot(meanTest, 'LineWidth', 2);
    meanTrainPlot.Color = TRAIN_PLOT_COLOR;
    meanTestPlot.Color = TEST_PLOT_COLOR;
    hold off;
    hold on;
    legend('Train likelihood with learned \theta', 'Test likelihood with learned \theta', ...
           'Train likelihood with true \theta', 'Test likelihood with true \theta', ...
           'Location', 'southeast');
    ylabel('Log likelihood of sequences (higher is better)');
    xlabel('EM Iteration');
    xlim([1, conf.maxIters]);
    ylim([-conf.L * 2, -conf.L]);
    outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, 'DecError_Likelihood', '.jpg');
    saveas(gcf, outpath);
end


function likelihoodRMSEScatter(params, conf, runsRepetitions)
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    hold on;
    title('likelihood vs RMSE');
    for j = 1:length(runsRepetitions)
        likelihood(j) = runsRepetitions(j).trainLikelihood(end);
        rmse(j) = runsRepetitions(j).errorsTrain(end).rmse;
    end
    scatter(likelihood, rmse, 20, 'filled');
    ylabel('Error (RMSE)');
    xlabel('Train Log Likelihood');
    outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, ...
                             'DecError_rmse_vs_likelihood', '.jpg');
    saveas(gcf, outpath);
end


function likelihoodRMSEScatterConvergence(params, conf, runsRepetitions)
    DOT_SIZE = 20;
    END_INTENSITY = 0.55;
    START_INTENSITY = 0.9;

    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    hold on;
    title('likelihood vs RMSE');
    for i = 1:conf.maxIters
        for j = 1:length(runsRepetitions)
            likelihood(i, j) = runsRepetitions(j).trainLikelihood(i);
            rmse(i, j) = runsRepetitions(j).errorsTrain(i).rmse;
        end
        gradColor = END_INTENSITY + (START_INTENSITY - END_INTENSITY) * (conf.maxIters - i) / conf.maxIters;
        if i == conf.maxIters
            gradColor = 0;
        end
        scatter(likelihood(i, :), rmse(i, :), DOT_SIZE, [gradColor, gradColor, 1], 'filled');
    end
    plot(likelihood(i, :), rmse(i, :), 'Color', [START_INTENSITY, START_INTENSITY, 1]);
    yMax = max(rmse(:)) * 1.1;
    ylim([0, yMax]);
    ylabel('Error (RMSE)');
    xlabel('Train Log Likelihood');
    outpath = misc.pathMaker(params, conf.N, conf.L, conf.background_g_noise, ...
                             'DecError_rmse_vs_likelihood_convergence', '.jpg');
    saveas(gcf, outpath);
    keyboard;
end



function showTwoThetasOverTime(conf, params, thetaOrig, thetaEsts, withExponent, outpath)
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    DOT_SIZE = 20;
    END_INTENSITY = 0.55;
    START_INTENSITY = 0.9;
    hold on;
    minVal = inf;
    maxVal = -inf;
    for i = 1:conf.maxIters
        thetaEst = thetaEsts(i);
        % m x length(theta(:)) / m
        thetaOrigMat = misc.thetaToMat(params, thetaOrig, false);
        thetaEstMat = misc.thetaToMat(params, thetaEst, false);
        % if isempty(strfind(outpath, 'tmp'))
            % fig = figure('units', 'pixels', 'Position', [0 0 1000 1000]);
        % end
        if withExponent
            thetaOrigMat = exp(thetaOrigMat);
            thetaEstMat = exp(thetaEstMat);
        else
            minVal = floor(min(minVal, min([thetaEstMat(:); thetaOrigMat(:)])));
            maxVal = ceil(max(maxVal, max([thetaEstMat(:); thetaOrigMat(:)])));
        end
        % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
        gradColor = END_INTENSITY + (START_INTENSITY - END_INTENSITY) * (conf.maxIters - i) / conf.maxIters;
        if i == conf.maxIters
            gradColor = 0;
        end
        scatter(thetaOrigMat(:), thetaEstMat(:), DOT_SIZE, [gradColor, gradColor, 1], 'filled');
    end
    if withExponent
        minVal = 0;
        maxVal = 1;
    end
    rmse = sqrt(mean((exp(thetaOrigMat(:)) - exp(thetaEstMat(:))) .^ 2, 1));
    plot([minVal, maxVal], [minVal, maxVal]);
    title(sprintf('Learned \\theta vs True \\theta Comparison'));
    xlabel('True \theta probabilities');
    ylabel('Estimated \theta probabilities');
    % text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), ...
    %      sprintf('Final EM Iteration RMSE = %.2e', rmse), 'FontSize', 14);
    saveas(gcf, outpath);
end