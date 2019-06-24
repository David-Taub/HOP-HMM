
function mainDecErrorPlot()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.maxIters = 100;
    conf.canCrossLayer = true;
    conf.patience = 10;
    conf.L = 1500;
    conf.N = 500;
    conf.withExponent = false;
    conf.repeat = 1;
    conf.order = 3;
    conf.m = 6;
    conf.k = 15;
    conf.backgroundAmount = 1;
    conf.doBound = false;
    conf.doResample = false;
    main(conf);

end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order, conf.doESharing);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    subtitle = sprintf('m=%d, k=%d, %d%% of data', conf.m, conf.k);

    [trainDataset, testDataset] = misc.crossValidationSplit(params, mergedPeaksMin, 0.15)
    % dataset.title = subtitle;
    % dataset.X = mergedPeaksMin.seqs;
    % dataset.theta = mergedPeaksMin.theta;
    % dataset.Y = mergedPeaksMin.Y;
    % dataset.Y2 = mergedPeaksMin.Y2;
    % dataset.pcPWMp = misc.preComputePWMp(mergedPeaksMin.seqs, params);

    thetaOrig = trainDataset.theta;
    thetaEst(1) = misc.genTheta(params, true);
    errorsTrain(1) = rateTheta(params, thetaOrig, thetaEst(1), trainDataset);
    errorsTest(1) = rateTheta(params, thetaOrig, thetaEst(1), testDataset);
    [~, ~, pX, ~, ~, ~] = EM.EStep(params, thetaEst(1), trainDataset.X, trainDataset.pcPWMp);
    trainLikelihood(1) = matUtils.logMatSum(pX, 1);
    [~, ~, pX, ~, ~, ~] = EM.EStep(params, thetaEst(1), testDataset.X, testDataset.pcPWMp);
    testLikelihood(1) = matUtils.logMatSum(pX, 1);
    matUtils.logMatSum(pX, 1);
    for i = 1:conf.maxIters
        [thetaEst(i + 1), trainLikelihood(i + 1)] = EM.EMIteration(params, trainDataset, thetaEst(i), conf.doBound, conf.doResample);
        thetaEst(i + 1) = misc.permThetaByAnother(params, thetaOrig, thetaEst(i + 1));
        [~, ~, pX, ~, ~, ~] = EM.EStep(params, thetaEst(i + 1), testDataset.X, testDataset.pcPWMp);
        testLikelihood(i + 1) = matUtils.logMatSum(pX, 1);

        errorsTrain(i + 1) = rateTheta(params, thetaOrig, thetaEst(i + 1) , trainDataset);
        errorsTest(i + 1) = rateTheta(params, thetaOrig, thetaEst(i + 1), testDataset);
    end
    close all
    show.showTheta(thetaEst(end));
    show.showTheta(thetaOrig);
    figure('units','normalized','outerposition',[0 0 1 1]);

    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    title('Viterbi Error')
    % plot([errorsTrain(:).totalError]);
    % plot([errorsTrain(:).PWMError]);
    plot([errorsTrain(:).layerError]);
    % plot([errorsTest(:).totalError]);
    % plot([errorsTest(:).PWMError]);
    plot([errorsTest(:).layerError]);
    % legend('Total train', 'PWM train', 'Layer train', 'Total test', 'PWM test', 'Layer test')
    legend('Layer train error', 'Layer test error');
    outpath = sprintf('DecError1_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    saveas(gcf, outpath);

    figure('units','normalized','outerposition',[0 0 1 1]);
    outpath = sprintf('DecError2_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    showTwoThetasOverTime(params, thetaOrig, thetaEst, 'false', subtitle, outpath);
    % show.showTwoThetas(params, thetaOrig, thetaEst, 'false', subtitle, 'tmp.jpg');

    figure('units','normalized','outerposition',[0 0 1 1]);
    outpath = sprintf('DecError3_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    title('MSLE')
    hold on;
    plot([errorsTrain(:).msle])
    saveas(gcf, outpath);

    outpath = sprintf('DecError4_m%dk%do%db%dr%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, conf.doResample, conf.N, conf.L);
    figure('units','normalized','outerposition',[0 0 1 1]);
    title('Log Likelihood')
    hold on;
    plot(trainLikelihood)
    plot(testLikelihood)
    legend('train likelihood', 'test likelihood')
    saveas(gcf, outpath);
    keyboard
end


function showTwoThetasOverTime(params, thetaOrig, thetaEsts, withExponent, subtitle, outpath)
    DOT_SIZE = 20;
    % colors = ['b', 'r', 'g', 'm'];
    hold on;
    BASE_INTENSITY = 0.8;
    minVal = inf;
    maxVal = -inf;
    for j = 1:length(thetaEsts)
        thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEsts(j));
        % m x length(theta(:)) / m
        thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
        thetaEstMat = misc.thetaToMat(params, thetaEst, true);
        % if isempty(strfind(outpath, 'tmp'))
            % fig = figure('units','normalized','outerposition',[0 0 1 1]);
        % end
        if withExponent
            thetaOrigMat = exp(thetaOrigMat);
            thetaEstMat = exp(thetaEstMat);
        end
        minVal = floor(min(minVal, min([thetaEstMat(:); thetaOrigMat(:)])));
        maxVal = ceil(max(maxVal, max([thetaEstMat(:); thetaOrigMat(:)])));
        % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
        gradColor = BASE_INTENSITY + (1 - BASE_INTENSITY) * (length(thetaEsts) - j) / length(thetaEsts);
        if j > 1
            for t = 1: length(thetaOrigMat(:))
                plot([thetaOrigMat(t), thetaOrigMat(t)], [prevTheta(t), thetaEstMat(t)], 'color', [gradColor, gradColor, 1]);
            end
        end
        if j == length(thetaEsts)
            gradColor = 0;
        end
        scatter(thetaOrigMat(:), thetaEstMat(:), DOT_SIZE, [gradColor, gradColor, 1], 'filled');
        prevTheta = thetaEstMat;
    end
    msle = mean((log(exp(thetaOrigMat(:)) + 1) - log(exp(thetaEstMat(:)) + 1)) .^ 2, 1);
    plot([minVal, maxVal], [minVal, maxVal]);
    % legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \\theta vs True \\theta (%s)', subtitle));
    xlabel('True \theta');
    ylabel('Estimated \theta');
    text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), sprintf('MSLE=%.2e', msle), 'FontSize', 14)
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
