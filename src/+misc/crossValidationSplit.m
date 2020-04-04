
% testTrainRatio - bigger means more more test
function [test, train] = crossValidationSplit(params, mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    X = mergedPeaksMin.seqs; % N x L x [1|3]
    % N x k x L
    pcPWMp = misc.preComputePWMp(X, params);
    N = size(X, 1);
    trainMask = true(N, 1);
    trainMask(1: floor(N * testTrainRatio)) = false;
    train.title = 'Train';
    train.X = X(trainMask, :);
    test.title = 'Test';
    test.X = X(~trainMask, :);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
    if isfield(mergedPeaksMin, 'starts')
        % real data train
        train.starts = mergedPeaksMin.starts(trainMask);
        train.chrs = mergedPeaksMin.chrs(trainMask);
        train.samplesCount = mergedPeaksMin.samplesCount(trainMask);
        % N x 1
        [~, train.Y] = max(mergedPeaksMin.overlaps(trainMask, :), [], 2);

        % train.Y = (train.Y == mergedPeaksMin.backgroundInd) + 1;

        train.backgroundInd = mergedPeaksMin.backgroundInd;
        train.tissueNames = mergedPeaksMin.tissueNames;
        train.tissueList = mergedPeaksMin.tissueList;
        train.tissueEIDs = mergedPeaksMin.tissueEIDs;

        % real data test
        test.starts = mergedPeaksMin.starts(~trainMask);
        test.chrs = mergedPeaksMin.chrs(~trainMask);
        test.samplesCount = mergedPeaksMin.samplesCount(~trainMask);
        [~, test.Y] = max(mergedPeaksMin.overlaps(~trainMask, :), [], 2);

        % test.Y = (test.Y == mergedPeaksMin.backgroundInd) + 1;

        test.backgroundInd = mergedPeaksMin.backgroundInd;
        test.tissueNames = mergedPeaksMin.tissueNames;
        test.tissueList = mergedPeaksMin.tissueList;
        test.tissueEIDs = mergedPeaksMin.tissueEIDs;
    end
    if isfield(mergedPeaksMin, 'Y')
        train.Y = mergedPeaksMin.Y(trainMask, :);
        test.Y = mergedPeaksMin.Y(~trainMask, :);
    end
    if isfield(mergedPeaksMin, 'Y2')
        train.Y2 = mergedPeaksMin.Y2(trainMask, :);
        test.Y2 = mergedPeaksMin.Y2(~trainMask, :);
    end
    if isfield(mergedPeaksMin, 'theta')
        train.theta = mergedPeaksMin.theta;
        test.theta = mergedPeaksMin.theta;
    end
end