function [test, train] = crossValidationSplit(params, mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    X = mergedPeaksMin.seqs; % N x L x [1|3]
    % N x k x L
    pcPWMp = misc.preComputePWMp(X, params);
    N = size(X, 1);
    trainMask = true(N, 1);
    trainMask(randperm(N, floor(N * testTrainRatio))) = false;
    train.title = 'Train';
    train.X = X(trainMask, :);
    test.title = 'Test';
    test.X = X(~trainMask, :);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
    if isfield(mergedPeaksMin, 'starts')
        train.starts = mergedPeaksMin.starts(trainMask);
        train.chrs = mergedPeaksMin.chrs(trainMask);
        test.starts = mergedPeaksMin.starts(~trainMask);
        test.chrs = mergedPeaksMin.chrs(~trainMask);
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