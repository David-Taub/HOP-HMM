function josh(dataset)
    VAR_THRESHOLD = 0.04;

    l = length(dataset.tissueType);
    [N, ~] = size(dataset.samples{1});

    varPeopleMask = false(N, m);
    dataset.samplesAvg = zeros(N, l*2);
    for i = 1: 2 * l
        varPeopleMask(:, i) = var(dataset.samples{i},0,2) > VAR_THRESHOLD;
        dataset.samplesAvg(:, i) = mean(dataset.samples{i},2);
    end

end

