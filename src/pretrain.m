function [E, G] = pretrain(params, dataset)
    dbstop if error
    paramsCopy = params;
    paramsCopy.m = 1;
    paramsCopy.doGTBound = 0;
    paramsCopy.backgroundAmount = 0;
    paramsCopy.enhancerAmount = 0;

    L = size(dataset.X, 2);
    NEW_L = 100;
    maxIter = 10;
    patience = 3;
    repeat = 2;

    E = zeros(params.m, params.n ^ params.order);
    G = zeros(params.m, params.k);

    for i = 1:params.m
        mask = dataset.Y == i;
        subDataset.X = dataset.X(mask, round(L / 2 - NEW_L / 2): round(L / 2 + NEW_L / 2));
        subDataset.pcPWMp = dataset.pcPWMp(mask, :, round(L / 2 - NEW_L / 2): round(L / 2 + NEW_L / 2));

        [theta, ~, ~] = EM.EM(dataset, params, maxIter, patience, repeat);

        E(i, :) = theta.E(1, :);
        G(i, :) = theta.G(1, :);
    end
end