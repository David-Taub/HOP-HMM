function [X, Y] = mainGenSequences()
    params.m = 5;
    params.order = 3;
    params.N = 1000;
    params.L = 500;
    params.tEpsilon = 0.01;
    [PWMs, ~] = BaumWelchPWM.PWMs();
    [params.k, params.n, params.J] = size(PWMs);

    [theta] = BaumWelchPWM.genThetaJ(params);
    [X, Y] = BaumWelchPWM.genSequencesJ(theta, params);

    clear PWMs;
    save(fullfile('data', 'dummyDNA.mat'));
end