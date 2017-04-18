function [X, Y] = mainGenSequences()
    m = 5;
    order = 3;
    L = 500; N = 1000;
    n = 4;
    [PWMs, lengths] = BaumWelchPWM.PWMs();
    k = size(PWMs, 1);
    [startT, T, E, M, F] = BaumWelchPWM.genRandParamsJ(m, n, order, k);
    [X, Y] = BaumWelchPWM.genSequencesJ(startT, T, E, M, F, L, N, PWMs, lengths);

    clear PWMs lengths;
    save(fullfile('data', 'dummyDNA.mat'));
end