% PWM - k x n x J
% the JASPAR dataset has k=519 PWMs with max length of J=20, median length of 10
function [PWMs, lengths, names] = PWMs()
    if ~isfile('../data/Jaspar/PWMs.mat')
        JasparDataProcessing.JasparTxtToMat();
    end
    load('../data/Jaspar/PWMs.mat');
end