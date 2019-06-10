% PWM - k x n x J
function [PWMs, lengths, names] = PWMs()
    if ~isfile('../data/Jaspar/PWMs.mat')
        JasparDataProcessing.JasparTxtToMat()
    end
    load('../data/Jaspar/PWMs.mat');
end