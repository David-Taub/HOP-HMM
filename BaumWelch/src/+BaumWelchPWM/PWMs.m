function [PWM, lengths, names] = PWMs()
    load('data/PWMs.mat');
    PWM = PWM(1:50, :, :);
    lengths = lengths(1:50);
    names = names{1:50};
end