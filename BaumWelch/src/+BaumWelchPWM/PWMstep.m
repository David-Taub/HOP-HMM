% slice - N x m x k
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% PWMsRep - N x J x n x k: N replication of emission matrix of m Jaspar PWMs with length J 
%        true length of i'th PWM< J and given in lengths(i) if a PWM is 
%        shorter than j, it is aligned to the end of the 3rd dimension.
% res - N x m
function res = PWMstep(slice, PWMsRep, Xs1H, Y, y, t)
    k = size(Y, 2);

    % N x m x k -> m x k x N
    sliceRep = permute(slice, [2,3,1]);
    % probability to get into the pwm submode
    % N x k x m
    R = permute(bsxfun(@times, sliceRep, Y), [3,2,1]);
    % N x k
    PWMProb = getPWMp(PWMstep, Xs1H, t);
    % N x k x m, N x k
    res = sumDim(bsxfun(@times, R, PWMProb), 2);
end