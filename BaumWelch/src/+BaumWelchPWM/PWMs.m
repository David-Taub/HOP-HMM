% k x n x J
function [PWM, lengths, names] = PWMs()
	S = load('data/PWMs.mat');
	PWM = S.PWMs;
	names = S.names;
	% n x J x k -> k x J x n
	PWM = permute(PWM, [3, 2, 1]);
	lengths = sum(sum(PWM, 3) > 0, 2)';
	emptyMap = repmat(sum(PWM, 3) == 0, [1,1,4]);
	PWM = PWM + 1;
	PWM = bsxfun(@times, PWM, 1 ./ sum(PWM, 3));
	PWM(emptyMap) = 0;
	% k x J x N -> k x n x J
	PWM = permute(PWM, [1, 3, 2]);
end